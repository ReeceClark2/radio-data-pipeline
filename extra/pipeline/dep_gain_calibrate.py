# Third-party libraries
import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.table import Column
from astropy.time import Time
from scipy.stats import linregress
import rcr
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


class Gain_Calibrate:
    def __init__(self, filepath):
        self.filepath = filepath
        self.indices = None
        self.data_subsets = None

        self.gain_start = None
        self.gain_end = None

        try:
            with fits.open(filepath) as hdul:
                self.header = hdul[0].header
                self.data = Table(hdul[1].data)
            
            # TODO add logging here
        except Exception as e:
            # TODO add logging here
            pass


    def find_calibrations(self):
        ifnums = np.unique(self.data['IFNUM'])
        plnums = np.unique(self.data['PLNUM'])

        channel_count = len(ifnums) * len(plnums)

        data_start_ind = None
        post_cal_start_ind = None

        counter = 0

        cal_started = False
        pre_cal_complete = False

        for ind, i in enumerate(self.data):
            if i['CALSTATE'] == 1:
                cal_started = True
            
            if cal_started and i["CALSTATE"] == 0 and i["SWPVALID"] == 1 and not pre_cal_complete:
                data_start_ind = ind
                pre_cal_complete = True

            if ind > 0 and pre_cal_complete and i["SWPVALID"] == 0 and self.data[ind - 1]["SWPVALID"] == 0:
                if post_cal_start_ind is None:
                    post_cal_start_ind = ind - 1
            else:
                post_cal_start_ind = None

            if pre_cal_complete and i['CALSTATE'] == 0 and i['SWPVALID'] == 1:
                counter += 1

            if counter <= 3 * channel_count and i['SWPVALID'] == 0 and data_start_ind:
                print('counter', data_start_ind)
                data_start_ind = None
                pre_cal_complete = False

            if pre_cal_complete and i['SWPVALID'] == 0 and i['CALSTATE'] == 1:
                break


        self.indices = [data_start_ind, post_cal_start_ind]

        if self.header['OBSMODE'] == 'onoff':
            for ind, i in enumerate(self.data):
                target = 'onoff:off'

                if target in i['OBSMODE']:
                    offstart = ind
                    self.indices = [data_start_ind, offstart-1, offstart, post_cal_start_ind]

                    break


    def linear(self, x, params): # model function
        return params[0] + x * params[1]


    def d_linear_1(self, x, params): # first model parameter derivative
        return 1


    def d_linear_2(self, x, params): # second model parameter derivative
        return x


    def rcr(self, array):
        '''
        Perform functional Robust Chauvenet Rejection on a given dataset.
        '''

        x = array[0]
        x -= np.average(x)

        y = array[1]

        result = linregress(x, y)
        m = result.slope
        b = result.intercept
        guess = [m, b]
        model = rcr.FunctionalForm(self.linear,
            x,
            y,
            [self.d_linear_1, self.d_linear_2],
            guess
        )
        
        r = rcr.RCR(rcr.SS_MEDIAN_DL) 
        r.setParametricModel(model)
        r.performBulkRejection(y)

        indices = r.result.indices

        x = np.array([x[i] for i in indices])
        y = np.array([y[i] for i in indices])

        best_fit_parameters = model.result.parameters
        
        sigma = (1 / (len(x) - 2)) * np.sum((y - best_fit_parameters[1] * x - best_fit_parameters[0]) ** 2)
        m_sd = np.sqrt(sigma / np.sum((x - np.mean(x)) ** 2))
        b_sd = np.sqrt(sigma * ((1 / len(x)) + ((np.mean(x) ** 2) / np.sum((x - np.mean(x)) ** 2))))
        uncertainties = (b_sd, m_sd)

        return best_fit_parameters, uncertainties


    def compute_gain_deltas(self):
        def average_neighbors(values, index, window=0):
            if window == 0:
                return values[index] if values[index] is not None else None
            neighbors = [
                v for v in values[max(0, index-window):index+window+1]
                if v is not None
            ]
            return np.mean(neighbors) if neighbors else None

        ifnums = np.unique(self.data['IFNUM'])
        plnums = np.unique(self.data['PLNUM'])

        first_len = len(self.data[0]['DATA'])
        if not all(len(row['DATA']) == first_len for row in self.data):
            return

        t0 = Time(self.header['DATE'], format='isot')

        # Cache: (idx, row) pairs + parsed time
        cached_data = [
            {
                'idx': idx,
                'time': Time(row['DATE-OBS'], format='isot'),
                'SWPVALID': row['SWPVALID'],
                'CALSTATE': row['CALSTATE'],
                'IFNUM': row['IFNUM'],
                'PLNUM': row['PLNUM'],
                'DATA': row['DATA']
            }
            for idx, row in enumerate(self.data)
        ]

        # Partition into pre and post ranges
        pre_data = [row for row in cached_data if row['idx'] < self.indices[0]]

        data_array = [row for row in cached_data if self.indices[0] < row['idx'] < self.indices[1]]
        post_data = [row for row in cached_data if row['idx'] > self.indices[1]]

        def extract(data_subset, swpvalid_val, calstate_val, i, j):
            rows = [
                ((row['time'] - t0).sec, np.mean(row['DATA']))
                for row in data_subset
                if row['SWPVALID'] == swpvalid_val and
                row['CALSTATE'] == calstate_val and
                row['IFNUM'] == i and
                row['PLNUM'] == j
            ]
            if not rows:
                return np.array([]), np.array([])
            times, intensities = zip(*rows)
            return np.array(times), np.array(intensities, dtype=float)

        pre_cal_deltas = []
        post_cal_deltas = []
        pre_cal_uncs = []
        post_cal_uncs = []

        for i in ifnums:
            for j in plnums:
                data = extract(data_array, 1, 0, i, j)

                pre_cal_delta = None
                pre_cal_unc = None
                try:
                    # Extract pre arrays
                    pre_on_times, pre_on_vals = extract(pre_data, 0, 1, i, j)
                    pre_off_times, pre_off_vals = extract(pre_data, 0, 0, i, j)

                    if len(pre_on_times) == 0 or len(pre_off_times) == 0:
                        pre_cal_deltas.append(None)
                        pre_cal_uncs.append(None)
                    
                    pre_on_params, pre_on_unc = self.rcr([pre_on_times, pre_on_vals])
                    pre_off_params, pre_off_unc = self.rcr([pre_off_times, pre_off_vals])
                    # Midpoint time
                    pre_cal_time = np.mean([np.mean(pre_on_times), np.mean(pre_off_times)])
                    pre_cal_time = np.clip(pre_cal_time,
                                (pre_on_times[0] + pre_off_times[-1]) / 2,
                                (pre_on_times[-1] + pre_off_times[0]) / 2)

                    pre_cal_delta = np.abs(
                        (pre_on_params[1] * (pre_cal_time - np.mean(pre_on_times)) + pre_on_params[0]) -
                        (pre_off_params[1] * (pre_cal_time - np.mean(pre_off_times)) + pre_off_params[0])
                    )
                    pre_cal_unc = np.sqrt(
                        pre_on_unc[0]**2 + pre_off_unc[0]**2 +
                        (pre_on_unc[1] * pre_cal_time)**2 + (pre_off_unc[1] * pre_cal_time)**2
                    )

                    pre_cal_deltas.append(pre_cal_delta)
                    pre_cal_uncs.append(pre_cal_unc)
                except:
                    pre_cal_deltas.append(None)
                    pre_cal_uncs.append(None)


                post_cal_delta = None
                post_cal_unc = None
                try:
                    # Extract post arrays
                    post_on_times, post_on_vals = extract(post_data, 0, 1, i, j)
                    post_off_times, post_off_vals = extract(post_data, 0, 0, i, j)
                    
                    if len(post_on_times) == 0 or len(post_off_times) == 0:
                        post_cal_deltas.append(None)
                        post_cal_uncs.append(None)

                    post_on_params, post_on_unc = self.rcr((post_on_times, post_on_vals))
                    post_off_params, post_off_unc = self.rcr((post_off_times, post_off_vals))

                    post_cal_time = np.mean([np.mean(post_on_times), np.mean(post_off_times)])
                    post_cal_time = np.clip(post_cal_time,
                                (post_on_times[0] + post_off_times[-1]) / 2,
                                (post_on_times[-1] + post_off_times[0]) / 2)

                    post_cal_delta = np.abs(
                        (post_on_params[1] * (post_cal_time - np.mean(post_on_times)) + post_on_params[0]) -
                        (post_off_params[1] * (post_cal_time - np.mean(post_off_times)) + post_off_params[0])
                    )
                    post_cal_unc = np.sqrt(
                        post_on_unc[0]**2 + post_off_unc[0]**2 +
                        (post_on_unc[1] * post_cal_time)**2 + (post_off_unc[1] * post_cal_time)**2
                    )

                    post_cal_deltas.append(post_cal_delta)
                    post_cal_uncs.append(post_cal_unc)
                except:
                    post_cal_deltas.append(None)
                    post_cal_uncs.append(None)

            
                if pre_cal_delta and post_cal_delta:
                    print("Using pre and post calibration!")
                    # pre_cal_delta = average_neighbors(pre_cal_deltas, h, 5)
                    # post_cal_delta = average_neighbors(post_cal_deltas, h, 5)

                    z_value = abs(pre_cal_delta - post_cal_delta) / np.sqrt(pre_cal_unc**2 + post_cal_unc**2)

                    if z_value < 0.6745:
                        weights =  np.array([1 / pre_cal_unc, 1 / post_cal_unc])
                        delta = np.average([pre_cal_delta, post_cal_delta], weights=weights)

                        self.data['DATA'][:] /= delta
                        continue

                    else:
                        for ind, _ in enumerate(data):
                            delta = pre_cal_delta + (post_cal_delta - pre_cal_delta) * (data[0][ind] - pre_cal_time) / (post_cal_time - pre_cal_time)
                            
                            self.data['DATA'][ind] /= delta
                            continue

                elif pre_cal_delta and not post_cal_delta:
                    print("Using pre calibration only!")
                    
                    self.data['DATA'][:] /= pre_cal_delta
                    continue

                elif not pre_cal_delta and post_cal_delta:
                    print("Using post calibration only!")
                    
                    self.data['DATA'][:] /= post_cal_delta
                    continue

                print("No calibration found!")

        # print(np.median(pre_cal_deltas), np.median(post_cal_deltas))
        mask = (
            (self.data['IFNUM'] == 0) &
            (self.data['PLNUM'] == 0)
        )
        subset_data = self.data[mask]

        continuum = []
        calstate_vals = []
        swpvalid_vals = []

        for row in subset_data:
            continuum.append(np.mean(row['DATA']))
            calstate_vals.append(row['CALSTATE'])
            swpvalid_vals.append(row['SWPVALID'])

        # X-axis values
        x_vals = list(range(len(continuum)))

        # Create figure and primary y-axis
        fig, ax1 = plt.subplots(figsize=(10, 5))

        # Plot continuum
        ax1.plot(x_vals, continuum, color='black', label='Continuum')
        ax1.set_xlabel('Index')
        ax1.set_ylabel('Mean Intensity', color='black')
        ax1.tick_params(axis='y', labelcolor='black')
        ax1.set_xlim(0, 150)
        ax1.set_ylim(min(continuum), max(continuum))

        # Mark indices
        ax1.axvline(x=self.indices[0] / 2, color='green', linestyle='--', label='Start Index')
        ax1.axvline(x=self.indices[1] / 2, color='red', linestyle='--', label='Stop Index')

        # Twin y-axis for CALSTATE and SWPVALID
        ax2 = ax1.twinx()
        ax2.step(x_vals, calstate_vals, where='mid', color='blue', label='CALSTATE')
        ax2.step(x_vals, swpvalid_vals, where='mid', color='orange', label='SWPVALID', linestyle='--')
        ax2.set_ylabel('CALSTATE / SWPVALID', color='gray')
        ax2.tick_params(axis='y', labelcolor='gray')
        ax2.set_ylim(-0.1, 2.1)  # Assuming values 0 or 1

        # Combine legends
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

        plt.tight_layout()
        plt.savefig('Testing.png', dpi=400)
        plt.close()





    def save(self, output_path=None):
        if output_path is None:
            base, ext = os.path.splitext(self.filepath)
            output_path = f"{base}_gain_calibrated{ext}"

        # Create Primary HDU with updated header
        primary_hdu = fits.PrimaryHDU(header=self.header)

        # Convert Table back to FITS HDU
        table_hdu = fits.BinTableHDU(data=self.data)

        # Combine HDUs
        hdulist = fits.HDUList([primary_hdu, table_hdu])

        # Write to file
        hdulist.writeto(output_path, overwrite=True)
               

    def gain_calibrate(self):
        self.find_calibrations()
        self.compute_gain_deltas()
        
        self.save()



if __name__ == "__main__":
    filepath = "C:/Users/starb/Downloads/0117683.fits"
    file = Gain_Calibrate(filepath)

    file.gain_calibrate()