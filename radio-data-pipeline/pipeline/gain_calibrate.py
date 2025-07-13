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
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

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

                # Filter to only IFNUM == 0 and PLNUM == 0
                self.data = self.data[(self.data['IFNUM'] == 0) & (self.data['PLNUM'] == 0)]

                num_points = len(self.data[(self.data['CALSTATE'] == 0) & (self.data['SWPVALID'] == 0)])
                print(num_points)
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

            if pre_cal_complete and i["SWPVALID"] == 0 and self.data[ind - 1]["SWPVALID"] == 0:
                if post_cal_start_ind is None:
                    post_cal_start_ind = ind
            else:
                post_cal_start_ind = None

            if pre_cal_complete and i['CALSTATE'] == 0 and i['SWPVALID'] == 1:
                counter += 1

            if counter <= 3 * channel_count and i['SWPVALID'] == 1:
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


    def compute_gain_deltas(self, window):
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

        def extract(data_subset, swpvalid_val, calstate_val, i, j, h):
            rows = [
                ((row['time'] - t0).sec, row['DATA'][h])
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

        for i in ifnums:
            for j in plnums:
                pre_cal_deltas = []
                pre_cal_uncs = []
                post_cal_deltas = []
                post_cal_uncs = []

                # Step 1: Collect calibration deltas
                for h in range(first_len):
                    data = extract(data_array, 1, 0, i, j, h)

                    try:
                        pre_on_times, pre_on_vals = extract(pre_data, 0, 1, i, j, h)
                        pre_off_times, pre_off_vals = extract(pre_data, 0, 0, i, j, h)

                        if len(pre_on_times) == 0 or len(pre_off_times) == 0:
                            pre_cal_deltas.append(None)
                            pre_cal_uncs.append(None)
                            continue

                        pre_on_params, pre_on_unc = self.rcr((pre_on_times, pre_on_vals))
                        pre_off_params, pre_off_unc = self.rcr((pre_off_times, pre_off_vals))

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

                    except:
                        pre_cal_delta = None
                        pre_cal_unc = None

                    pre_cal_deltas.append(pre_cal_delta)
                    pre_cal_uncs.append(pre_cal_unc)

                    try:
                        post_on_times, post_on_vals = extract(post_data, 0, 1, i, j, h)
                        post_off_times, post_off_vals = extract(post_data, 0, 0, i, j, h)

                        if len(post_on_times) == 0 or len(post_off_times) == 0:
                            post_cal_deltas.append(None)
                            post_cal_uncs.append(None)
                            continue

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

                    except:
                        post_cal_delta = None
                        post_cal_unc = None

                    post_cal_deltas.append(post_cal_delta)
                    post_cal_uncs.append(post_cal_unc)

                t0 = Time(self.data[0]['DATE-OBS'], format='isot')
                time_array = np.array([(Time(row['DATE-OBS'], format='isot') - t0).sec for row in self.data])

                # Step 2: Apply averaged calibration deltas
                for h in range(first_len):
                    if pre_cal_deltas[h] is None or post_cal_deltas[h] is None:
                        continue

                    # Use 2 deltas on each side, ignoring None
                    def average_neighbors(values, index, window=window):
                        if window == 0:
                            return values[index] if values[index] is not None else None
                        neighbors = [
                            v for v in values[max(0, index - window): index + window + 1]
                            if v is not None
                        ]
                        return np.mean(neighbors) if neighbors else None

                    pre_avg = average_neighbors(pre_cal_deltas, h)
                    post_avg = average_neighbors(post_cal_deltas, h)
                    pre_unc_avg = average_neighbors(pre_cal_uncs, h)
                    post_unc_avg = average_neighbors(post_cal_uncs, h)

                    if pre_avg is None or post_avg is None or pre_unc_avg is None or post_unc_avg is None:
                        continue

                    z_value = abs(pre_avg - post_avg) / np.sqrt(pre_unc_avg**2 + post_unc_avg**2)

                    if z_value < 0.6745:
                        weights = np.array([1 / pre_unc_avg, 1 / post_unc_avg])
                        delta = np.average([pre_avg, post_avg], weights=weights)
                        self.data['DATA'][:, h] /= delta
                    else:
                        for ind, _ in enumerate(self.data):
                            # Linear interpolation of delta across time (fallback: midpoint time)
                            t = time_array[ind]
                            delta = pre_avg + (post_avg - pre_avg) * (t - pre_cal_time) / (post_cal_time - pre_cal_time)
                            self.data['DATA'][ind, h] /= delta

                    print(f"IFNUM={i}, PLNUM={j}, CH={h} → pre={pre_avg}±{pre_unc_avg}, post={post_avg}±{post_unc_avg}")

                # Convert to 2D array with shape (N, 2)
                deltas = np.column_stack((pre_cal_deltas, post_cal_deltas))

                # Use boolean mask to preserve structured array behavior
                mask = (self.data['IFNUM'] == 0) & (self.data['PLNUM'] == 0)
                filtered_data = self.data[mask]

                # Compute spectrum over time (average across rows)
                spectrum = np.mean(filtered_data['DATA'], axis=0)
                ch_indices = np.arange(len(spectrum))
                delta_diff = abs(np.array(pre_cal_deltas) - np.array(post_cal_deltas))
                self._pre_cal_deltas = pre_cal_deltas
                self._post_cal_deltas = post_cal_deltas





    def animate_gain_deltas(self, obj, window_range=(5, 105, 10), output_file="gain_deltas_window_scan.mp4"):
        """
        Animate the effect of different window sizes on gain delta smoothing.

        Parameters:
            obj: an instance of the class containing `compute_gain_deltas(window)`
            window_range: tuple (start, stop, step) for the smoothing window size
            output_file: filename for the output animation
        """

        # Store figures for reuse in animation
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        def update(frame_window):
            print(f"[INFO] Processing window = {frame_window}")
            ax1.clear()
            ax2.clear()

            # Call compute_gain_deltas with the current window
            try:
                obj.compute_gain_deltas(window=frame_window)
            except Exception as e:
                print(f"[ERROR] Failed at window={frame_window}: {type(e).__name__}: {e}")

            # Structured array filtering
            mask = (obj.data['IFNUM'] == 0) & (obj.data['PLNUM'] == 0)
            filtered_data = obj.data[mask]

            # Plot 1: Mean Spectrum
            spectrum = np.mean(filtered_data['DATA'], axis=0)
            ch_indices = np.arange(len(spectrum))
            ax1.plot(ch_indices, spectrum, color='black')
            ax1.set_ylim(0, 300)
            ax1.set_xlabel("Channel Index")
            ax1.set_ylabel("Mean Intensity")
            ax1.set_title(f"Mean Spectrum\nWindow={frame_window}")
            ax1.grid(True)

            # Plot 2: Delta Curve
            pre_deltas = obj._pre_cal_deltas  # Assume you store them as attributes
            post_deltas = obj._post_cal_deltas
            delta_diff = abs(np.array(pre_deltas) - np.array(post_deltas))

            ax2.plot(ch_indices, pre_deltas, label="Pre-Cal", color='blue', alpha=0.6, linewidth=0.4)
            ax2.plot(ch_indices, post_deltas, label="Post-Cal", color='orange', alpha=0.6, linewidth=0.4)
            ax2.plot(ch_indices, delta_diff, label="|Pre - Post|", color='red', alpha=0.8, linewidth=0.4)
            ax2.set_xlabel("Channel Index")
            ax2.set_ylabel("Delta")
            ax2.set_title("Calibration Deltas")
            ax2.legend()
            ax2.grid(True)

            ax1.relim()
            ax1.autoscale_view()

            ax2.relim()
            ax2.autoscale_view()
            

        # Create list of window sizes
        window_sizes = list(range(*window_range))

        ani = animation.FuncAnimation(fig, update, frames=window_sizes, repeat=False)
        ani.save(output_file, writer='ffmpeg', dpi=200)
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
        # self.compute_gain_deltas()
        self.animate_gain_deltas(self, window_range=(5, 105, 10), output_file="deltas_vs_window.mp4")



if __name__ == "__main__":
    filepath = "C:/Users/starb/Downloads/0136376.fits"
    file = Gain_Calibrate(filepath)

    file.gain_calibrate()