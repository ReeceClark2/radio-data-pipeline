# Third-party libraries
import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.table import Column
from astropy.time import Time
from scipy.stats import linregress
import rcr
import time
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import copy
import math

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

                self.data_start_index = self.header["DATAIND"]
                self.post_calibration_start_index = self.header["POSTIND"]
                self.off_start_index = self.header["ONOFFIND"]

                # Filter to only IFNUM == 0 and PLNUM == 0
                # self.data = self.data[(self.data['IFNUM'] == 0) & (self.data['PLNUM'] == 0)]

            # TODO add logging here
        except Exception as e:
            # TODO add logging here
            pass


    def linear(self, x, params): # model function
        return params[0] + x * params[1]


    def d_linear_1(self, x, params): # first model parameter derivative
        return 1


    def d_linear_2(self, x, params): # second model parameter derivative
        return x


    def perform_rcr(self, array):
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


    def get_continuum(self, data):
        intensities = np.array(data['DATA']) 
        count = intensities.shape[1]
        channel_means = np.sum(intensities, axis=1) / count

        times = Time(data["DATE-OBS"], format='isot')
        t0 = Time(self.header["DATE"], format="isot")
        time_rel = (times - t0)

        return (time_rel.sec, channel_means)
    

    def get_spectrum(self, data):
        intensities = np.array(data['DATA']) 
        count = intensities.shape[0]
        channel_means = np.sum(intensities, axis=0) / count

        return channel_means
    

    def get_cal_arrays(self, start, end, ifnum, plnum):
        try:
            pre_cal_on_mask = self.data[:self.data_start_index][
                (self.data['CALSTATE'][:self.data_start_index] == 1) &
                (self.data['SWPVALID'][:self.data_start_index] == 0) & 
                (self.data['IFNUM'][:self.data_start_index] == ifnum) &
                (self.data['PLNUM'][:self.data_start_index] == plnum)
            ]
            pre_cal_on_mask['DATA']  = [row[start:end] for row in pre_cal_on_mask['DATA']]

            pre_cal_on_array = self.get_continuum(pre_cal_on_mask)
        except Exception:
            pre_cal_on_array = None

        try:
            pre_cal_off_mask = self.data[:self.data_start_index][
                (self.data['CALSTATE'][:self.data_start_index] == 0) &
                (self.data['SWPVALID'][:self.data_start_index] == 0) & 
                (self.data['IFNUM'][:self.data_start_index] == ifnum) &
                (self.data['PLNUM'][:self.data_start_index] == plnum)
            ]
            pre_cal_off_mask['DATA'] = [row[start:end] for row in pre_cal_off_mask['DATA']]

            pre_cal_off_array = self.get_continuum(pre_cal_off_mask)
        except Exception:
            pre_cal_off_array = None

        try:
            post_cal_on_mask = self.data[self.post_calibration_start_index:][
                (self.data['CALSTATE'][self.post_calibration_start_index:] == 1) &
                (self.data['SWPVALID'][self.post_calibration_start_index:] == 0) & 
                (self.data['IFNUM'][self.post_calibration_start_index:] == ifnum) &
                (self.data['PLNUM'][self.post_calibration_start_index:] == plnum)
            ]
            post_cal_on_mask['DATA'] = [row[start:end] for row in post_cal_on_mask['DATA']]

            post_cal_on_array = self.get_continuum(post_cal_on_mask)
        except Exception:
            post_cal_on_array = None

        try:
            post_cal_off_mask = self.data[self.post_calibration_start_index:][
                (self.data['CALSTATE'][self.post_calibration_start_index:] == 0) &
                (self.data['SWPVALID'][self.post_calibration_start_index:] == 0) &
                (self.data['IFNUM'][self.post_calibration_start_index:] == ifnum) &
                (self.data['PLNUM'][self.post_calibration_start_index:] == plnum)
            ]
            post_cal_off_mask['DATA'] = [row[start:end] for row in post_cal_off_mask['DATA']]

            post_cal_off_array = self.get_continuum(post_cal_off_mask)
        except Exception:
            post_cal_off_array = None


        return [[pre_cal_on_array, pre_cal_off_array], [post_cal_on_array, post_cal_off_array]]
    

    def gain_calibrate_channel(self, channel):
        ifnum = channel['IFNUM'][0]
        plnum = channel['PLNUM'][0]

        channel_count = len(channel['DATA'][0])
        t0 = (Time(self.header["DATE"], format="isot")).sec
        time_array = (Time(channel['DATE-OBS'], format='isot')).sec

        bin_count = 1
        bin_size = channel_count / bin_count  # use integer division

        pre_cal_deltas = np.zeros(bin_count)
        pre_cal_uncertainties = np.zeros(bin_count)

        post_cal_deltas = np.zeros(bin_count)
        post_cal_uncertainties = np.zeros(bin_count)

        for i in range(bin_count):
            start = round(i * bin_size)
            end = round((i + 1) * bin_size)
            
            bin_data = channel.copy()
            bin_data['DATA'] = [row[start:end] for row in channel['DATA']]
            bin_cal_arrays = self.get_cal_arrays(start, end, ifnum, plnum)

            for ind, j in enumerate(bin_cal_arrays):
                try:
                    cal_on_params, on_uncertainty = self.perform_rcr(j[0])
                    cal_off_params, off_uncertainty = self.perform_rcr(j[1])

                    cal_time = np.mean([np.mean(j[0][0]), np.mean(j[1][0])])
                    cal_time = np.clip(
                        cal_time,
                        (j[0][0][0] + j[1][0][-1]) / 2,
                        (j[0][0][-1] + j[1][0][0]) / 2
                    )

                    cal_delta = np.abs(
                        (cal_on_params[1] * (cal_time - np.mean(j[0][0])) + cal_on_params[0]) -
                        (cal_off_params[1] * (cal_time - np.mean(j[1][0])) + cal_off_params[0])
                    )
                    cal_uncertainty = np.sqrt(
                        on_uncertainty[0]**2 + off_uncertainty[0]**2 +
                        (on_uncertainty[1] * cal_time)**2 + (off_uncertainty[1] * cal_time)**2
                    )
                except:
                    cal_delta = None
                    cal_uncertainty = None
                
                if ind == 0:
                    pre_cal_deltas[i] = cal_delta
                    pre_cal_uncertainties[i] = cal_uncertainty
                elif ind == 1:
                    post_cal_deltas[i] = cal_delta
                    post_cal_uncertainties[i] = cal_uncertainty

        def is_valid(x):
            return x is not None and not np.isnan(x)
        
        plt.plot(pre_cal_deltas)
        plt.plot(post_cal_deltas)
        # plt.ylim(0, 50)
        plt.savefig('delta valuess')
        print('gain calibrated successfully')
        plt.close()


        num_points = len(channel['DATA'][0])  # number of data columns
        for i in range(num_points):
            bin_index = int(i // bin_size)
            frac = (i % bin_size) / bin_size  # interpolation fraction within bin

            pre_cal_delta = None
            post_cal_delta = None

            # interpolate pre_cal_delta if possible
            if bin_index < len(pre_cal_deltas) - 1 and is_valid(pre_cal_deltas[bin_index]) and is_valid(pre_cal_deltas[bin_index + 1]):
                # linear interpolation between two bins
                pre_cal_delta = pre_cal_deltas[bin_index] * (1 - frac) + pre_cal_deltas[bin_index + 1] * frac
            elif bin_index < len(pre_cal_deltas):
                pre_cal_delta = pre_cal_deltas[bin_index]

            # interpolate post_cal_delta if possible
            if bin_index < len(post_cal_deltas) - 1 and is_valid(post_cal_deltas[bin_index]) and is_valid(post_cal_deltas[bin_index + 1]):
                post_cal_delta = post_cal_deltas[bin_index] * (1 - frac) + post_cal_deltas[bin_index + 1] * frac
            elif bin_index < len(post_cal_deltas):
                post_cal_delta = post_cal_deltas[bin_index]

            # now apply the logic
            if is_valid(pre_cal_delta) and is_valid(post_cal_delta):
                z_value = abs(pre_cal_delta - post_cal_delta) / np.sqrt(
                    pre_cal_uncertainties[i] ** 2 + post_cal_uncertainties[i] ** 2
                )

                if z_value < 0.6745:
                    weights = np.array([1 / pre_cal_uncertainties[i], 1 / post_cal_uncertainties[i]])
                    delta = np.average([pre_cal_delta, post_cal_delta], weights=weights)
                    channel['DATA'][:, i] /= delta
                else:
                    for ind, _ in enumerate(channel['DATA']):
                        t = time_array[ind]
                        # linear interpolation over time
                        delta = pre_cal_delta + (post_cal_delta - pre_cal_delta) * (t - t0)
                        channel['DATA'][ind, i] /= delta

            elif is_valid(pre_cal_delta):
                channel['DATA'][:, i] /= pre_cal_delta

            elif is_valid(post_cal_delta):
                channel['DATA'][:, i] /= post_cal_delta

            else:
                continue

        
        spectrum = self.get_spectrum(channel)

        plt.plot(spectrum)
        # plt.ylim(0, 50)
        plt.savefig('delta valuess applied')
        print('gain calibrated successfully')
        plt.close()
        exit()



    def gain_calibrate(self):
        ifnums = np.unique(self.data['IFNUM'])
        plnums = np.unique(self.data['PLNUM'])
    
        for i in ifnums:
            for j in plnums:
                mask = (self.data['IFNUM'] == i) & (self.data['PLNUM'] == j)
                channel = self.data[mask]
                self.gain_calibrate_channel(channel)

            

if __name__ == "__main__":
    filepath = "C:/Users/starb/Downloads/0136873_validate.fits"
    file = Gain_Calibrate(filepath)

    
    start = time.time()
    file.gain_calibrate()
    end = time.time()

    print(f"Elapsed time: {end - start:.6f} seconds")
