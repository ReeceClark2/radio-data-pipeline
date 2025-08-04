# Third-party libraries
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from scipy.stats import linregress
import rcr
import numpy as np

import time

class Gain_Calibrate:
    def __init__(self, filepath):
        # filepath must route to a validate file so the required index fields exist!
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
        except Exception as e:
            raise

    def linear(self, x, params): # model function
        return params[0] + x * params[1]

    def d_linear_1(self, x, params): # first model parameter derivative
        return 1

    def d_linear_2(self, x, params): # second model parameter derivative
        return x

    def perform_rcr(self, array):
        # Perform functional Robust Chauvenet Rejection (RCR) on a given dataset.
        
        # Parse 2D array into x and y arrays
        x = array[0]
        x -= np.average(x)

        y = array[1]

        # Using linear regression, create a guess function for RCR
        result = linregress(x, y)
        m = result.slope
        b = result.intercept
        guess = [m, b]

        # Perform RCR
        model = rcr.FunctionalForm(self.linear,
            x,
            y,
            [self.d_linear_1, self.d_linear_2],
            guess
        )
        
        r = rcr.RCR(rcr.SS_MEDIAN_DL) 
        r.setParametricModel(model)
        r.performBulkRejection(y)

        # Fetch indices
        indices = r.result.indices

        # Keep on valid indices
        x = np.array([x[i] for i in indices])
        y = np.array([y[i] for i in indices])

        best_fit_parameters = model.result.parameters
        
        # Using the best fit parameters, calculate a slope and y-intercept uncertainty
        sigma = (1 / (len(x) - 2)) * np.sum((y - best_fit_parameters[1] * x - best_fit_parameters[0]) ** 2)
        m_sd = np.sqrt(sigma / np.sum((x - np.mean(x)) ** 2))
        b_sd = np.sqrt(sigma * ((1 / len(x)) + ((np.mean(x) ** 2) / np.sum((x - np.mean(x)) ** 2))))
        uncertainties = (b_sd, m_sd)

        return best_fit_parameters, uncertainties

    def get_continuum(self, data):
        # Get intensities and shape of provided Astropy data table snippet
        intensities = np.array(data['DATA']) 
        count = intensities.shape[1]
        channel_means = np.sum(intensities, axis=1) / count

        # Use the headers start time to 0 the time array of the observation
        times = Time(data["DATE-OBS"], format='isot')
        t0 = Time(self.header["DATE"], format="isot")
        time_rel = (times - t0)

        return (time_rel.sec, channel_means)

    def get_cal_arrays(self, ifnum, plnum):
        # Filter out the pre cal on and off arrays
        try:
            pre_cal_on_mask = self.data[:self.data_start_index][
                (self.data['CALSTATE'][:self.data_start_index] == 1) &
                (self.data['SWPVALID'][:self.data_start_index] == 0) & 
                (self.data['IFNUM'][:self.data_start_index] == ifnum) &
                (self.data['PLNUM'][:self.data_start_index] == plnum)
            ]

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

            pre_cal_off_array = self.get_continuum(pre_cal_off_mask)
        except Exception:
            pre_cal_off_array = None

        # Filter out the post cal on and off arrays
        try:
            post_cal_on_mask = self.data[self.post_calibration_start_index:][
                (self.data['CALSTATE'][self.post_calibration_start_index:] == 1) &
                (self.data['SWPVALID'][self.post_calibration_start_index:] == 0) & 
                (self.data['IFNUM'][self.post_calibration_start_index:] == ifnum) &
                (self.data['PLNUM'][self.post_calibration_start_index:] == plnum)
            ]

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

            post_cal_off_array = self.get_continuum(post_cal_off_mask)
        except Exception:
            post_cal_off_array = None

        # Return all arrays to be iterated over
        return [[pre_cal_on_array, pre_cal_off_array], [post_cal_on_array, post_cal_off_array]]

    def gain_calibrate_channel(self, channel, ifnum, plnum):
        # Get the start time and array of times through the channel
        t0 = Time(self.header["DATE"], format="isot")
        time_array = (Time(channel['DATE-OBS'], format='isot') - t0).sec

        cal_arrays = self.get_cal_arrays(ifnum, plnum)
        for ind, j in enumerate(cal_arrays):
            # For pre and post cal in cal_arrays try to compute the delta
            try:
                # j[0] is the on array for a calibration spike while j[1] is the off array. j[0][0] are the times of the on array and j[0][1] are the intensities
                if j[0][0].size != 0 and j[0][1].size != 0 and j[1][0].size != 0 and j[1][1].size != 0:
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
                else:
                    cal_delta = None
                    cal_uncertainty = None
            except:
                cal_delta = None
                cal_uncertainty = None

            if ind == 0:
                pre_cal_delta = cal_delta
                pre_cal_uncertainty = cal_uncertainty
                
            elif ind == 1:
                post_cal_delta = cal_delta
                post_cal_uncertainty = cal_uncertainty
    

        if pre_cal_delta and post_cal_delta:
            z_value = abs(pre_cal_delta - post_cal_delta) / np.sqrt(
                pre_cal_uncertainty ** 2 + post_cal_uncertainty ** 2
            )

            if z_value < 0.6745:
                weights = np.array([1 / pre_cal_uncertainty, 1 / post_cal_uncertainty])
                delta = np.average([pre_cal_delta, post_cal_delta], weights=weights)
                channel['DATA'] /= delta
            else:
                for ind, _ in enumerate(channel['DATA']):
                    t = time_array[ind]
                    delta = pre_cal_delta + (post_cal_delta - pre_cal_delta) * (t - t0)
                    channel['DATA'][ind] /= delta

        elif pre_cal_delta and not post_cal_delta:
            channel['DATA'] /= pre_cal_delta

        elif post_cal_delta and not pre_cal_delta:
            channel['DATA'] /= post_cal_delta

        mask = (self.data['IFNUM'] == ifnum) & (self.data['PLNUM'] == plnum)
        self.data['DATA'][mask] = channel['DATA']


    def gain_calibrate(self):
        ifnums = np.unique(self.data['IFNUM'])
        plnums = np.unique(self.data['PLNUM'])
    
        for i in ifnums:
            for j in plnums:
                mask = (self.data['IFNUM'] == i) & (self.data['PLNUM'] == j)
                channel = self.data[mask]
                self.gain_calibrate_channel(channel, i, j)
                
            
if __name__ == "__main__":
    filepath = "C:/Users/starb/Downloads/0136873_validate.fits"
    file = Gain_Calibrate(filepath)
    
    start = time.time() 
    file.gain_calibrate()
    end = time.time()

    print(f"Elapsed time: {end - start:.6f} seconds")