import numpy as np
from scipy.stats import linregress
from datetime import datetime
import rcr
from astropy.time import Time
from file_exception import MyException
from file_init import Mike
from val import Val
from sort import Sort
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


class Cal:
    def __init__(self, file):
        '''
        Initialize file.
        '''

        self.file = file
    

    def average(self, data, axis):
        '''
        Average across time (axis = 0) or integrate across frequency (axis = 1).
        '''

        intensities = np.array([row[6] for row in data])
        count = intensities.shape[axis]
        channel_means = np.sum(intensities, axis=axis) / count

        return channel_means
    

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
        y = array[1]

        if len(x) > 1 and len(y) > 1:
            result = linregress(x, y)
            m = result.slope
            b = result.intercept
        else:
            m = 0
            b = y[0]

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

        x = x[indices]
        y = y[indices]

        best_fit_parameters = model.result.parameters
        
        sigma = (1 / (len(x) - 2)) * np.sum((y - best_fit_parameters[1] * x - best_fit_parameters[0]) ** 2)
        m_sd = np.sqrt(sigma / np.sum((x - np.mean(x)) ** 2))
        b_sd = np.sqrt(sigma * ((1 / len(x)) + ((np.mean(x) ** 2) / np.sum((x - np.mean(x)) ** 2))))
        uncertainties = (b_sd, m_sd)

        return best_fit_parameters, uncertainties


    def sdfits_to_array(self, data):
        '''
        Convert sdfits data to more accessible, lighter arrays.

        Params:
        file: Mike class file
        data: some astropy FITS loaded data

        Return:
        2D array: times and frequencies
        '''

        freq = self.average(data, axis=1)

        times = Time(data["DATE-OBS"], format='isot')
        t0 = Time(self.file.header["DATE"], format='isot')

        # Time deltas in seconds
        time_rel = (times - t0).sec

        result = (np.array(time_rel), np.array(freq))

        return result
    

    def compute_gain_deltas(self, ind):
        '''
        Compute gain deltas with pre calibration and post calibration.

        Params:
        file: Mike class file
        ind: index of channel being processed
        '''
        
        for ind, i in enumerate(self.file.data):
            subset_data = self.file.data[ind]
            subset_indices = self.file.data_indicies[ind]

            try:
                pre_cal = subset_data[
                    (np.arange(len(subset_data)) < subset_indices[0]) &
                    (subset_data["SWPVALID"] == 0)
                ]
                try:
                    delta1, t1 = get_delta(pre_cal)
                    self.file.gain_start.append([delta1, t1])
                except:
                    self.file.gain_start.append(None)
            except:
                pre_cal = None
                self.file.gain_start.append(None)

            try:
                post_cal = subset_data[
                    (np.arange(len(subset_data)) >= subset_indices[1]) &
                    (subset_data["SWPVALID"] == 0)
                ]
                try:
                    delta2, t2 = get_delta(post_cal)
                    self.file.gain_end.append([delta2, t2])
                except:
                    self.file.gain_end.append(None)
            except:
                post_cal = None
                self.file.gain_end.append(None)
                return


            def get_delta(cal):
                cal_on_array = self.sdfits_to_array(cal[cal["CALSTATE"] == 1])
                cal_on_params, cal_on_uncertainties = self.rcr(cal_on_array)

                if len(cal[cal["CALSTATE"] == 0]):
                    cal_off_array = self.sdfits_to_array(cal[cal["CALSTATE"] == 0])
                else:
                    return None

                cal_off_params, cal_off_uncertainties = self.rcr(cal_off_array)

                time = (np.mean(cal_on_array[0]) + np.mean(cal_off_array[0])) / 2
                if time < (cal_on_array[0][0] + cal_off_array[0][-1]) / 2:
                    time = (cal_on_array[0][0] + cal_off_array[0][-1]) / 2
                elif time > (cal_on_array[0][-1] + cal_off_array[0][0]) / 2:
                    time = (cal_on_array[0][-1] + cal_off_array[0][0]) / 2

                delta = (cal_on_params[1] * time + cal_on_params[0]) - (cal_off_params[1] * time + cal_off_params[0])

                return delta, time
            
        return


if __name__ == "__main__":
    '''
    Test function to implement calibration.
    '''

    file = Mike("C:/Users/starb/Downloads/0136376.fits")
    v = Val(file)
    # v.validate_primary_header()
    # v.validate_data()
    
    s = Sort(file)
    s.sort()

    c = Cal(file)
    c.compute_gain_deltas()
    
    print(file.gain_start)
    print(file.gain_end)
