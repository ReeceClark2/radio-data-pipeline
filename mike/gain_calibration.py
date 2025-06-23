# Standard library
import os

# Third-party libraries
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from scipy.stats import linregress

# Local application imports
import rcr
from file_init import Radio_File
from sort import Sort
from utils import sdfits_to_array
from val import Val
from file_exception import MyException



class Gain_Cal:
    """
    Class to handle continuum data from FITS files and perform gain calibration.
    Attributes:
        file: An object representing the FITS file, expected to have a 'file_path' attribute and other calibration-related attributes.
        dataH: The header of the FITS file's primary data extension.
        data: The data table extracted from the FITS file.

    Methods:
        __init__(file):
            Initializes the gain_calibration object by loading the FITS file header and data table.
        calib_Heights(file):
            Calibrates the continuum data using gain start and end values.
            For each data entry, applies gain calibration based on the presence of gain start and/or gain end values.
            Updates the data in-place and sets the 'Gain_Calibrated' flag accordingly.
    """
    """Class to handle continuum data from FITS files."""
    def __init__(self, file):
        """
        Initialize the Gain_Calibration object with a FITS file.
        """

        self.file = file
    

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
            guess = [m, b]
            error_y = np.full(len(x), 10**-6)
            model = rcr.FunctionalForm(self.linear,
                x,
                y,
                [self.d_linear_1, self.d_linear_2],
                guess, error_y = error_y
            )
            
            r = rcr.RCR(rcr.SS_MEDIAN_DL) 
            r.setParametricModel(model)
            r.performBulkRejection(y)

            indicies = r.result.indicies

            x = np.array([x[i] for i in indicies])
            y = np.array([y[i] for i in indicies])

            best_fit_parameters = model.result.parameters
            best_fit_uncertainties = model.result.parameter_uncertainties    
            rcr_uncertainty = best_fit_uncertainties
            
            sigma = (1 / (len(x) - 2)) * np.sum((y - best_fit_parameters[1] * x - best_fit_parameters[0]) ** 2)
            m_sd = np.sqrt(sigma / np.sum((x - np.mean(x)) ** 2))
            b_sd = np.sqrt(sigma * ((1 / len(x)) + ((np.mean(x) ** 2) / np.sum((x - np.mean(x)) ** 2))))
            uncertainties = (b_sd, m_sd)

            return best_fit_parameters, uncertainties

        else:
            m = 0
            b = y[0]

        guess = [m, b]
        error_y = np.full(len(x), 10**-6)
        model = rcr.FunctionalForm(self.linear,
            x,
            y,
            [self.d_linear_1, self.d_linear_2],
            guess , error_y = error_y
        )
        
        r = rcr.RCR(rcr.SS_MEDIAN_DL) 
        r.setParametricModel(model)
        r.performBulkRejection(y)

        indices = r.result.indices

        x = np.array([x[i] for i in indicies])
        y = np.array([y[i] for i in indicies])

        best_fit_parameters = model.result.parameters
        sigma = (1 / (len(x) - 2)) * np.sum((y - best_fit_parameters[1][0] * x - best_fit_parameters[0][0]) ** 2)
        m_sd = np.sqrt(sigma / np.sum((x - np.mean(x)) ** 2))
        b_sd = np.sqrt(sigma * ((1 / len(x)) + ((np.mean(x) ** 2) / np.sum((x - np.mean(x)) ** 2))))
        uncertainties = (b_sd, m_sd)
        
        return best_fit_parameters, uncertainties
    

    def compute_gain_deltas(self):
        '''
        Compute gain deltas with pre calibration and post calibration.

        Params:
        file: Radio_File class file
        ind: index of channel being processed
        '''
        def get_delta(cal):
            t0 = Time(self.header["DATE"], format="isot")
            cal_on_array = sdfits_to_array(cal[cal["CALSTATE"] == 1], t0)
            cal_on_params, cal_on_uncertainties = self.rcr(cal_on_array)
            
            if len(cal[cal["CALSTATE"] == 0]):
                cal_off_array = sdfits_to_array(cal[cal["CALSTATE"] == 0], t0)
            else:
                return None
            
            cal_off_params, cal_off_uncertainties = self.rcr(cal_off_array)
            
            time = (np.mean(cal_on_array[0]) + np.mean(cal_off_array[0])) / 2
            if time < (cal_on_array[0][0] + cal_off_array[0][-1]) / 2:
                time = (cal_on_array[0][0] + cal_off_array[0][-1]) / 2
            elif time > (cal_on_array[0][-1] + cal_off_array[0][0]) / 2:
                time = (cal_on_array[0][-1] + cal_off_array[0][0]) / 2
            
            delta = np.abs((cal_on_params[1] * (time - np.mean(cal_on_array[0])) + cal_on_params[0]) - (cal_off_params[1] * (time - np.mean(cal_off_array[0])) + cal_off_params[0]))
            
            #Sigma accounting for uncertainty in slope and y intercept from least squares
            delta_uncertainty = np.sqrt(np.square(cal_on_uncertainties[0]) + np.square(cal_off_uncertainties[0]) + np.square(cal_on_uncertainties[1]) * np.square(time) + np.square(cal_off_uncertainties[1]) * np.square(time))
            
            return delta, time, delta_uncertainty
        
                        
        for ind, i in enumerate(self.file.data):
            subset_data = self.file.data[ind]
            subset_indicies = self.file.data_indicies[ind]
            try:
                pre_cal = subset_data[
                    (np.arange(len(subset_data)) < subset_indicies[0]) &
                    (subset_data["SWPVALID"] == 0)
                ]
                
                #Check if there are enough data points for calibration
                if len(pre_cal) < 3:
                    self.file.gain_start.append(None)
                    raise MyException(f"Pre calibration data for channel {ind} does not have enough points for calibration.")
                    
                else:
                    delta1, t1, sigma1 = get_delta(pre_cal)
                    self.file.gain_start.append([delta1, t1, sigma1])
            except:
                self.file.gain_start.append(None)

            try:
                post_cal = subset_data[
                    (np.arange(len(subset_data)) >= subset_indicies[-1]) &
                    (subset_data["SWPVALID"] == 0)
                ]
                #Check if there are enough data points for calibration
                if len(post_cal) < 3:
                    self.file.gain_end.append(None)
                    raise MyException(f"Post calibration data for channel {ind} does not have enough points for calibration.")
                       
                else:
                    delta2, t2, sigma2 = get_delta(post_cal)
                    self.file.gain_end.append([delta2, t2, sigma2])

            except:
                self.file.gain_end.append(None)


        return


    def cal_heights(self, feeds=None):
        """
        calibrate the continuum data.

        param file: Radio_File class file

        returns: adds calibrated height data to the file's continuum
        """

        t0 = Time(self.file.header["DATE"], format='isot')
        # Go through each data channel and calibrate the heights
        
        for ind, i in enumerate(self.file.data):
            feednum = np.unique(self.file.data[ind]["IFNUM"])[0]
            pol = np.unique(self.file.data[ind]["PLNUM"])[0]
            if feeds is not None and feednum not in feeds:
                continue  

            calib_height_data = []
            # First check if the gain start and end values are present
            
            if self.file.gain_start[ind] is not None and self.file.gain_end[ind] is not None:
                # Get all the gain values
                delta1 = self.file.gain_start[ind][0]
                delta2 = self.file.gain_end[ind][0]
                time1 = self.file.gain_start[ind][1]
                time2 = self.file.gain_end[ind][1]
                sigma1 = self.file.gain_start[ind][2]
                sigma2 = self.file.gain_end[ind][2]
                
                #Get z value
                z_value = abs(delta1 - delta2) / np.sqrt(sigma1**2 + sigma2**2)

                # Get an array of the continuum data
                data = self.file.data[ind]
                data = sdfits_to_array(data, t0)

                # For the time array in the data find the calibrated height for each intensity
                #If z_value is below 0.6745 average the cal heights, otherwise interpolate between them
                if z_value < 0.6745:

                    for ind, i in enumerate(data[0]):
                        delta = (delta1 + delta2) / 2
                        # Add each calibration height to the calib_height list
                        data[1][ind] = (data[1][ind] / delta)

                else:
                    for ind, i in enumerate(data[0]):
                        delta = delta1 + (delta2 - delta1) * (data[0][ind] - time1) / (time2 - time1)
                        
                        # Add each calibration height to the calib_height list
                        data[1][ind] = (data[1][ind] / delta)
                
                calib_height_data = data
                self.file.gain_calibrated.append(ind)
                
            # If gain_start is None and gain_end is not None, use gain_end for calibration
            elif self.file.gain_start[ind] is None and self.file.gain_end[ind] is not None:
                delta = self.file.gain_end[ind][0]

                data = self.file.data[ind]
                data = sdfits_to_array(data, t0)

                temp = data[1] / delta
                calib_height_data = [data[0], temp]
                self.file.gain_calibrated.append(ind)

            # If gain_start is present but gain_end is not, use gain_start for calibration
            elif self.file.gain_start[ind] is not None and self.file.gain_end[ind] is None:
                delta = self.file.gain_start[ind][0]

                data = self.file.data[ind]
                data = sdfits_to_array(data, t0)

                data[1] = data[1] / delta
                calib_height_data = data
                self.file.gain_calibrated.append(ind)

            # If gain_start is present but gain_end is not, use gain_start for calibration
            elif self.file.gain_start[ind] is None and self.file.gain_end[ind] is None:
                data = self.file.data[ind]
                data = sdfits_to_array(data, t0)

                calib_height_data = data

                # If both gain_start and gain_end are None, just copy the original data
                # But gotta turn it into a list of time, intensity lists

                self.file.continuum[ind] = (calib_height_data)
                self.file.gain_calibrated.append(ind)

                continue

            # Add the calibrated height data to continuum
        
            self.file.continuum.append(calib_height_data)

        self.file.logger.info(f"Gain calibrated {len(self.file.gain_calibrated)} of {len(self.file.data)} channels.")

        return
    

    def gain_cal(self):
        self.compute_gain_deltas()
        self.cal_heights()
        return


if __name__ == "__main__":
    """Test function to implement continuum data handling."""

    keep_times = [[8,12], [1404, 1410], [1412, 1420]]  # Specify the indicies you want to keep
    feed= [0]  # Specify the feeds you want to keep

    file = Radio_File("C://Users//leesnow//Downloads//0136376.fits")
    data = Gain_Cal(file)

    v = Val(file)
    # v.validate_primary_header()
    # v.validate_data()

    s = Sort(file)
    s.sort()

    #c.compute_gain_deltas()

    g = Gain_Cal(file)
    g.gain_cal()
    

    plt.savefig("testplot")
"""     if keep_times != []:
        child = Radio_Child_File(file)
        sortchild = Sort(child)
        contChild = Gain_Cal(child)

        sortchild.user_cuts(keep_times, "continuum", "cut", feed)
        g.gain_cal(file, feed)

        print(file.data[0]["DATA"].shape)
        print(child.data[0]["DATA"].shape) """
