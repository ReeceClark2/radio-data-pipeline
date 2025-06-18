
from scipy.stats import linregress
import matplotlib as plt
from calendar import day_abbr
from astropy.io import fits
from astropy.table import Table
import numpy as np
from file_exception import MyException
import warnings
from file_init import Mike
from astropy.table import Column
import os
import time
from cal import Cal
from val import Val
from sort import Sort
from gain_calibration import Gain_Cal
from astropy.time import Time

class Flux_Cal:
    def __init__(self, file, calfile):
        """
        Initialize the flux_Calibration object with a FITS file.
            """
        
        self.file = file
        self.calfile = calfile

    
    def get_radio_calibrate_params(self, source):
        """
        Get parameters for calculating true flux in janskys of the calibration source

        param file: Mike class file of calibration source

        returns: parameters for equations 14 & 15 for particular calibration source
        """ 

        sources = {
            "CAS A": {
                "t_ref": 2006.9, "t_0": 2005.64,
                "logF_0": 3.2530, "a_1": -0.732, "nu_ref": 1477,
                "a_2": -0.0094, "a_3": 0.0053,
                "mnu_0": -0.00350, "mdeltlog": 0.00124, "nu_0": 1315,
                "var_logF_0": 0.0051**2, "var_a_1": 0.011**2,
                "var_a_2": 0.0058**2, "var_a_3": 0.0058**2,
                "var_mnu_0": 0.00022**2, "var_mdeltlog": 0.00018**2
            },
            "CYG A": {
                "t_ref": 0, "t_0": 0,
                "logF_0": 3.1861, "a_1": -1.038, "nu_ref": 1416,
                "a_2": -0.1457, "a_3": 0.0170,
                "mnu_0": 0, "mdeltlog": 0, "nu_0": 1,
                "var_logF_0": 0.0046**2, "var_a_1": 0.011**2,
                "var_a_2": 0.0075**2, "var_a_3": 0.0075**2,
                "var_mnu_0": 0, "var_mdeltlog": 0
            },
            "TAU A": {
                "t_ref": 2009.05, "t_0": 2009.05,
                "logF_0": 2.9083, "a_1": -0.226, "nu_ref": 1569,
                "a_2": -0.0113, "a_3": -0.0275,
                "mnu_0": -0.00044, "mdeltlog": 0, "nu_0": 1,
                "var_logF_0": 0.0044**2, "var_a_1": 0.014**2,
                "var_a_2": 0.0081**2, "var_a_3": 0.0077**2,
                "var_mnu_0": 0.00019**2, "var_mdeltlog": 0
            },
            "VIR A": {
                "t_ref": 0, "t_0": 0,
                "logF_0": 2.3070, "a_1": -0.876, "nu_ref": 1482,
                "a_2": -0.047, "a_3": -0.073,
                "mnu_0": 0, "mdeltlog": 0, "nu_0": 1,
                "var_logF_0": 0.0045**2, "var_a_1": 0.017**2,
                "var_a_2": 0.0031**2, "var_a_3": 0.0030**2,
                "var_mnu_0": 0, "var_mdeltlog": 0
            }
        }

        return sources.get(source.upper(), None)
    

    def average(self, data, axis):
        '''
        Average across time (axis = 0) or integrate across frequency (axis = 1).
        '''

        intensities = np.array([row[6] for row in data])
        count = intensities.shape[axis]
        channel_means = np.sum(intensities, axis=axis) / count

        return channel_means


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
        
        time_rel = (times - t0).sec

        return [time_rel, freq]
    
    

    def on_off_background_subtract(self, onoff_file): 
        """
        Get flux value of the calibration source in noise calibration units

        onoff_file: Mike class on off file of calibration source

        returns: list of fluxes in noise calibration units, each list value corresponds to the value in a particular channel
        """
     
        index = onoff_file.data_indicies
        calint= []
        for channel, i in enumerate(self.file.data):
            onstart = index[channel][0]
            onend = index[channel][1]
            offstart = index[channel][2]
            offend = index[channel][3]

            ondata = self.sdfits_to_array(onoff_file.data[channel][onstart:onend])
            offdata = self.sdfits_to_array(onoff_file.data[channel][offstart:offend])

            onmedian = np.median(ondata)
            offmedian = np.median(offdata)

            calint.append(onmedian - offmedian) # Calibration source flux in noise calibration units

        return calint
    

    def radio_calibrate(self, start_freq, end_freq, date, source):
        """
        Get flux value of the observed source in Janskys

        start_freq: start frequency of observation
        end_freq: end frequency of observation
        date: date of observation, in datetime format
        source: calibration source, all caps, eg. CYG A

        returns: 
        flux: average flux in janskyes of the calibration source
        fluxError: uncertainty in flux value
        effectiveFrequency: effective frequency of source        
        """
        
        if start_freq >= end_freq or start_freq < 0 or end_freq < 0:
            return {"flux": None, "fluxError": None, "effectiveFrequency": None}

        params = self.get_radio_calibrate_params('CYG A') # TEMPORARY
        if params is None:
            return {"flux": None, "fluxError": None, "effectiveFrequency": None}

        t = date.year + date.month / 12 + date.day / 30
        nu = np.linspace(start_freq, end_freq, 100001)
        log_nu_nuref = np.log(nu / params["nu_ref"])
        log_nu_nu0 = np.log(nu / params["nu_0"])

        eq14 = (
            params["logF_0"] + params["a_1"] * log_nu_nuref +
            params["a_2"] * log_nu_nuref**2 + params["a_3"] * log_nu_nuref**3 +
            (params["mnu_0"] * (t - params["t_ref"]) + params["mdeltlog"] * (t - params["t_0"]) * log_nu_nu0)
        )

        eq15 = np.sqrt(
            params["var_logF_0"] + params["var_a_1"] * log_nu_nuref**2 +
            params["var_a_2"] * log_nu_nuref**4 + params["var_a_3"] * log_nu_nuref**6 +
            params["var_mnu_0"] * (t - params["t_ref"])**2 +
            params["var_mdeltlog"] * log_nu_nu0**2 * (t - params["t_0"])**2
        )

        flux = 10 ** eq14
        sigma_flux = 10 ** (eq14 + eq15)
        effective_freq = nu * flux

        # Trapezoidal integration
        final_avg_flux = np.trapezoid(flux, nu) / (end_freq - start_freq)
        final_sigma_flux = np.trapezoid(sigma_flux, nu) / (end_freq - start_freq)
        uncertainty = final_sigma_flux - final_avg_flux
        final_effective_freq = np.trapezoid(effective_freq, nu) / (final_avg_flux * (end_freq - start_freq))
        
        return {
            "flux": final_avg_flux,
            "fluxError": uncertainty,
            "effectiveFrequency": final_effective_freq
        }
              
            
    def flux_cal(self):
        """
        Overwrites the continuum field with flux calibrated intensities (in Janskys)
        
        file: file path of the observation source
        calfile: file path of calibration source

        sets Flux_Calibrated flag to True
        """
        # List of median intensities of the calibration source
        calints = self.on_off_background_subtract(self.calfile)

        # Go through each data channel and calibrate the flux
        calib_flux_data = []

        for i, feed in enumerate(self.file.continuum):
            
            # Get an array of the continuum data for source
            obsint = feed[1]
            
            # Get flux of calibration source in noise source units
            calint = calints[i]
            
            # Get flux of calibration source in janskys
            date = Time(self.calfile.header["DATE-OBS"], format='isot')
            calintj = self.radio_calibrate(1355, 1435, date, self.calfile.header["object"].upper())["flux"] 

            # Conversion factor
            j_conversion = calintj / calint  

            # For the time array in the data find the calibrated height for each intensity
            new_intensities = [obsint * j_conversion for obsint in obsint]
            calib_flux_data.append(new_intensities)
                
        for ind, i in enumerate(self.file.continuum):
            old_continuum = self.file.continuum[ind]
            new_continuum = (old_continuum[0], calib_flux_data[ind])
            self.file.continuum[ind] = new_continuum

        self.file.flux_calibrated = True

        return


if __name__ == "__main__":
    '''
    Test function to implement calibration.
    '''    

    file = Mike("C://Users//leesnow//Downloads//0136376.fits")
    calfile = Mike("C://Users//leesnow//Downloads//0117613.fits")

    v = Val(file)
    vc = Val(calfile)
    
    s = Sort(file)
    s.sort()
    
    sc = Sort(calfile)
    sc.sort
    
    c = Cal(file)
    c.compute_gain_deltas()

    cc = Cal(calfile)
    cc.compute_gain_deltas()
    
    Data = Gain_Cal(file)
    Data.gain_cal(file)

    Datac = Gain_Cal(calfile)
    Datac.gain_cal(calfile)

    FData = Flux_Cal(file, calfile)
    FData.flux_cal(file, calfile)

