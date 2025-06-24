# Third-party libraries
import numpy as np
from astropy.time import Time


def average(data, axis):
    '''
    Average across time (axis = 0) or integrate across frequency (axis = 1).
    '''

    intensities = np.array([row[6] for row in data])
    count = intensities.shape[axis]
    channel_means = np.sum(intensities, axis=axis) / count

    return channel_means


def sdfits_to_array(data, t0):
    '''
    Convert sdfits data to more accessible, lighter arrays.

    Params:
    file: Radio_File class file
    data: some astropy FITS loaded data

    Return:
    2D array: times and frequencies
    '''

    freq = average(data, axis=1)

    times = Time(data["DATE-OBS"], format='isot')
    
    time_rel = (times - t0).sec

    return [time_rel, freq]