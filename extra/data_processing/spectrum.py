# Third-party libraries
import numpy as np

# Local application imports
from .child_init import Radio_Child_File
from .file_init import Radio_File
from .gain_calibration import Gain_Cal
from .sort import Sort
from .utils import average
from .val import Val


class Spectrum:
    def __init__(self, file):
        '''
        Initialize file.
        '''

        self.file = file


    def make_spec(self):
        '''
        Generate a spectrum plot for each channel in the data.

        params: file: Radio_File class file

        returns: populates the file's spectrum field with frequency and intensity data
        '''
        # Go through each channel in the data
        for ind, i in enumerate(self.file.data):
            # Only get actual data, not the calibration data
            data = self.file.data[ind][self.file.data_indices[ind][0]:self.file.data_indices[ind][-1]]

            # Get the summed intensities for each frequency
            result = average(data, 0)

            # Reverse the list becaise the data is in reverse order
            result = result[::-1]

            # Get the frequency range for the channel based on the feed
            # if the spectrum field has not been fully populated yet create the frequency range from the original start and stop frequencies
            frequencies = np.linspace(self.file.freqs[ind][0], self.file.freqs[ind][1], len(result))

            # Add the frequency and result to the file's spectrum field
            self.file.spectrum[ind] = ([np.array(frequencies), np.array(result)])

        return


if __name__ == "__main__":
    '''
    Test function to implement calibration.
    '''
    keep_indices = [[1300, 1400], [1406, 1410], [1412, 1420]]  # Specify the indices you want to keep
    feed= [1]  # Specify the feeds you want to keep
    fits_path = "TrackingHighRes/0136484.fits"
    file = Radio_File(fits_path)
    v = Val(file)
    # v.validate_primary_header()
    # v.validate_data()
    
    s = Sort(file)
    s.split_slp_feed()
    s.sort_data()
    s.get_startend_freqs()
    s.get_startstop_channels()

    spec = Spectrum(file)
    spec.make_spec()

    if keep_indices != []:
        child = Radio_Child_File(file)
        sortchild = Sort(child)
        specchild = Spectrum(child)

        sortchild.user_cuts(keep_indices, "spectrum", "cut", feed)
        specchild.make_spec(keep_indices, feed)
        print(file.data[0]["DATA"].shape)
        print(child.data[0]["DATA"].shape)

    c = Gain_Cal(file)
    c.compute_gain_deltas()
