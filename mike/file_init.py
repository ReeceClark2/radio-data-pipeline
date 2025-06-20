# Third-party libraries
import numpy as np
from astropy.io import fits
from astropy.table import Table

# Local application imports
from file_exception import MyException
from logger import Log_Collector


class Radio_File:
    def __init__(self, file_path):
        '''
        Initialize binary fits cube by reading its header and data from a given file.
        
        File fields:
        missing values: empty primary header values
        file_path: path to file in directory

        validated_header: check whether header has been validated
        validated_data: check whether data has been validated

        data_indicies: list of arrays of shape (1,2) (track)
            1st value is the first data point in the corresponding indexed data array
            2nd value is the first post calibration point in the corresponding indexed data array
        data_indicies: list of arrays of shape (1,2,3,4) (on/off)
            1st value is the first on data point in the corresponding indexed data array
            2nd value is the first transitioning point in the corresponding indexed data array
            3rd value is the first off data point in the corresponding indexed data array
            4th value is the first post calibration point in the corresponding indexed data array
        gain_start: list of arrays of shape (1,2)
            1st value is the gain delta for the given index
            2nd value is the time for the gain delta
        gain_end: list of arrays of shape (1,2)
            1st value is the gain delta for the given index
            2nd value is the time for the gain delta
        continuum:
        '''

        try:
            with fits.open(file_path) as hdul:
                self.header = hdul[0].header
                self.data = Table(hdul[1].data)
            self.file.logger.info(f"Error reading FITS file: {e}")
        except Exception as e:
            self.file.logger.error(f"Error reading FITS file: {e}")
            raise MyException(f"Error reading FITS file: {e}")
        
        self.file_path = file_path
        self.logger = Log_Collector(name=f"logger_{file_path}")

        self.validated_header = False
        self.validated_data = False

        self.gain_calibrated = []
        self.flux_calibrated = []

        self.labels = []
        self.data_indicies = []
        self.gain_start = []
        self.gain_end = []

        self.freqs = []

        self.continuum = []
        self.spectrum = []


if __name__ == "__main__":
    '''
    Test function to test uploaded file.
    '''

    file = Radio_File("C:/Users/starb/Downloads/0136645.fits")

    np.set_printoptions(threshold=100000)

    print(repr(file.header))
    print(repr(file.data['TAMBIENT']))
