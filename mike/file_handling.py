# Third-party libraries
import numpy as np

# Local application imports
from data_processing.file_init import Radio_File
from data_processing.flux_calibration import Flux_Cal
from data_processing.gain_calibration import Gain_Cal
from data_processing.sort import Sort
from data_processing.spectrum import Spectrum
from data_processing.val import Val
from data_processing.weather import Weather


file = Radio_File("C:/Users/starb/Downloads/0136379.fits")

v = Val(file)
v.validate_primary_header()
v.validate_data()


class Track:
    def __init__(self, file):
        self.file = file

    def process(self):
        s = Sort(file) # Needs to know if on/off
        s.sort()

        w = Weather(file)
        w.weather_correction()

        gc = Gain_Cal(file)
        gc.gain_cal()

        fc = Flux_Cal(file)
        fc.flux_cal()


class OnOff:
    def __init__(self, file):
        self.file = file

    def process(self):
        s = Sort(file) # Needs to know if on/off
        s.sort()

        w = Weather(file)
        w.weather_correction()

        gc = Gain_Cal(file)
        gc.gain_cal()

        fc = Flux_Cal(file)
        fc.flux_cal()


class Map:
    def __init__(self, file):
        self.file = file
        self.cal = False

        if self.file.file_type == 'cal':
            self.cal = True

    def process(self):
        s = Sort(file) # Needs to know if on/off
        s.sort()

        w = Weather(file)
        w.weather_correction()

        gc = Gain_Cal(file)
        gc.gain_cal()

        if self.cal == True:
            pass
            # TODO Handle flux calibration for calibration files here
        
        else:
            pass
            # TODO Write normal process
        


if file.file_type == 'tracking':
    Track(file)
elif file.file_type == 'onoff':
    OnOff(file)
elif file.file_type == 'map' or file.file_type == 'cal':
    Map(file)