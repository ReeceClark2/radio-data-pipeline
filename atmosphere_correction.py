import numpy as np
from astropy.io import fits
from astropy.table import Table
import itur
import utils


class Atmosphere_Correction:
    def __init__(self, file_path: str):
        '''
        Initialization function for provided file. Responsible for 
        opening the SDFITS file's header and data and initializing 
        necessary params.
        '''

        self.filepath = file_path

        with fits.open(self.filepath) as hdul:            
            # Use astropy's built in verification methods
            hdul.verify('exception')

            self.header = hdul[0].header
            self.data = Table(hdul[1].data)    
    
    def get_water_vapor_density(self, temperature, relative_humidity):
        '''
        Calculates water vapor density using the Buck equations for freezing and 
        non-freezing conditions.
        '''

        # If the temperature is above freezing then use the above freezing equation
        if temperature + 273.15 >= 0:
            # Buck equations require Celsius instead of Kelvin!
            e_s = (1.0007 + (3.46e-6)) * 6.1121 * np.exp((17.502 * (temperature - 273.15)) / ((temperature - 273.15) + 240.97)) 
        # Similary, if the temperature is below freezing then use the below freezing equation
        elif temperature + 273.15 < 0:
            e_s = (1.0003 + (4.18e-6)) * 6.1115 * np.exp((22.452 * (temperature - 273.15)) / ((temperature - 273.15) + 272.55))
        # If this fails then correct for oxygen attenuation lines only
        else:
            return 0
        
        e = (relative_humidity / 100) * e_s
        water_vapor_density = (216.7 * e) / temperature

        return water_vapor_density

    def gaseous_attenuation_correction(self, frequencies, elevation, water_vapor_density, pressure, temperature):
        '''
        Calculates the transmissions to correct signals for water vapor and 
        oxygen attenuation lines. This uses the International Telecommunication Union's 
        Radiocommunication standard library.
        '''

        g = itur.models.itu676.gaseous_attenuation_slant_path(frequencies, elevation, water_vapor_density, pressure, temperature, V_t=None, h=None, mode='approx')
        transmission = 10 ** (-(g.value) / 10.0)

        return transmission

    def atmosphere_correction(self):
        '''
        Loop through each spectrum and apply the atmosphere correction. 
        Weather parameters and elevation change during an observation, so 
        they must be corrected for time dependent.
        '''

        for i in self.data:
            frequencies = utils.get_frequency_range(self.header, i["IFNUM"])
            frequencies = np.linspace(frequencies[1] / 1000, frequencies[0] / 1000, frequencies[2])

            # Pull relevant parameters
            elevation = i["ELEVATIO"]
            temperature = i["TAMBIENT"] + 273.15 # SDFITS provide ambient temperature, so it is converted to Kelvin
            pressure = i["PRESSURE"]
            relative_humidity = i["HUMIDITY"]
            water_vapor_density = self.get_water_vapor_density(temperature, relative_humidity)

            gaseous_transmission = self.gaseous_attenuation_correction(frequencies, elevation, water_vapor_density, pressure, temperature)

            # This has been structured to accomodate additional transmission functions like 
            # cloud_attenuation if the appropriate instruments are installed. The transmissions 
            # would the be multiplicative.

            # Inversely apply the transmission to model initial signal.
            i["DATA"] *= (1 / gaseous_transmission)

        utils.save(self.filepath, self.header, self.data, "corrected")


if __name__ == "__main__":
    filepath = "C:/Users/starb/Downloads/0144767_validated.fits"
    
    ac = Atmosphere_Correction(filepath)
    ac.atmosphere_correction()
    