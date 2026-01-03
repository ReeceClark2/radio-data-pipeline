This repository serves as the data reduction process for SDFITS files provided to the Skynet Robotic Telescope Network. The pipeline is designed to handle radio tracking, ON/OFF, and (soon!) mapping files. While integrated into Skynet, access to regularly scheduled calibration observations allow for all returned data to be flux calibrated into Jy. The pipeline supports radio astronomy for students, educators, and astronomers!

Validation
Validation is the first stage of the pipeline. This stage ensures that data is physical and it ensures the integrity of the file. Validation also makes sure the time stamps are valid for the associated data and that no nan values exist.

Atmosphere Correction
Skynet will employ radio telescopes across the globe with L-, S-, C-, X-, Ku-, and K-band receivers. L-, S-, and C- bands are not dramatically attenuated by the atmosphere, but X-, Ku-, and K-bands begin to suffer signifcant transmission loss due to water vapor and oxygen absorption lines with K-band containing a water vapor line. The model uses Buck equations for approximating saturation vapor pressure and ITUR's suggestion for slant path attenuations.

Continuum
Continuum creates a gain calibrated and flux calibrated observation using pre- and post-calibration spikes in the observation to place the observation into calibration units and a calibration observation to convert calibration units to Jy. The pipeline ensures a robust data reduction process for the continuum.

Spectrum
Spectrum creates a spectrum for the expected frequency range while also handling ON/OFF files when relevant.

Additional files, file_corruption.py and file_merge.py, are used for testing with main.py used to run the entire pipeline from start to end for a provided SDFITS file.
