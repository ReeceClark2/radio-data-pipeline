This repository serves as the data reduction process for SDFITS files provided to the **Skynet Robotic Telescope Network**. The pipeline is designed to handle **radio tracking**, **ON/OFF**, and **mapping** observations. While integrated into Skynet, access to regularly scheduled calibration observations allows for all returned data to be flux calibrated into **Janskys (Jy)**. The pipeline supports radio astronomy for students, educators, and astronomers! 

---

## Data Processing Pipeline

The processes required to properly handle radio observation data follow a generalized architecture designed to accommodate diverse instrumentation across a global network.

### Validation

Validation is the first stage of the pipeline. This stage ensures that data is physical and maintains the integrity of the file. The pipeline employs a **Python-based system** that checks for proper formatting (e.g., strings vs. floats) and flags nonphysical values, such as negative temperatures in Kelvin. It also ensures that calibration diode flags (**CALSTATE**) and data collection flags (**SWPVALID**) are present to accurately identify calibration spikes.

### Atmosphere Correction

Skynet employs radio telescopes across the globe with **L-, S-, C-, X-, Ku-, and K-band** receivers. While lower frequency bands are not dramatically attenuated, **X-, Ku-, and K-bands** suffer significant transmission loss due to water vapor and oxygen absorption lines.

* 
**Modeling**: The pipeline uses **Buck equations** to estimate saturation vapor pressure based on embedded weather data.


* 
**Attenuation**: It utilizes the **ITU-R** Python library to calculate power loss (dB) based on water vapor density and frequency, converting this to a transmission fraction inversely applied to observed intensities.


* 
**User Control**: Observers performing atmospheric science can opt to skip this correction to preserve atmospheric data.



### Gain Calibration

Radio receivers initially record data in arbitrary, nonphysical units.

* 
**Calibration Spikes**: The pipeline identifies artificial noise diode spikes triggered before and after observations.


* 
**Precision**: Using **Robust Chauvenet Rejection (RCR)**, the system performs linear fits on diode data to calculate a conversion factor into "calibration units".


* 
**Drift Handling**: If a significant drift is detected between pre- and post-calibration spikes (determined by a **z-score threshold** of 1.96), the pipeline applies a time-dependent interpolation to the data.



### Flux Calibration (Continuum)

To convert calibration units into physical **Janskys (Jy)**, Skynet observes bright radio sources—**Virgo A, Taurus A, and Cygnus A**—every two hours.

* 
**Automated Photometry**: These sources are photometered using **Radio Cartographer (RC)** to measure calibration units precisely.


* 
**Final Reduction**: The pipeline interpolates between these secondary conversion factors to provide a science-ready flux calibrated continuum.



### Spectrum

The pipeline generates a spectrum by integrating the data cube along the time axis.

* 
**ON/OFF Support**: For ON/OFF observations, the pipeline subtracts the background sky to remove the receiver's frequency response and recover the true intensity of the source.


* 
**Applications**: This provides clean data for measuring **redshifts** or **HI broadening**. While currently returned in arbitrary units for precision frequency work, bandpass calibration methods for the spectrum are in development.


---

## File Structure

* `main.py`: The entry point to run the entire pipeline from start to end for a provided SDFITS file.
* 
`file_corruption.py`: Testing utility used to manually corrupt files (e.g., nonphysical values) to verify validation effectiveness.


* `file_merge.py`: Utility for data management and testing.
* 
`radio_cartographer`: Integration with RC for generating radio images and performing photometry.



Would you like me to draft a section for the README on how to use the `Observer Controls` for RFI mitigation?
