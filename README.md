# Skynet 2: Radio Data Processing Pipeline

This repository serves as the data reduction process for SDFITS files provided to the **Skynet Robotic Telescope Network**. The pipeline is designed to handle **radio tracking**, **ON/OFF**, and **mapping** observations. While integrated into Skynet, access to regularly scheduled calibration observations allows for all returned data to be flux calibrated into **Janskys (Jy)**.

The pipeline supports radio astronomy for students, educators, and astronomers through programs like ERIRA, MWU!, and OPIS!.

---

## Data Processing Pipeline

The processes required to properly handle radio observation data follow a generalized architecture designed to accommodate diverse instrumentation across a global network.

### 1. Validation

Validation is the first stage of the pipeline. This stage ensures that data is physical and maintains the integrity of the file.

**Formatting**: The pipeline employs a Python-based system to check data types for proper formatting, such as strings and floats.

**Physicality**: It ensures the absence of nonphysical values, such as negative temperatures in Kelvin.

**Metadata**: The system identifies calibration spikes by checking the CALSTATE (diode status) and SWPVALID (data validity) flags, which are then recorded in the SDFITS header.


### 2. Atmosphere Correction

Skynet employs radio telescopes across the globe with L-, S-, C-, X-, Ku-, and K-band receivers. While lower frequency bands are not dramatically attenuated, higher bands like X-, Ku-, and K-band can experience significant transmission loss due to water vapor and oxygen absorption lines.

**Vapor Pressure**: The pipeline uses Buck equations to estimate saturation vapor pressure for both above and below-freezing conditions.

**Attenuation**: It utilizes the ITU-R library to calculate power loss in decibels based on water vapor density and target frequencies.

**Correction**: Modeled transmission rates are inversely applied to observed intensities to recover unattenuated values.


### 3. Gain Calibration

Radio receivers initially record data in arbitrary, nonphysical units.

**Noise Diode**: Gain calibration uses artificial noise diode spikes triggered before and after the observation to establish a conversion to "calibration units".

**Precision**: Linear fits are applied to the diode data, and outliers are rejected using Robust Chauvenet Rejection (RCR) to increase the precision of the conversion factor.

**Drift Handling**: If a significant drift is detected between pre- and post-calibration spikes—determined by a convolved Gaussian z-score exceeding a 1.96 threshold—a time-dependent interpolation is applied.


### 4. Flux Calibration

Data in calibration units is converted to physical Janskys (Jy) using observations of standard radio sources, including Virgo A, Taurus A, and Cygnus A.

**Scheduling**: A calibration source is observed every two hours to account for time-dependent diode responses.

**Photometry**: Sources are automatically photometered using Radio Cartographer (RC) to provide precise measurements.

**Database**: Secondary conversion factors are cataloged in the Skynet 2 database and interpolated to convert the observation into Jy.


### 5. Observation Modes

The pipeline supports multiple observation modes to meet diverse scientific and educational objectives.

**Radio Tracking**: Allows for time-dependent flux-calibrated intensities to be recorded for measuring pulsar periods or intensities.

**ON/OFF**: Subtracts background sky coordinates to remove the receiver's frequency response and recover true source intensity for measuring redshifts or HI broadening.

**Mapping**: Supports raster and daisy maps, utilizing Radio Cartographer to create images that can be photometered to perform radio astronomy.


## Helper Files

* `main.py`: The entry point to run the entire pipeline from start to end for a provided SDFITS file.
* `file_corruption.py`: Testing utility used to manually corrupt files to verify the effectiveness of the validation pipeline.
* `file_merge.py`: Utility for data management and testing.


## Observer Controls

Observers are provided with tools to ensure clean, calibrated data is returned.

**RFI Mitigation**: If Radio Frequency Interference is present, observers can cut these portions, and the pipeline will regenerate the continuum and spectrum to reflect the changes.

**Atmospheric Science**: Observers performing atmospheric science can choose to opt-out of the atmosphere correction process.
