Skynet 2: Radio Data Processing Pipeline

This repository serves as the data reduction process for SDFITS files provided to the Skynet Robotic Telescope Network. The pipeline is designed to handle radio tracking, ON/OFF, and mapping observations. While integrated into Skynet, access to regularly scheduled calibration observations allows for all returned data to be flux calibrated into Janskys (Jy).

The pipeline supports radio astronomy for students, educators, and astronomers through programs like ERIRA, MWU!, and OPIS!.

Pipeline Overview

The processes required to properly handle radio observation data follow a generalized architecture designed to accommodate diverse instrumentation across a global network.

1. Validation

Validation is the first stage of the pipeline. This stage ensures that data is physical and maintains the integrity of the file.

    Formatting: Employs a Python-based pipeline to check data types (e.g., strings, floats).

Physicality: Checks for nonphysical values, such as negative temperatures in Kelvin.

Metadata: Employs an algorithm to identify calibration spikes using the CALSTATE (diode status) and SWPVALID (data validity) flags.

2. Atmosphere Correction

Skynet employs radio telescopes across the globe with L-, S-, C-, X-, Ku-, and K-band receivers. While lower frequency bands are not dramatically attenuated, higher bands suffer significant transmission loss due to water vapor and oxygen absorption lines.

    Vapor Pressure: Uses Buck equations to estimate saturation vapor pressure for above and below-freezing conditions.

Attenuation: Utilizes the ITU-R library to calculate power loss in decibels based on water vapor density and frequency.

Correction: The transmission fraction is inversely applied to observed intensities to find unattenuated values.

3. Gain Calibration

Radio receivers initially record data in arbitrary, nonphysical units.

    Noise Diode: Uses artificial noise diode spikes located in the receiver to establish a conversion to "calibration units".

Precision: Linear fits are applied to the diode's "on" and "off" data, with outliers rejected using Robust Chauvenet Rejection (RCR).

Drift Handling: If a significant drift is detected between pre- and post-calibration spikes (determined by a convolved Gaussian z-score > 1.96), a time-dependent interpolation is applied.

4. Flux Calibration (Continuum)

Data in calibration units is converted to physical Janskys (Jy) using observations of standard radio sources: Virgo A, Taurus A, and Cygnus A.

    Scheduling: A calibration source is observed every two hours to account for time-dependent diode responses.

Photometry: Sources are automatically photometered using Radio Cartographer (RC) to provide precise measurements.

Database: Secondary conversion factors are cataloged and interpolated to convert any radio continuum into Jy.

5. Spectrum & Observation Modes

The pipeline generates a spectrum by integrating the data cube along the time axis.

    Tracking: Provides time-dependent flux-calibrated intensities for pulsar timing or intensity measurements.

ON/OFF: Subtracts background sky coordinates to remove the receiver's frequency response and recover true source intensity.

Mapping: Supports raster and daisy maps, utilizing Radio Cartographer to create images that can be further photometered.

Repository Structure

    main.py: The entry point to run the entire pipeline for a provided SDFITS file.

    file_corruption.py: Testing utility used to manually corrupt files to verify validation effectiveness.

file_merge.py: Utility for data management and testing.

rmac_draft.pdf: Technical proceeding detailing the implemented methods.

Observer Controls

Processed observations provide a flux-calibrated continuum and a non-flux-calibrated spectrum. Observers have the ability to:

    RFI Mitigation: Manually cut Radio Frequency Interference portions; the pipeline then regenerates the continuum and spectrum.

Atmospheric Science: Opt-out of atmosphere correction if the goal is to study the atmosphere rather than astronomical sources.
