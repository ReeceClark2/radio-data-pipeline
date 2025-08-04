# Radio Data Pipeline

This repository provides a data processing pipeline for radio astronomy FITS files. It's designed to prepare raw observation data for scientific analysis by performing a series of validation and calibration steps. The core of the repository consists of a **pipeline** for data validation and a **workflow** for processing different types of observation files.

---

## Pipeline

The pipeline is responsible for the initial processing and cleaning of the raw FITS files. The process ensures data quality and prepares the files for subsequent analysis.

### `validate.py`

This script handles the validation of FITS files, correcting common issues and extracting key information. The validation steps include:

* **Header Validation**: Confirms that the FITS header adheres to the standard format, checking for correct size, required keywords like `SIMPLE` and `BITPIX`, and ensuring no duplicate or empty cards.
* **Data Validation**: Ensures the integrity of the data table. This involves checking for `NaN` values and ensuring consistent data types within each column. Columns containing time-related data are also converted to a standard format.
* **Calibration Indexing**: Automatically identifies and marks the start and end indices for pre- and post-calibration sections within the data. It also finds the transition point for "on/off" observations, which is crucial for background subtraction.
* **Channel Cleaning**: Removes invalid or corrupted data channels to ensure only reliable data is used.

The `validate.py` script saves the processed file with a `_validate` suffix, containing the cleaned data and updated header with the newly identified calibration indices.

### `weather.py`

This is a planned future installment that will perform weather attenuation corrections on the FITS files, further enhancing the quality of the data.

---

## Workflow

The workflow directory contains scripts tailored to different observation types. These scripts use the validated FITS files to produce calibrated scientific data products like spectra and continuums.

### `radio_tracking_continuum.py`

This script processes radio tracking files to produce a continuum plot. A continuum observation measures the average radio emission across a broad frequency range over time. The script performs the following steps:

* **Data Filtering**: Filters the observation data based on specified time and frequency ranges, allowing users to isolate specific parts of the observation for analysis.
* **Gain Calibration**: Calculates and applies gain corrections using the calibration sections identified during the validation stage. This corrects for fluctuations in the receiver's sensitivity over the course of the observation. The script can use a single calibration spike or interpolate between pre- and post-calibration spikes if both are available and stable.
* **Continuum Generation**: Computes a single intensity value for each time step by averaging the channel data, resulting in a continuum data product.
* **Flux Calibration**: Applies a conversion factor to the continuum data to convert arbitrary intensity units into meaningful flux units (e.g., Janskys).

### `radio_tracking_spectrum.py`

This script processes radio tracking files to produce a spectrum. A spectrum shows the intensity of radio emission as a function of frequency. The script performs the following steps:

* **Data Filtering**: Similar to the continuum script, this filters the observation data based on specified time and frequency ranges.
* **Spectrum Generation**: Averages the intensity across the time axis for each frequency channel. This creates a spectrum that represents the average emission profile of the observed target.

---

### Other Workflow Files

Future workflow scripts will be developed to handle other observation modes:

* **On/Off**: These files contain data from a target ("on") and a nearby blank sky region ("off"). The workflow for these files will involve subtracting the "off" data from the "on" data to isolate the target's signal from background noise.
* **Map**: This will process data from raster or daisy-style mapping observations to create two-dimensional images of the observed sky region.
* **Calibration**: These files are specific observations of a known radio source used to determine the receiver's conversion factor. The workflow for these files will automate this process and produce a reliable factor for flux calibration.
