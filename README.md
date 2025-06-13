# Radio Data Pipeline
This repository provides Skynet 2 with the ability to process radio FITS files. The FITS files are wrapped in a file class with several fields to control the calibration of the file. This pipeline can take in four different file types: tracking, on/off, map, and calibration. All files are treated with validation, channel sorting, weather correction, gain calibration, and flux calibration.

## Tracking
Radio tracking files are files where the observation is fixed onto a set of coordinates on the World Coordinate System (WCS). The data recorded includes the intensity across a frequency range for some time step that continues for the duration of the observation. This data is used to create spectra and continuum results.

## On/Off
On/off files have the first part of the observation focused onto some target with the second part of the observation removed from the target. By using the on and off parts of the data, the background can be subtracted by on - off.

## Map
Map files are created using Radio Cartographer. Maps can be created using raster maps or daisy maps. Raster maps scan the selected region by moving across the whole width of the map, changing height slightly, and repeating, effectively mapping the sky. Daisies follow a rhodonea curve with some number of leaves.

## Calibration
These files are taken every 2 hours. The files are automatically photometered with Radio Cartographer and used for flux calibration.
