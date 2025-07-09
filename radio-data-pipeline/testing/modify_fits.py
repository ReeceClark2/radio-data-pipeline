from pathlib import Path
from astropy.io import fits
import numpy as np


def update_obsfreq_in_folder(input_dir: str, output_dir: str, new_obsfreq: float):
    """
    Updates the OBSFREQ header value in all FITS files in a folder and saves them to a new directory.

    Args:
        input_dir (str): Absolute path to the folder containing input FITS files.
        output_dir (str): Absolute path to the folder where modified FITS files will be saved.
        new_obsfreq (float): New value to assign to the OBSFREQ header keyword.
    """
    input_path = Path(input_dir).resolve()
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)

    fits_files = list(input_path.glob("*.fits"))
    if not fits_files:
        print(f"[INFO] No FITS files found in {input_path}")
        return

    for file in fits_files:
        with fits.open(file) as hdul:
            print('BEFORE:')
            print('  OBSFREQ (header):', type(hdul[0].header["OBSFREQ"]), hdul[0].header["OBSFREQ"])
            print('  BMAJ:', type(hdul[0].header["BMAJ"]), hdul[0].header["BMAJ"])
            print('  BMIN:', type(hdul[0].header["BMIN"]), hdul[0].header["BMIN"])
            print('  CRVAL1 (data):', type(hdul[1].data['CRVAL1']), hdul[1].data['CRVAL1'])
            print('  OBSFREQ (data):', type(hdul[1].data['OBSFREQ']), hdul[1].data['OBSFREQ'])
            print('  CRVAL1 shape:', hdul[1].data['CRVAL1'].shape)

            hdul[0].header["OBSFREQ"] = float(new_obsfreq)
            hdul[0].header["BMAJ"] = float(np.rad2deg(1.22 * ((3e8 / (new_obsfreq * 1e6)) / 20)))
            hdul[0].header["BMIN"] = float(np.rad2deg(1.22 * ((3e8 / (new_obsfreq * 1e6)) / 20)))

            # This loop is unnecessary unless you're doing something per element
            hdul[1].data['CRVAL1'][:] = new_obsfreq * 1e6
            hdul[1].data['OBSFREQ'][:] = new_obsfreq * 1e6

            output_file = output_path / file.name
            hdul.writeto(output_file, overwrite=True)

            print('AFTER:')
            print('  OBSFREQ (header):', type(hdul[0].header["OBSFREQ"]), hdul[0].header["OBSFREQ"])
            print('  BMAJ:', type(hdul[0].header["BMAJ"]), hdul[0].header["BMAJ"])
            print('  BMIN:', type(hdul[0].header["BMIN"]), hdul[0].header["BMIN"])
            print('  CRVAL1 (data):', type(hdul[1].data['CRVAL1']), hdul[1].data['CRVAL1'])
            print('  OBSFREQ (data):', type(hdul[1].data['OBSFREQ']), hdul[1].data['OBSFREQ'])
            print('  CRVAL1 shape:', hdul[1].data['CRVAL1'].shape)

            print(f"[INFO] Updated {file.name} and saved to {output_file}")


update_obsfreq_in_folder("C:/Users/starb/OneDrive/Desktop/SUMMER-2025/mike/radio-data-pipeline/data/cyg_a_low_hI", "C:/Users/starb/OneDrive/Desktop/SUMMER-2025/mike/radio-data-pipeline/data/cyg_a_low_HI_mod", 1395)
