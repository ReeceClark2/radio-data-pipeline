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
            hdul[0].header["OBSFREQ"] = new_obsfreq
            hdul[0].header["BMAJ"] = np.rad2deg(1.22 * ((3e8 / (new_obsfreq * 1e6)) / 20))
            hdul[0].header["BMIN"] = np.rad2deg(1.22 * ((3e8 / (new_obsfreq * 1e6)) / 20))
            output_file = output_path / file.name
            hdul.writeto(output_file, overwrite=True)
            print(f"[INFO] Updated {file.name} and saved to {output_file}")


update_obsfreq_in_folder("C:/Users/starb/OneDrive/Desktop/SUMMER-2025/mike/radio-data-pipeline/data/cyg_a_hi", "C:/Users/starb/OneDrive/Desktop/SUMMER-2025/mike/radio-data-pipeline/data/cyg_a_hi_mod", 1421.875)
