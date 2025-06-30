# Standard library
import os
import shutil
import subprocess
from pathlib import Path

# Third-party libraries
import numpy as np
from astropy.io import fits

# Local application imports
# (Assumes this file lives in the correct package relative to others like Gain_Cal)
# e.g., from .some_module import some_function  # if needed


def run_rc(input_file: str,
           container_name: str = "rc",
           script_path: str = "/skynet/radio-cartographer/radio-cartographer/run2.sh",
           parameters: dict = None,
           output_dir: Path = None):
    """
    Run the run2.sh script inside the given Docker container, passing parameters as environment variables.
    After execution, copy the _phot results folder from the container to the host output directory.

    Args:
        input_file (str): Path to the input FITS file inside the container.
        container_name (str): Name of the running Docker container.
        script_path (str): Absolute path to the run2.sh script in the container.
        parameters (dict): Dictionary of environment variables to pass to the container.
        output_dir (Path): Directory on the host where results will be copied.
    """
    env_vars = {"FILES_TO_PROCESS": input_file}
    if parameters:
        env_vars.update(parameters)

    input_basename = Path(input_file).stem
    container_output_dir = f"/skynet/radio-cartographer/radio-cartographer/{input_basename}_phot"

    cmd = ["docker", "exec"]
    for key, value in env_vars.items():
        cmd += ["-e", f"{key}={value}"]
    cmd += [container_name, "bash", script_path]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Docker execution failed:\n{e.stderr}")
        return

    # Define where to copy results to locally
    if output_dir is None:
        output_dir = Path(__file__).resolve().parent.parent / "results"
    output_dir.mkdir(parents=True, exist_ok=True)

    local_output = output_dir / f"{input_basename}_phot"
    if local_output.exists():
        shutil.rmtree(local_output)

    try:
        subprocess.run([
            "docker", "cp",
            f"{container_name}:{container_output_dir}",
            str(local_output)
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to copy results from container:\n{e.stderr}")
        return

    # Convert .txt output to FITS
    create_fits_from_txt_folder(local_output, local_output / f"{input_basename}_output.fits")


def get_data(filepath: Path) -> np.ndarray:
    """
    Load a .txt file as a NumPy array, then delete the file.

    Args:
        filepath (Path): Path to the .txt file.

    Returns:
        np.ndarray: Loaded data array.
    """
    data = np.genfromtxt(filepath, delimiter='\t')
    filepath.unlink(missing_ok=True)
    return data


def create_fits_from_txt_folder(input_dir: Path, output_path: Path):
    """
    Create a FITS file from a set of .txt data files in a directory.

    Args:
        input_dir (Path): Folder containing the .txt files.
        output_path (Path): Path to the output FITS file.
    """
    if not input_dir.exists():
        raise FileNotFoundError(f"[ERROR] Input directory does not exist: {input_dir}")

    files = ['main', 'path', 'scale', 'weight', 'correlation', 'raw']
    hdu_list = fits.HDUList()

    for name in files:
        txt_file = input_dir / f"{name}.txt"
        if not txt_file.exists():
            continue

        data = get_data(txt_file)
        hdu = fits.PrimaryHDU(data) if name == "main" else fits.ImageHDU(data)
        hdu_list.append(hdu)

    hdu_list.writeto(output_path, overwrite=True)
    print(f"[INFO] FITS file written to: {output_path.resolve()}")



if __name__ == "__main__":
    file = "0137018"
    local_result_dir = Path(__file__).resolve().parent.parent / "results" / f"{file}_phot"

    # Clean up previous results if they exist
    if local_result_dir.exists() and local_result_dir.is_dir():
        print(f"[INFO] Removing existing result directory: {local_result_dir}")
        shutil.rmtree(local_result_dir)

    # Option A: Run inside container and extract results
    if True:
        run_rc(
            input_file=f"/skynet/radio-cartographer/radio-cartographer/testing/{file}.fits",
            parameters={
                "COORD": "equatorial",
                "PROCCOORD": "equatorial",
                "MINFREQ": "1355.0",
                "MAXFREQ": "1435.0",
            }
        )

    # Option B: Just convert .txt output to FITS
    else:
        output_fits = local_result_dir / f"{file}_output.fits"
        create_fits_from_txt_folder(local_result_dir, output_fits)

