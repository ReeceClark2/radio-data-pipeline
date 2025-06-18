from file_init import Mike
from val import Val
from sort import Sort
from weather import Weather
from gain_calibration import Gain_Cal
from flux_calibration import Flux_Cal
from spectrum import Spectrum
from child_init import Sully

import numpy as np
import matplotlib
import os
import time
matplotlib.use("Agg")
import matplotlib.pyplot as plt


if __name__ == "__main__":
    fits_folder = "EpicONoFFHiRes"
    fits_files = [f for f in os.listdir(fits_folder) if f.endswith(".fits")]

    for fits_filename in fits_files:
        file_path = os.path.join(fits_folder, fits_filename)

        file = Mike(file_path)
        cal_file = Mike("TrackingHighRes/0136483.fits")
            
        v = Val(file)
        v.validate_primary_header()
        v.validate_data()

        s = Sort(file)
        s.sort()

        # w = Weather(file)
        # w.weather_correction()

        gc = Gain_Cal(file)
        gc.gain_cal()


        v = Val(cal_file)
        v.validate_primary_header()
        v.validate_data()

        s = Sort(cal_file)
        s.sort()

        w = Weather(cal_file)
        w.weather_correction()

        gc = Gain_Cal(cal_file)
        gc.gain_cal()

        # fc = Flux_Cal(file, cal_file)
        # fc.flux_cal()

        spec = Spectrum(file)
        spec.make_spec()

        start_time = time.time()
        print(f"Processing file: {fits_filename}")
        keep_indices = [[7,13], [15,21]]  # Specify the indices you want to keep
        feed= [0,1]  # Specify the feeds you want to keep
        if keep_indices != []:
            child = Sully(file, file_path, keep_indices, "continuum", "cut", feed)
            child.user_cuts([[1300, 1400], [1402, 1404], [1412, 1420]], "spectrum", "keep", feed)
            child.specMaker()
            gcChild = Gain_Cal(child)
            gcChild.gain_cal()

            # print(file.data[0]["DATA"].shape)
            # print(child.data[0]["DATA"].shape)

        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Time taken for this file: {elapsed_time:.4f} seconds")

        # Store elapsed times for averaging
        if 'elapsed_times' not in locals():
            elapsed_times = []
        elapsed_times.append(elapsed_time)
    print(f"Average time taken so far: {np.mean(elapsed_times):.4f} seconds")


    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    file_path = "TrackingHighRes/0136764.fits"
    file = Mike(file_path)
    cal_file = Mike("TrackingHighRes/0136483.fits")


    v = Val(file)
    v.validate_primary_header()
    v.validate_data()

    s = Sort(file)
    s.sort()
    # w = Weather(file)
    # w.weather_correction()

    gc = Gain_Cal(file)
    gc.gain_cal()

    v = Val(cal_file)
    v.validate_primary_header()
    v.validate_data()

    s = Sort(cal_file)
    s.sort()

    w = Weather(cal_file)
    w.weather_correction()

    gc = Gain_Cal(cal_file)
    gc.gain_cal()

    # fc = Flux_Cal(file, cal_file)
    # fc.flux_cal()

    spec = Spectrum(file)
    spec.make_spec()

    keep_indices = [[7,13], [15,21]]  # Specify the indices you want to keep
    feed= [0,1]  # Specify the feeds you want to keep
    if keep_indices != []:
        child = Sully(file, file_path, keep_indices, "continuum", "cut", feed)
        child.user_cuts([[1300, 1400], [1402, 1404], [1412, 1420]], "spectrum", "keep", feed)
        child.specMaker()
        print('hi')
        gcChild = Gain_Cal(child)
        gcChild.gain_cal()

        print(file.data[0]["DATA"].shape)
        print(child.data[0]["DATA"].shape)


    fig, axs = plt.subplots(2, 2, figsize=(12, 8))

    # Plot first two spectra (top-left)
    for i in range(2):
        x, y = child.spectrum[i]
        axs[0, 0].plot(x, y, label=child.labels[i], linewidth=0.5)
    axs[0, 0].set_xlabel('Frequency (MHz)')
    axs[0, 0].set_ylabel('Intensity')
    axs[0, 0].set_title('Spectra 1 & 2')
    axs[0, 0].legend()

    # Plot last two spectra (bottom-left)
    for i in range(2, 4):
        x, y = child.spectrum[i]
        axs[1, 0].plot(x, y, label=child.labels[i], linewidth=0.5)
    axs[1, 0].set_xlabel('Frequency (MHz)')
    axs[1, 0].set_ylabel('Intensity')
    axs[1, 0].set_title('Spectra 3 & 4')
    axs[1, 0].legend()

    # Plot first two continuum (top-right)
    for i in range(2):
        x, y = child.continuum[i]
        axs[0, 1].plot(x, y, label=child.labels[i], linewidth=0.5)
    axs[0, 1].set_xlabel('Time (Seconds)')
    axs[0, 1].set_ylabel('Intensity (Jy)')
    axs[0, 1].set_title('Continuum 1 & 2')
    axs[0, 1].legend()

    # Plot last two continuum (bottom-right)
    for i in range(2, 4):
        x, y = child.continuum[i]
        axs[1, 1].plot(x, y, label=child.labels[i], linewidth=0.5)
    axs[1, 1].set_xlabel('Time (Seconds)')
    axs[1, 1].set_ylabel('Intensity (Jy)')
    axs[1, 1].set_title('Continuum 3 & 4')
    axs[1, 1].legend()

    plt.tight_layout()
    plt.savefig('GB Ours2', dpi=400)
    print("Plot saved as 'GB Ours2.png'")


