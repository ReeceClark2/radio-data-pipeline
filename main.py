import os
from time import time

import utils
from validate import Validation
from atmosphere_correction import Atmosphere_Correction
from continuum import Continuum
from spectrum import Spectrum

import matplotlib.pyplot as plt


if __name__ == "__main__":
    # These ranges are specified as a lists of lists such as [[start1, end1],[start2, end2]]
    including_frequency_ranges = None
    excluding_frequency_ranges = None 
    including_time_ranges = None
    excluding_time_ranges = None

    start_time = time()

    filepath = "C:/Users/starb/Downloads/Raw/0144717daisy_merge.fits"
    root, extension = os.path.splitext(filepath)

    v = Validation(filepath)
    v.validate()

    time0 = time()
    print("Validation time:", round((time0 - start_time),3), " seconds")

    # ac = Atmosphere_Correction(root + "_validated" + extension)
    # ac.atmosphere_correction()

    time1 = time()
    print("Atmosphere correction time:", round((time1 - time0),3), " seconds")

    c = Continuum(root + "_validated" + extension, 0, 1, including_frequency_ranges, excluding_frequency_ranges, including_time_ranges, excluding_time_ranges)
    continuum = c.continuum()

    time2 = time()
    print("Continuum creation time:", round((time2 - time1),3), " seconds")

    s = Spectrum(root + "_validated" + extension, 0, 1, including_frequency_ranges, excluding_frequency_ranges, including_time_ranges, excluding_time_ranges)
    spectrum = s.spectrum()

    time3 = time()
    print("Spectrum creation time:", round((time3 - time2),3), " seconds")

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))

    axes[0].plot(continuum[0], continuum[1], color="black")
    axes[0].set_xlim(min(continuum[0]), max(continuum[0]))
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Intensity")
    axes[0].set_title("Continuum")

    axes[1].plot(spectrum[0], spectrum[1], color="black")
    axes[1].set_xlim(min(spectrum[0]), max(spectrum[0]))
    axes[1].set_xlabel("Frequency (MHz)")
    axes[1].set_ylabel("Intensity")
    axes[1].set_title("Spectrum")

    plt.tight_layout()
    plt.show()
