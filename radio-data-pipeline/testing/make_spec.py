import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def average_intensity(data, axis):
    """
    Average across time (axis=0) or integrate across frequency (axis=1).
    """
    intensities = np.array([row[6] for row in data])
    return np.sum(intensities, axis=axis) / intensities.shape[axis]


def sdfits_to_spectrum(data, freq_start=1414, freq_end=1430):
    """
    Convert filtered SDFITS data to frequency and intensity arrays.
    """
    intensities = average_intensity(data, axis=0)
    freqs = np.linspace(freq_start, freq_end, len(intensities))
    return freqs, intensities


def plot_calibrated_spectrum(ax, freqs, intensities):
    ax.plot(freqs, intensities[::-1], linewidth=0.6)
    ax.set_ylim(0, 1)
    ax.set_title("Calibrated Spectrum")
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("Intensity (arb. units)")


def plot_gain_deltas(ax, df):
    pre = df['PreCalDelta']
    post = df['PostCalDelta']
    diff = np.abs(pre - post)

    ax.plot(pre, label="Pre Cal", color='tab:blue', alpha=0.6, linewidth=0.4)
    ax.plot(post, label="Post Cal", color='tab:orange', alpha=0.6, linewidth=0.4)
    ax.plot(diff, label="Difference", color='tab:red', linewidth=0.4)

    ax.invert_xaxis()
    ax.set_title("Gain Deltas")
    ax.set_xlabel("Channel")
    ax.set_ylabel("Delta")
    ax.legend()


def plot_dual_overlay(ax, freqs, intensities, df_diff):
    # Interpolate gain difference to match frequency resolution
    interp_func = interp1d(
        np.linspace(0, len(df_diff) - 1, len(df_diff)),
        df_diff,
        kind='linear',
        bounds_error=False,
        fill_value="extrapolate"
    )
    resampled_diff = interp_func(np.linspace(0, len(df_diff) - 1, len(freqs)))

    ax_left = ax
    ax_right = ax_left.twinx()

    ax_left.plot(freqs, intensities[::-1], color='black', linewidth=0.6, label="Spectrum")
    ax_right.plot(freqs, resampled_diff[::-1], color='tab:red', linewidth=0.6, alpha=0.8, label="ΔGain")

    ax_left.set_title("Overlay: Spectrum vs. ΔGain")
    ax_left.set_xlabel("Frequency (MHz)")
    ax_left.set_ylabel("Intensity", color='black')
    ax_right.set_ylabel("Gain Δ", color='tab:red')

    ax_left.tick_params(axis='y', labelcolor='black')
    ax_right.tick_params(axis='y', labelcolor='tab:red')


def main():
    # Filepaths
    fits_path = "C:/Users/starb/Downloads/0136376_gain_calibrated.fits"
    csv_path = "C:/Users/starb/OneDrive/Desktop/SUMMER-2025/mike/radio-data-pipeline/results/gain_deltas.csv"

    # Load FITS data
    with fits.open(fits_path) as hdul:
        data_table = Table(hdul[1].data)

    # Filter IFNUM=0 and PLNUM=0
    filtered_data = data_table[(data_table['IFNUM'] == 0) & (data_table['PLNUM'] == 0)]
    freqs, intensities = sdfits_to_spectrum(filtered_data)

    # Load gain deltas CSV
    df = pd.read_csv(csv_path)
    gain_diff = np.abs(df['PreCalDelta'] - df['PostCalDelta'])

    # Create plots
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))

    plot_calibrated_spectrum(axs[0], freqs, intensities)
    plot_gain_deltas(axs[1], df)
    plot_dual_overlay(axs[2], freqs, intensities, gain_diff)

    plt.tight_layout()
    plt.savefig("calibrated_spectrum_with_gain_deltas.png", dpi=400)


if __name__ == "__main__":
    main()
