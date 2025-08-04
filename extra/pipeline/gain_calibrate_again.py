# Third-party libraries
import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from scipy.stats import linregress
import rcr
import time


class GainCalibrator:
    """
    A class to perform gain calibration on telescope FITS data.

    Workflow:
    - Load FITS file and extract data table and header indices.
    - For each (IFNUM, PLNUM) channel combination:
      * Extract pre- and post-calibration ON/OFF data subsets.
      * Fit linear models using Robust Chauvenet Rejection.
      * Compute gain delta and uncertainty.
      * Scale the DATA column appropriately.
    """

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.data = None
        self.header = None

        # Header indices for different data regions
        self.data_start_index = None
        self.post_calibration_start_index = None
        self.off_start_index = None

        self._load_fits(filepath)

    def _load_fits(self, filepath: str):
        """Load FITS file, header, and data table."""
        try:
            with fits.open(filepath) as hdul:
                self.header = hdul[0].header
                self.data = Table(hdul[1].data)

                # Read important indices
                self.data_start_index = self.header.get("DATAIND")
                self.post_calibration_start_index = self.header.get("POSTIND")
                self.off_start_index = self.header.get("ONOFFIND")

        except Exception as e:
            raise RuntimeError(f"Failed to load FITS file {filepath}: {e}")

    # ----------------- Model and RCR utilities -----------------

    @staticmethod
    def _linear(x, params):
        """Simple linear model y = b + m*x."""
        return params[0] + x * params[1]

    @staticmethod
    def _d_linear_1(x, params):
        """Derivative of linear model wrt first parameter (intercept)."""
        return 1

    @staticmethod
    def _d_linear_2(x, params):
        """Derivative of linear model wrt second parameter (slope)."""
        return x

    def perform_rcr(self, x: np.ndarray, y: np.ndarray):
        """
        Perform functional Robust Chauvenet Rejection on a dataset.
        Returns:
            (params, uncertainties) -> best-fit [b, m] and their std devs.
        """
        # Center x values
        x = np.array(x) - np.average(x)
        y = np.array(y)

        # Initial linear regression guess
        result = linregress(x, y)
        guess = [result.intercept, result.slope]

        model = rcr.FunctionalForm(
            self._linear,
            x,
            y,
            [self._d_linear_1, self._d_linear_2],
            guess
        )

        r = rcr.RCR(rcr.SS_MEDIAN_DL)
        r.setParametricModel(model)
        r.performBulkRejection(y)

        indices = r.result.indices
        x = x[indices]
        y = y[indices]

        # Final fit parameters
        best_fit_parameters = model.result.parameters  # [b, m]

        # Uncertainties
        sigma = (1 / (len(x) - 2)) * np.sum((y - (best_fit_parameters[1] * x + best_fit_parameters[0])) ** 2)
        m_sd = np.sqrt(sigma / np.sum((x - np.mean(x)) ** 2))
        b_sd = np.sqrt(sigma * ((1 / len(x)) + (np.mean(x) ** 2 / np.sum((x - np.mean(x)) ** 2))))

        return best_fit_parameters, (b_sd, m_sd)

    # ----------------- Data extraction utilities -----------------

    def _get_time_and_mean(self, table_subset):
        """
        Given a subset of data rows, return (relative_times, mean_intensities).
        """
        intensities = np.array(table_subset['DATA'])
        if intensities.size == 0:
            raise ValueError("Empty data subset.")

        channel_means = np.sum(intensities, axis=1) / intensities.shape[1]

        times = Time(table_subset["DATE-OBS"], format='isot')
        t0 = Time(self.header["DATE"], format="isot")
        time_rel = (times - t0).sec

        return time_rel, channel_means

    def _mask_cal_data(self, start, ifnum, plnum, cal_state):
        """
        Return a Table subset filtered by CALSTATE, SWPVALID, IFNUM, PLNUM,
        starting from a given index.
        """
        subset = self.data[start:]
        mask = (
            (subset['CALSTATE'] == cal_state) &
            (subset['SWPVALID'] == 0) &
            (subset['IFNUM'] == ifnum) &
            (subset['PLNUM'] == plnum)
        )
        return subset[mask]

    def get_cal_arrays(self, ifnum, plnum):
        """
        Extract pre/post calibration ON/OFF arrays.
        Returns:
            [[pre_on, pre_off], [post_on, post_off]] as tuples or None.
        """
        results = []

        for start, label in [(0, "pre"), (self.post_calibration_start_index, "post")]:
            on_subset = self._mask_cal_data(start, ifnum, plnum, cal_state=1)
            off_subset = self._mask_cal_data(start, ifnum, plnum, cal_state=0)

            try:
                on_array = self._get_time_and_mean(on_subset)
            except Exception:
                on_array = None

            try:
                off_array = self._get_time_and_mean(off_subset)
            except Exception:
                off_array = None

            results.append([on_array, off_array])

        return results

    # ----------------- Core calibration -----------------

    def _compute_delta(self, on_array, off_array):
        """
        Given ON and OFF arrays, compute calibration delta and uncertainty.
        """
        if on_array is None or off_array is None:
            return None, None
        if on_array[0].size == 0 or off_array[0].size == 0:
            return None, None

        cal_on_params, on_uncertainty = self.perform_rcr(on_array[0], on_array[1])
        cal_off_params, off_uncertainty = self.perform_rcr(off_array[0], off_array[1])

        # Mean calibration time
        cal_time = np.mean([np.mean(on_array[0]), np.mean(off_array[0])])
        cal_time = np.clip(
            cal_time,
            (on_array[0][0] + off_array[0][-1]) / 2,
            (on_array[0][-1] + off_array[0][0]) / 2
        )

        # Compute delta
        cal_delta = np.abs(
            (cal_on_params[1] * (cal_time - np.mean(on_array[0])) + cal_on_params[0]) -
            (cal_off_params[1] * (cal_time - np.mean(off_array[0])) + cal_off_params[0])
        )

        # Uncertainty propagation
        cal_uncertainty = np.sqrt(
            on_uncertainty[0] ** 2 + off_uncertainty[0] ** 2 +
            (on_uncertainty[1] * cal_time) ** 2 + (off_uncertainty[1] * cal_time) ** 2
        )

        return cal_delta, cal_uncertainty

    def gain_calibrate_channel(self, channel):
        """
        Perform gain calibration for a single IFNUM/PLNUM channel subset.
        """
        ifnum = channel['IFNUM'][0]
        plnum = channel['PLNUM'][0]

        # Time array relative to first observation in channel
        t0 = Time(channel['DATE-OBS'], format='isot')
        time_array = (Time(channel['DATE-OBS'], format='isot') - t0).sec

        # Get calibration arrays
        pre_arrays, post_arrays = self.get_cal_arrays(ifnum, plnum)

        # Compute deltas
        pre_cal_delta, pre_cal_unc = self._compute_delta(pre_arrays[0], pre_arrays[1])
        post_cal_delta, post_cal_unc = self._compute_delta(post_arrays[0], post_arrays[1])

        # Apply calibration scaling
        if pre_cal_delta and post_cal_delta:
            # Compare deltas
            z_value = abs(pre_cal_delta - post_cal_delta) / np.sqrt(pre_cal_unc ** 2 + post_cal_unc ** 2)
            if z_value < 0.6745:
                # Weighted average
                weights = np.array([1 / pre_cal_unc, 1 / post_cal_unc])
                delta = np.average([pre_cal_delta, post_cal_delta], weights=weights)
                channel['DATA'] /= delta
            else:
                # Interpolate between pre and post
                for idx, t in enumerate(time_array):
                    delta = pre_cal_delta + (post_cal_delta - pre_cal_delta) * (t - t0.sec)
                    channel['DATA'][idx] /= delta
        elif pre_cal_delta:
            channel['DATA'] /= pre_cal_delta
        elif post_cal_delta:
            channel['DATA'] /= post_cal_delta

        # Write calibrated channel back
        mask = (self.data['IFNUM'] == ifnum) & (self.data['PLNUM'] == plnum)
        self.data[mask] = channel

    def gain_calibrate(self):
        """
        Run calibration for all (IFNUM, PLNUM) combinations.
        """
        if self.data is None:
            raise RuntimeError("No data loaded.")

        ifnums = np.unique(self.data['IFNUM'])
        plnums = np.unique(self.data['PLNUM'])

        for ifnum in ifnums:
            for plnum in plnums:
                mask = (self.data['IFNUM'] == ifnum) & (self.data['PLNUM'] == plnum)
                channel = self.data[mask]
                self.gain_calibrate_channel(channel)


if __name__ == "__main__":
    filepath = "C:/Users/starb/Downloads/0136873_validate.fits"
    calibrator = GainCalibrator(filepath)
    start = time.time() 
    calibrator.gain_calibrate()
    end = time.time()

    print(f"Elapsed time: {end - start:.6f} seconds")