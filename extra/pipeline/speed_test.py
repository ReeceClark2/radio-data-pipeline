import os
import copy
from collections import defaultdict

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from scipy.stats import linregress
import rcr
from tqdm import tqdm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class Gain_Calibrate:
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.indices = None
        self._pre_cal_deltas = None
        self._post_cal_deltas = None

        try:
            with fits.open(filepath) as hdul:
                self.header = hdul[0].header
                data = Table(hdul[1].data)
                # Filter to only IFNUM == 0 and PLNUM == 0
                data = data[(data['IFNUM'] == 0) & (data['PLNUM'] == 0)]
                self.data = data
                self.old_data = copy.deepcopy(data)
                num_points = len(data[(data['CALSTATE'] == 0) & (data['SWPVALID'] == 0)])
                print(f"[INFO] Points after filtering: {num_points}")
        except Exception as e:
            print(f"[ERROR] Failed to load FITS: {e}")
            self.data = Table()
            self.old_data = Table()

    # --------------------------
    # Calibration index finding
    # --------------------------
    def find_calibrations(self):
        data = self.data
        data_start_ind = None
        post_cal_start_ind = None
        cal_started = False
        pre_cal_complete = False
        counter = 0
        channel_count = len(np.unique(data['IFNUM'])) * len(np.unique(data['PLNUM']))

        for ind, row in enumerate(data):
            if row['CALSTATE'] == 1:
                cal_started = True

            if cal_started and row["CALSTATE"] == 0 and row["SWPVALID"] == 1 and not pre_cal_complete:
                data_start_ind = ind
                pre_cal_complete = True

            if pre_cal_complete and row["SWPVALID"] == 0 and data[ind - 1]["SWPVALID"] == 0:
                if post_cal_start_ind is None:
                    post_cal_start_ind = ind
            else:
                post_cal_start_ind = None

            if pre_cal_complete and row['CALSTATE'] == 0 and row['SWPVALID'] == 1:
                counter += 1
            if counter <= 3 * channel_count and row['SWPVALID'] == 1:
                pre_cal_complete = False
            if pre_cal_complete and row['SWPVALID'] == 0 and row['CALSTATE'] == 1:
                break

        self.indices = [data_start_ind, post_cal_start_ind]

        if self.header.get('OBSMODE') == 'onoff':
            for ind, row in enumerate(data):
                if 'onoff:off' in row['OBSMODE']:
                    self.indices = [data_start_ind, ind - 1, ind, post_cal_start_ind]
                    break

    # --------------------------
    # Linear fit helpers
    # --------------------------
    def linear(self, x, params):
        return params[0] + x * params[1]

    def d_linear_1(self, x, params):
        return 1

    def d_linear_2(self, x, params):
        return x

    def rcr(self, arr):
        x = np.array(arr[0], dtype=float)
        y = np.array(arr[1], dtype=float)
        x -= np.average(x)

        result = linregress(x, y)
        guess = [result.slope, result.intercept]
        model = rcr.FunctionalForm(
            self.linear, x, y,
            [self.d_linear_1, self.d_linear_2],
            guess
        )
        r = rcr.RCR(rcr.SS_MEDIAN_DL)
        r.setParametricModel(model)
        r.performBulkRejection(y)

        indices = r.result.indices
        x = x[indices]
        y = y[indices]
        best_fit_parameters = model.result.parameters

        sigma = (1 / (len(x) - 2)) * np.sum((y - best_fit_parameters[1]*x - best_fit_parameters[0])**2)
        m_sd = np.sqrt(sigma / np.sum((x - np.mean(x))**2))
        b_sd = np.sqrt(sigma * ((1/len(x)) + ((np.mean(x)**2)/np.sum((x-np.mean(x))**2))))
        uncertainties = (b_sd, m_sd)

        return best_fit_parameters, uncertainties

    # --------------------------
    # Main calibration computation
    # --------------------------
    def compute_gain_deltas(self, window: int):
        # Reset data
        self.data = self.old_data.copy()

        if len(self.data) == 0:
            return

        first_len = len(self.data[0]['DATA'])
        if not all(len(row['DATA']) == first_len for row in self.data):
            return

        # Precompute times and arrays
        t0 = Time(self.header['DATE'], format='isot')
        times = np.array([(Time(row['DATE-OBS'], format='isot') - t0).sec for row in self.data])
        IFNUM = np.array(self.data['IFNUM'])
        PLNUM = np.array(self.data['PLNUM'])
        SWPVALID = np.array(self.data['SWPVALID'])
        CALSTATE = np.array(self.data['CALSTATE'])
        DATA = np.stack(self.data['DATA'])

        # Pre-slice into pre and post indices
        idxs = np.arange(len(self.data))
        pre_idx = idxs < self.indices[0]
        post_idx = idxs > self.indices[1]

        # Group indices by IFNUM/PLNUM
        grouped = defaultdict(list)
        for idx, (i, j) in enumerate(zip(IFNUM, PLNUM)):
            grouped[(i, j)].append(idx)
        grouped = {k: np.array(v, dtype=int) for k, v in grouped.items()}

        all_pre_deltas = np.full(first_len, np.nan)
        all_post_deltas = np.full(first_len, np.nan)

        # iterate IFNUM/PLNUM
        for (i, j), gidx in grouped.items():
            gidx_pre = gidx[pre_idx[gidx]]
            gidx_post = gidx[post_idx[gidx]]

            def get_times_vals(indices, swp, cal):
                mask = (SWPVALID[indices] == swp) & (CALSTATE[indices] == cal)
                if not np.any(mask):
                    return np.array([]), np.array([])
                sub_idx = indices[mask]
                return times[sub_idx], DATA[sub_idx, :]

            # pre calibration
            pre_on_t, pre_on_val = get_times_vals(gidx_pre, 0, 1)
            pre_off_t, pre_off_val = get_times_vals(gidx_pre, 0, 0)
            post_on_t, post_on_val = get_times_vals(gidx_post, 0, 1)
            post_off_t, post_off_val = get_times_vals(gidx_post, 0, 0)

            pre_deltas, pre_uncs = [], []
            post_deltas, post_uncs = [], []

            for h in range(first_len):
                try:
                    if len(pre_on_t) and len(pre_off_t):
                        pre_on_params, pre_on_unc = self.rcr((pre_on_t, pre_on_val[:, h]))
                        pre_off_params, pre_off_unc = self.rcr((pre_off_t, pre_off_val[:, h]))
                        pre_cal_time = np.mean([np.mean(pre_on_t), np.mean(pre_off_t)])
                        pre_cal_time = np.clip(pre_cal_time,
                                               (pre_on_t[0] + pre_off_t[-1]) / 2,
                                               (pre_on_t[-1] + pre_off_t[0]) / 2)
                        pre_delta = np.abs(
                            (pre_on_params[1]*(pre_cal_time-np.mean(pre_on_t))+pre_on_params[0]) -
                            (pre_off_params[1]*(pre_cal_time-np.mean(pre_off_t))+pre_off_params[0])
                        )
                        pre_unc = np.sqrt(pre_on_unc[0]**2 + pre_off_unc[0]**2 +
                                          (pre_on_unc[1]*pre_cal_time)**2 + (pre_off_unc[1]*pre_cal_time)**2)
                    else:
                        pre_delta, pre_unc = None, None
                except Exception:
                    pre_delta, pre_unc = None, None

                try:
                    if len(post_on_t) and len(post_off_t):
                        post_on_params, post_on_unc = self.rcr((post_on_t, post_on_val[:, h]))
                        post_off_params, post_off_unc = self.rcr((post_off_t, post_off_val[:, h]))
                        post_cal_time = np.mean([np.mean(post_on_t), np.mean(post_off_t)])
                        post_cal_time = np.clip(post_cal_time,
                                                (post_on_t[0] + post_off_t[-1]) / 2,
                                                (post_on_t[-1] + post_off_t[0]) / 2)
                        post_delta = np.abs(
                            (post_on_params[1]*(post_cal_time-np.mean(post_on_t))+post_on_params[0]) -
                            (post_off_params[1]*(post_cal_time-np.mean(post_off_t))+post_off_params[0])
                        )
                        post_unc = np.sqrt(post_on_unc[0]**2 + post_off_unc[0]**2 +
                                           (post_on_unc[1]*post_cal_time)**2 + (post_off_unc[1]*post_cal_time)**2)
                    else:
                        post_delta, post_unc = None, None
                except Exception:
                    post_delta, post_unc = None, None

                pre_deltas.append(pre_delta)
                pre_uncs.append(pre_unc)
                post_deltas.append(post_delta)
                post_uncs.append(post_unc)

            pre_deltas = np.array(pre_deltas, dtype=object)
            post_deltas = np.array(post_deltas, dtype=object)

            # Smooth neighbors
            def avg_neighbors(values, idx, w):
                if values[idx] is None:
                    return None
                lo = max(0, idx-w)
                hi = min(len(values), idx+w+1)
                arr = [v for v in values[lo:hi] if v is not None]
                return np.mean(arr) if arr else None

            for h in range(first_len):
                pre_avg = avg_neighbors(pre_deltas, h, window)
                pre_unc = avg_neighbors(pre_uncs, h, window)
                post_avg = avg_neighbors(post_deltas, h, window)
                post_unc = avg_neighbors(post_uncs, h, window)

                if pre_avg and post_avg:
                    z = abs(pre_avg-post_avg)/np.sqrt(pre_unc**2+post_unc**2)
                    if z < 0.6745:
                        weights = np.array([1/pre_unc, 1/post_unc])
                        delta = np.average([pre_avg, post_avg], weights=weights)
                        DATA[:, h] /= delta
                    else:
                        for ind in range(len(DATA)):
                            t = times[ind]
                            delta = pre_avg + (post_avg-pre_avg)*(t-pre_on_t.mean())/(post_on_t.mean()-pre_on_t.mean())
                            DATA[ind, h] /= delta
                elif pre_avg and not post_avg:
                    DATA[:, h] /= pre_avg
                elif post_avg and not pre_avg:
                    DATA[:, h] /= post_avg

            all_pre_deltas = np.array([np.nan if v is None else v for v in pre_deltas])
            all_post_deltas = np.array([np.nan if v is None else v for v in post_deltas])

        # Update table
        self.data['DATA'] = DATA
        self._pre_cal_deltas = all_pre_deltas
        self._post_cal_deltas = all_post_deltas

    # --------------------------
    # Animation
    # --------------------------
    def animate_gain_deltas(self, window_range=(1, 100, 1), output_file="gain_deltas_animation_low.mp4"):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        window_sizes = list(range(*window_range))

        # Pre-create animation writer
        writer = animation.FFMpegWriter(fps=10, bitrate=1800)

        with writer.saving(fig, output_file, dpi=200):
            for frame_idx in tqdm(range(len(window_sizes)), desc="Rendering frames"):
                window = window_sizes[frame_idx]
                ax1.clear()
                ax2.clear()
                try:
                    self.compute_gain_deltas(window)
                    mask = (self.data['IFNUM'] == 0) & (self.data['PLNUM'] == 0)
                    filtered_data = self.data[mask]
                    spectrum = np.mean(filtered_data['DATA'], axis=0)
                    ch_indices = np.arange(len(spectrum))
                    pre_deltas = self._pre_cal_deltas
                    post_deltas = self._post_cal_deltas
                    delta_diff = np.abs(pre_deltas - post_deltas)

                    ax1.plot(ch_indices, spectrum, color='black', linewidth=0.6)
                    ax1.set_ylim(0, 300)
                    ax1.set_xlabel("Channel Index")
                    ax1.set_ylabel("Mean Intensity")
                    ax1.set_title(f"Mean Spectrum\nWindow={window}")
                    ax1.grid(True)

                    ax2.plot(ch_indices, pre_deltas, label="Pre-Cal", color='blue', alpha=0.6, linewidth=0.6)
                    ax2.plot(ch_indices, post_deltas, label="Post-Cal", color='orange', alpha=0.6, linewidth=0.6)
                    ax2.plot(ch_indices, delta_diff, label="|Pre-Post|", color='red', alpha=0.8, linewidth=0.6)
                    ax2.set_xlabel("Channel Index")
                    ax2.set_ylabel("Delta")
                    ax2.set_title("Calibration Deltas")
                    ax2.legend()
                    ax2.grid(True)

                except Exception as e:
                    ax1.text(0.5, 0.5, f"Error at window={window}\n{type(e).__name__}",
                            ha='center', va='center', transform=ax1.transAxes, fontsize=12, color='red')
                    ax1.set_axis_off()
                    ax2.set_axis_off()

                fig.tight_layout()
                writer.grab_frame()

        plt.close(fig)

    # --------------------------
    # Saving
    # --------------------------
    def save(self, output_path=None):
        if output_path is None:
            base, ext = os.path.splitext(self.filepath)
            output_path = f"{base}_gain_calibrated{ext}"
        primary_hdu = fits.PrimaryHDU(header=self.header)
        table_hdu = fits.BinTableHDU(data=self.data)
        fits.HDUList([primary_hdu, table_hdu]).writeto(output_path, overwrite=True)

    def gain_calibrate(self):
        self.find_calibrations()
        self.animate_gain_deltas()


if __name__ == "__main__":
    filepath = "C:/Users/starb/Downloads/0136872.fits"
    gc = Gain_Calibrate(filepath)
    gc.gain_calibrate()
