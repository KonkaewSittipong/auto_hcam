import logging
import warnings
import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
from pathlib import Path
from astropy import time, coordinates as co, units as u
import os
import glob


def weighted_mean(values, weights, axis=None):
    return np.nansum(values * weights, axis) / np.nansum(weights, axis)


class Ensemble():
    def __init__(self, file=None, save_path=None, solvwcs=False, config=None, diagnostics=False):
        self.file = file
        # FIX: guard against file=None
        self.save_path = (save_path if save_path is not None
                          else (os.path.dirname(file) if file else '.'))
        self.solvwcs = solvwcs
        self.diagnostics = diagnostics

        self.open_df = pd.DataFrame()
        self.lines = []
        self.location = co.EarthLocation.of_site('TNO')
        self.src_pos = co.SkyCoord('21 29 45.0469239408 -04 29 06.973600236',
                                   unit=(u.hourangle, u.deg))
        self.exc_frame = None
        self.sigma_clip = 3

        if config:
            self.read_config_file(config)

        if self.save_path:
            self.save_path, self.figs_dir = self.makedir()

    # ------------------------------------------------------------------
    def makedir(self):
        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path, exist_ok=True)

        base_name = "run"
        search_pattern = os.path.join(self.save_path, f"{base_name}[0-9][0-9][0-9]")
        existing_runs = glob.glob(search_pattern)

        if not existing_runs:
            next_run_num = 1
        else:
            run_numbers = [int(os.path.basename(r).replace(base_name, ""))
                           for r in existing_runs]
            next_run_num = max(run_numbers) + 1

        current_run_dir = os.path.join(self.save_path, f"{base_name}{next_run_num:03d}")
        figs_dir = os.path.join(current_run_dir, "figs")
        os.makedirs(figs_dir, exist_ok=True)
        print(f"[New Experiment] Created directory: {current_run_dir}")
        return current_run_dir, figs_dir

    # ------------------------------------------------------------------
    def read_config_file(self, config):
        # FIX: all parameter parsing moved inside the for loop
        try:
            with open(config, 'r') as conf:
                self.lines = conf.readlines()

            for line in self.lines:
                parts = line.split()
                if len(parts) == 0:
                    continue

                par = parts[0]

                if par.lower() in ('location', 'loc'):
                    try:
                        self.location = co.EarthLocation.of_site(parts[1])
                    except co.errors.UnknownSiteException:
                        try:
                            self.location = co.EarthLocation.from_geodetic(
                                parts[1], parts[2], height=parts[3])
                        except Exception:
                            # FIX: warn() undefined → warnings.warn
                            warnings.warn(
                                f"Location '{parts[1]}' not recognised, "
                                "proceeding without barycentering!")

                elif par.lower() == 'exc_frame':
                    if self.exc_frame is None:
                        self.exc_frame = {}
                    self.exc_frame[parts[1]] = parts[2:]

                elif par.lower() == 'sigma_clip':
                    self.sigma_clip = float(parts[1])

                elif par.lower() == 'pos':
                    self.src_pos = co.SkyCoord(
                        ' '.join(parts[1:]), unit=(u.hourangle, u.deg))

        except FileNotFoundError:
            print(f"Error: Config file {config} not found.")

    # ------------------------------------------------------------------
    def read_log_file(self, logfile, instrument="hcam"):
        headers = None
        try:
            with open(logfile, 'r') as logf:
                for line in logf:
                    line_strip = line.strip()
                    # Find Comment headers
                    if line_strip.startswith("#") and " = CCD" in line_strip:
                        headers = line_strip.split()[3:]
                        break
        except FileNotFoundError:
            print(f"Error: Log file {logfile} not found.")
            return

        if not headers:
            print("Could not find header line (starting with '#' and containing ' = CCD')")
            return
        # Target headers
        target_prefixes = ('counts_', 'countse_', 'sky_', 'flag_', 'x_', 'y_')
        base_cols = ['CCD', 'nframe', 'MJD', 'MJDok', 'Exptim']
        cols_to_keep = [c for c in headers
                        if c in base_cols or c.startswith(target_prefixes)]

        self.open_df = pd.read_csv(
            logfile, sep=r'\s+', comment='#',
            names=headers, header=None,
            usecols=cols_to_keep, float_precision='high',
        )
        self.open_df = self.open_df[self.open_df['MJDok'] != 0].copy()
        self.open_df.reset_index(drop=True, inplace=True)

        self.filter_data()
        self.get_airmass()
        # BJD computed in run() after solve_wcs so the source position is accurate

    # ------------------------------------------------------------------
    def filter_data(self):
        if self.open_df.empty:
            print("No data to filter.")
            return

        flag_cols = [c for c in self.open_df.columns if c.startswith('flag_')]
        for f_col in flag_cols:
            star_idx = f_col.split('_')[1]
            bad_mask = self.open_df[f_col] != 0
            for col in [f'counts_{star_idx}', f'countse_{star_idx}',
                        f'sky_{star_idx}']:
                if col in self.open_df.columns:
                    self.open_df.loc[bad_mask, col] = np.nan

    # ------------------------------------------------------------------
    def get_instrumental_mags(self, diagnostics=False):
        counts_cols = [c for c in self.open_df.columns if c.startswith('counts_')]
        exptim = self.open_df['Exptim']
        new_cols = {}
        for col in counts_cols:
            n = col.split('_')[1]
            counts  = self.open_df[f'counts_{n}']
            countse = self.open_df[f'countse_{n}']
            safe_counts = np.where(counts > 0, counts, np.nan)
            new_cols[f'instru_{n}']  = -2.5 * np.log10(safe_counts / exptim)
            new_cols[f'instrue_{n}'] = (2.5 / np.log(10)) * (countse / safe_counts)
        self.open_df = pd.concat(
            [self.open_df, pd.DataFrame(new_cols, index=self.open_df.index)], axis=1
        )

    # ------------------------------------------------------------------
    def barycenter_times(self):
        if self.open_df.empty or 'MJD' not in self.open_df.columns:
            return
        t = time.Time(self.open_df['MJD'], scale='utc', format='mjd',
                      location=self.location)
        ssbcorr = t.light_travel_time(self.src_pos)
        self.open_df['BJD'] = (t.tdb + ssbcorr).value

    # ------------------------------------------------------------------
    def get_airmass(self):
        if self.open_df.empty or 'MJD' not in self.open_df.columns:
            return
        times = time.Time(self.open_df['MJD'], format='mjd', scale='utc',
                          location=self.location)
        altaz_frame  = co.AltAz(obstime=times, location=self.location)
        target_altaz = self.src_pos.transform_to(altaz_frame)
        self.open_df['secz'] = target_altaz.secz.value
        return self.open_df

    # ------------------------------------------------------------------
    def fit_airmass_coeff1(self, ee, weights, airmass):
        exp_weight = np.nansum(weights, axis=1)
        M = np.zeros((2, 2))
        M[0, 0] = np.nansum(exp_weight)
        M[0, 1] = np.nansum(exp_weight * airmass)
        M[1, 0] = M[0, 1]
        M[1, 1] = np.nansum(exp_weight * (airmass ** 2.0))
        R = np.array([np.nansum(exp_weight * ee),
                      np.nansum(exp_weight * ee * airmass)])
        P = np.linalg.solve(M, R)
        return P

    # ------------------------------------------------------------------
    def solve_ensemble(self):
        # FIX: prefix corrected instru_ / instrue_
        mag_cols = [c for c in self.data.columns if c.startswith('instru_')]
        err_cols = [c for c in self.data.columns if c.startswith('instrue_')]
        m_obs   = self.data[mag_cols].values
        err_obs = self.data[err_cols].values

        W_all = np.where(np.isnan(err_obs) | (err_obs == 0),
                         0.0, 1.0 / (err_obs ** 2))

        mean_star_1 = weighted_mean(m_obs[:, 0], W_all[:, 0])
        if np.isnan(mean_star_1):
            mean_star_1 = 0.0

        m_obs -= mean_star_1
        M_obs_all = np.nan_to_num(m_obs, nan=0.0)

        E, S = M_obs_all.shape
        A = S - 1
        num_unknowns = A + E
        W_solve = W_all[:, 1:]
        M_solve = M_obs_all[:, 1:]

        M = np.zeros((num_unknowns, num_unknowns))
        R = np.zeros(num_unknowns)

        np.fill_diagonal(M[:A, :A], np.sum(W_solve, axis=0))
        np.fill_diagonal(M[A:, A:], np.sum(W_all,   axis=1))
        M[:A, A:] = W_solve.T
        M[A:, :A] = W_solve
        R[:A] = np.sum(M_solve * W_solve, axis=0)
        R[A:] = np.sum(M_obs_all * W_all, axis=1)

        theta, _, _, _ = np.linalg.lstsq(M, R, rcond=None)

        # star 0 is gauge-fixed at 0 in shifted system; stars 1..S-1 = theta[:A]
        mean_mags_shifted = np.zeros(S)
        mean_mags_shifted[1:] = theta[:A]
        meanmags = mean_mags_shifted + mean_star_1
        resids = m_obs - mean_mags_shifted[np.newaxis, :] - theta[A:][:, np.newaxis]
        self.data['exposure_corr'] = theta[A:]

        return resids
    # ------------------------------------------------------------------
 

    # ------------------------------------------------------------------
    # FIX: removed duplicate filter_by_sky; kept only the correct version
    # def filter_by_sky(self, tolerance=0.03):
    #     lims = [1.0 - tolerance, 1.0 + tolerance]
    #     sky_cols = [c for c in self.data.columns if c.startswith('sky_')]
    #     if not sky_cols: return []

    #     median_sky_rate = self.data[sky_cols].median(axis=1) / self.data['Exptim']
    #     diag = {}

    #     for col in sky_cols:
    #         sn = col.split('_')[1]
    #         ratio = (self.data[col] / self.data['Exptim']) / median_sky_rate
    #         med = ratio.median()
    #         diag[sn] = {'ser': ratio, 'med': med, 'bad': ~ratio.between(*lims), 'drop': not (lims[0] <= med <= lims[1])}

    #     if getattr(self, 'diagnostics', False):
    #         nc, ns = 4, len(sky_cols)
    #         nr = int(np.ceil(ns / nc))
    #         fig, axes = plt.subplots(nr, nc, figsize=(20, nr * 3.5), sharex=True)
    #         axes = np.atleast_1d(axes).flatten()

    #         for ax, (sn, info) in zip(axes, diag.items()):
    #             c = 'crimson' if info['drop'] else 'seagreen'
    #             ax.set_facecolor('#fff0f0' if info['drop'] else 'white')
    #             ax.plot(self.data.index, info['ser'], color=c, lw=0.8, alpha=0.6)
                
    #             if info['bad'].any():
    #                 ax.scatter(self.data.index[info['bad']], info['ser'][info['bad']], color='red', s=10, zorder=5)
                
    #             [ax.axhline(l, color='crimson', ls='-', lw=1, alpha=0.3) for l in lims]
    #             ax.set_title(f"Star {sn}\nRatio: {info['med']:.3f}{'' if not info['drop'] else ' REMOVED'}", 
    #                         fontsize=10, fontweight='bold', color='crimson' if info['drop'] else 'black')

    #         [fig.delaxes(ax) for ax in axes[ns:]]
    #         fig.suptitle(rf'Sky Diagnostics ($\pm {tolerance}$)', fontsize=22, y=1.02)
    #         plt.tight_layout()
    #         if hasattr(self, 'figs_dir'):
    #             plt.savefig(os.path.join(self.figs_dir, 'sky_diagnostics.png'), dpi=300, bbox_inches='tight')
    #         plt.show()

    #     return [sn for sn, info in diag.items() if info['drop']]

    def filter_by_sky(self, tolerance=0.03):
        lower_limit, upper_limit = 1.0 - tolerance, 1.0 + tolerance
        sky_cols = [c for c in self.data.columns if c.startswith('sky_')]
        if not sky_cols:
            return []

        median_sky_rate = (self.data[sky_cols].median(axis=1)
                           / self.data['Exptim'])
        star_diagnostics_data = {}

        for col in sky_cols:
            star_num    = col.split('_')[1]
            star_sky_rate = self.data[col] / self.data['Exptim']
            ratio_series  = star_sky_rate / median_sky_rate
            bad_mask      = ((ratio_series < lower_limit) |
                             (ratio_series > upper_limit))
            rel_sky_med   = ratio_series.median()
            star_diagnostics_data[star_num] = {
                'series':   ratio_series,
                'median':   rel_sky_med,
                'bad_mask': bad_mask,
                'drop':     rel_sky_med < lower_limit or rel_sky_med > upper_limit,
            }

        if getattr(self, 'diagnostics', False):
            n_stars = len(sky_cols)
            ncols   = 4
            nrows   = int(np.ceil(n_stars / ncols))
            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(20, nrows * 3.5), sharex=True)
            axes = np.atleast_1d(axes).flatten()

            for i, (star_num, info) in enumerate(star_diagnostics_data.items()):
                ax         = axes[i]
                is_removed = info['drop']
                line_color = 'crimson' if is_removed else 'seagreen'
                bg_color   = '#fff0f0' if is_removed else 'white'
                ax.set_facecolor(bg_color)
                ax.plot(self.data.index, info['series'],
                        color=line_color, lw=0.8, alpha=0.6)
                if info['bad_mask'].any():
                    ax.scatter(self.data.index[info['bad_mask']],
                               info['series'][info['bad_mask']],
                               color='red', s=10, zorder=5)
                ax.axhline(upper_limit, color='crimson', ls='-', lw=1, alpha=0.3)
                ax.axhline(lower_limit, color='crimson', ls='-', lw=1, alpha=0.3)
                title_text  = f"Star {star_num}\nRatio: {info['median']:.3f}"
                title_color = 'crimson' if is_removed else 'black'
                if is_removed:
                    title_text += "\n REMOVED"
                ax.set_title(title_text, fontsize=10, fontweight='bold',
                             color=title_color)

            for j in range(i + 1, len(axes)):
                fig.delaxes(axes[j])
            fig.suptitle(rf'Sky Diagnostics ($\pm {tolerance}$)',
                         fontsize=22, y=1.02)
            plt.tight_layout()
            if getattr(self, 'figs_dir', None):
                plt.savefig(os.path.join(self.figs_dir, 'sky_diagnostics.png'),
                            dpi=300, bbox_inches='tight')
            plt.show()

        return [sn for sn, info in star_diagnostics_data.items() if info['drop']]

    # ------------------------------------------------------------------
    # FIX: indentation error corrected; column name ecounts_ → countse_
    def filter_by_snr(self, min_snr=10.0):
        counts_cols = [c for c in self.data.columns if c.startswith('counts_')]
        if not counts_cols:
            return []

        snr_data = {}
        for col in counts_cols:
            star_num = col.split('_')[1]
            # FIX: was 'ecounts_' — correct column from log is 'countse_'
            err_col  = f'countse_{star_num}'
            if err_col in self.data.columns:
                snr_series = self.data[col] / (self.data[err_col] + 1e-10)
                median_snr = snr_series.median()
                snr_data[star_num] = {
                    'snr':  median_snr,
                    'drop': (median_snr < min_snr) or np.isnan(median_snr),
                }
            else:
                print(f"Warning: error column '{err_col}' not found "
                      f"for star {star_num}")

        if getattr(self, 'diagnostics', False):
            nums     = list(snr_data.keys())
            snr_vals = [info['snr']  for info in snr_data.values()]
            colors   = ['crimson' if info['drop'] else 'seagreen'
                        for info in snr_data.values()]
            plt.figure(figsize=(12, 6))
            plt.bar(nums, snr_vals, color=colors, alpha=0.7)
            plt.axhline(min_snr, color='red', ls='--', lw=2,
                        label=f'Min S/N ({min_snr})')
            plt.title(f'S/N Ratio | Min: {min_snr}',
                      fontsize=14, fontweight='bold')
            plt.ylabel('Signal-to-Noise Ratio')
            plt.xlabel('Star ID')
            plt.yscale('log')
            plt.grid(axis='y', ls=':', alpha=0.6)
            plt.legend()
            if getattr(self, 'figs_dir', None):
                plt.savefig(os.path.join(self.figs_dir, 'snr_diagnostics.png'),
                            dpi=300, bbox_inches='tight')
            plt.show()

        return [sn for sn, info in snr_data.items() if info['drop']]

    # ------------------------------------------------------------------
    # FIX: added save_folder param; fixed prefixes instrumag_→instru_,
    #      einstrumag_→instrue_; fixed ax3 indentation
    def plot_all_comparison_lr(self, all_stars=True, xlim=None, save_folder=None):
        mag_cols = ([c for c in self.df_keep.columns if c.startswith('instru_')]
                    if all_stars else self.surviving_cols)
        if not mag_cols:
            return

        has_cal = 'exposure_corr' in self.data.columns
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(20, 15))
        colors = plt.cm.turbo(np.linspace(0, 1, len(mag_cols)))

        exp_corr = (self.data['exposure_corr'].reindex(self.df_keep.index)
                    if has_cal else None)

        for i, col in enumerate(mag_cols):
            star_id    = col.split('_')[-1]
            star_color = colors[i]
            err_col    = f'instrue_{star_id}'
            valid_idx  = self.df_keep[col].dropna().index

            raw_data = self.df_keep.loc[valid_idx, col]
            raw_err  = (self.df_keep.loc[valid_idx, err_col]
                        if err_col in self.df_keep.columns else None)

            if len(raw_data) > 0:
                ax1.errorbar(raw_data.index, raw_data, yerr=raw_err, fmt='o',
                             markersize=2, color=star_color, ecolor=star_color,
                             alpha=0.5, elinewidth=0.8, capsize=0, zorder=2)
                ax1.axhline(np.nanmedian(raw_data), color=star_color,
                            linestyle='--', alpha=0.3, zorder=1)

            if has_cal:
                cal_data      = self.df_keep[col] - exp_corr
                valid_cal_idx = cal_data.dropna().index
                cal_data      = cal_data.loc[valid_cal_idx]
                cal_err       = (self.df_keep.loc[valid_cal_idx, err_col]
                                 if err_col in self.df_keep.columns else None)
                if len(cal_data) > 0:
                    ax2.errorbar(cal_data.index, cal_data, yerr=cal_err,
                                 fmt='o-', markersize=3,
                                 color=star_color, ecolor=star_color,
                                 alpha=0.6, elinewidth=0.8, capsize=0,
                                 label=f'Star {star_id}', zorder=2)
                    ax2.axhline(np.nanmedian(cal_data), color=star_color,
                                linestyle='--', alpha=0.3, zorder=1)

        ax1.set_title('BEFORE: Raw Instrumental Mag', fontsize=14, fontweight='bold')
        ax1.invert_yaxis()
        ax1.grid(True, ls=':', alpha=0.6)
        if xlim is not None:
            ax1.set_xlim(xlim)

        if has_cal:
            ax2.set_title('AFTER: Ensemble Calibrated Mag',
                          fontsize=14, fontweight='bold')
            ax2.invert_yaxis()
            ax2.grid(True, ls=':', alpha=0.6)
            ax2.legend(loc='upper left', bbox_to_anchor=(1.02, 1),
                       ncol=2, fontsize=8, markerscale=2, title="Star IDs")
            if xlim is not None:
                ax2.set_xlim(xlim)

            ax3.plot(self.data.index, self.data['exposure_corr'],
                     color='teal', marker='.', markersize=4,
                     linestyle='-', alpha=0.7,
                     label='Exposure Correction (Mean Level)')
            ax3.set_title('Atmospheric Correction', fontsize=14, fontweight='bold')
            ax3.set_ylabel('Correction (mag)')
            ax3.set_xlabel('Frame Number')
            ax3.invert_yaxis()
            ax3.grid(True, ls=':', alpha=0.6)
            ax3.legend(loc='upper right')
            if xlim is not None:
                ax3.set_xlim(xlim)

        fig.suptitle('Global Photometry Comparison', fontsize=18, y=0.98)
        plt.tight_layout(rect=[0, 0, 0.9, 1])

        out_dir = save_folder or getattr(self, 'figs_dir', self.save_path)
        if out_dir:
            save_name = os.path.join(out_dir, 'all_stars_comparison.png')
            plt.savefig(save_name, dpi=300, bbox_inches='tight')
            print(f"Saved: {save_name}")
        plt.show()

    # ------------------------------------------------------------------
    # FIX: prefixes instrumag_→instru_, einstrumag_→instrue_;
    #      reindex exposure_corr to df_keep index
    def save_results(self):
        if 'exposure_corr' not in self.data.columns:
            print("Warning: 'exposure_corr' not found. Run solve_ensemble() first.")
            return

        exp_corr  = self.data['exposure_corr'].reindex(self.df_keep.index)
        time_cols = [c for c in ['MJD', 'BJD', 'Exptim', 'secz']
                     if c in self.df_keep.columns]
        df_out = self.df_keep[time_cols].copy()

        new_cols = {'exposure_corr': exp_corr}
        mag_cols = [c for c in self.df_keep.columns if c.startswith('instru_')]
        for col in mag_cols:
            star_id = col.split('_')[-1]
            err_col = f'instrue_{star_id}'
            new_cols[f'mag_{star_id}']  = self.df_keep[col] - exp_corr
            if err_col in self.df_keep.columns:
                new_cols[f'emag_{star_id}'] = self.df_keep[err_col]
            # Append median RA/DEC for this aperture if WCS was solved
            ra_col  = f'ra_{star_id}'
            dec_col = f'dec_{star_id}'
            if ra_col in self.df_keep.columns:
                new_cols[f'RA_{star_id}']  = self.df_keep[ra_col]
                new_cols[f'DEC_{star_id}'] = self.df_keep[dec_col]
        df_out = pd.concat(
            [df_out, pd.DataFrame(new_cols, index=self.df_keep.index)], axis=1
        )

        csv_path = os.path.join(self.save_path, 'calibrated_lightcurves.txt')
        with open(csv_path, 'w') as f:
            if hasattr(self, 'surviving_stars') and self.surviving_stars:
                f.write("# Ensemble reference stars: "
                        + ", ".join(self.surviving_stars) + "\n")
            else:
                f.write("# Warning: No reference stars found.\n")
            df_out.to_csv(f, index=False)

        print(f"Saved: {csv_path}")

    # ------------------------------------------------------------------
    def plot_rms_history(self, target_rms=None, save_folder=None):
        if not hasattr(self, 'history') or not self.history:
            print("No RMS history found. Run the pipeline first.")
            return

        fig, ax    = plt.subplots(figsize=(10, 6))
        iterations = [h['iteration'] for h in self.history]
        rms_vals   = [h['rms']       for h in self.history]
        reject_ids = [h['reject']    for h in self.history]

        ax.plot(iterations, rms_vals, marker='o', linestyle='-',
                color='indigo', linewidth=2, markersize=8, alpha=0.8,
                label='Ensemble RMS')

        for it, rms, sid in zip(iterations, rms_vals, reject_ids):
            ax.annotate(f'drop {sid}',
                        xy=(it, rms),
                        xytext=(8, 6), textcoords='offset points',
                        fontsize=8, color='crimson')

        if target_rms is not None:
            ax.axhline(target_rms, color='crimson', linestyle='--', linewidth=2,
                       alpha=0.8, label=f'Target RMS ({target_rms})')

        ax.legend(fontsize=12)
        ax.set_title('Pipeline Convergence: Ensemble RMS per Iteration',
                     fontsize=16, fontweight='bold')
        ax.set_xlabel('Cleaning Iteration')
        ax.set_ylabel('Global RMS (mag)')
        ax.grid(True, linestyle=':', alpha=0.7)
        ax.set_xticks(iterations)
        plt.tight_layout()

        if save_folder:
            save_name = os.path.join(save_folder, 'rms_convergence_history.png')
            plt.savefig(save_name, dpi=300, bbox_inches='tight')
            print(f"Saved: {save_name}")
        plt.show()

    # ------------------------------------------------------------------
    def solve_wcs(self, ra_center=None, dec_center=None, radius=5.0,
                  follow_radius=0.005,
                  scale_low=None, scale_high=None,
                  astrometry_cache='astrometry_cache',
                  binning=1, verbose=False, single_frame=False):
        """
        Solve WCS per-frame from x_N / y_N columns in self.open_df using the
        astrometry Python package (no image file needed).
        Frame 0 is solved blind (or with position hint); subsequent frames
        use the previous frame's RA/Dec centre as a hint.
        On success, ra_N / dec_N columns are written into self.open_df.
        """
        try:
            import astrometry
        except ImportError:
            print("[WCS] 'astrometry' package not found. "
                  "Install with: pip install astrometry")
            return

        x_cols = sorted(c for c in self.open_df.columns if c.startswith('x_'))
        if not x_cols:
            print("[WCS] No x_ columns in open_df — run read_log_file first.")
            return

        ap_labels = [c.split('_')[1] for c in x_cols
                     if f"y_{c.split('_')[1]}" in self.open_df.columns]
        n_frames  = len(self.open_df)
        n_aps     = len(ap_labels)
        print(f"[WCS] {n_aps} apertures × {n_frames} frames")

        if n_aps < 6:
            print(f"[WCS] Need ≥6 apertures for plate solving, got {n_aps}.")
            return

        sl = (scale_low  if scale_low  is not None else 0.40) * binning
        sh = (scale_high if scale_high is not None else 0.50) * binning

        def _idx(series, scales, label):
            try:
                files   = series.index_files(
                    cache_directory=astrometry_cache, scales=scales)
                on_disk = [Path(f) for f in files if Path(f).exists()]
                print(f"[WCS]   {label}: {len(on_disk)} on disk")
                return on_disk
            except Exception as e:
                print(f"[WCS]   {label}: {e}")
                return []

        print("[WCS] Loading index files...")
        index_files = list(dict.fromkeys(
            _idx(astrometry.series_5200, {2, 3, 4}, "series_5200")
            + _idx(astrometry.series_4200, {2, 3, 4}, "series_4200")
        ))
        if not index_files:
            print(f"[WCS] No index files found in '{astrometry_cache}'.")
            return

        size_hint = astrometry.SizeHint(lower_arcsec_per_pixel=sl,
                                        upper_arcsec_per_pixel=sh)
        params = astrometry.SolutionParameters(
            logodds_callback=lambda ll: (
                astrometry.Action.STOP if ll[0] > 100.0
                else astrometry.Action.CONTINUE
            ),
        )

        if verbose:
            logging.basicConfig(level=logging.INFO,
                                format="%(message)s", force=True)
            logging.getLogger().setLevel(logging.INFO)

        ra_out  = {ap: np.full(n_frames, np.nan) for ap in ap_labels}
        dec_out = {ap: np.full(n_frames, np.nan) for ap in ap_labels}
        prev_ra, prev_dec = ra_center, dec_center
        prev_wcs = None
        field_solved = False   # shrink radius after first successful solve

        if single_frame:
            print(f"[WCS] single_frame=True — solving frame 0 only, applying to all {n_frames} frames...")
        else:
            print(f"[WCS] Solving {n_frames} frames...")

        with astrometry.Solver(index_files) as solver:
            for fi, row in enumerate(self.open_df.itertuples(index=False)):

                # ── single_frame mode: skip solving after frame 0 ──────────
                if single_frame and fi > 0:
                    if prev_wcs is None:
                        print(f"[WCS] Frame 0 solve failed — cannot apply to remaining frames.")
                        return
                    wcs = prev_wcs
                    for ap in ap_labels:
                        x_raw = getattr(row, f'x_{ap}', np.nan)
                        y_raw = getattr(row, f'y_{ap}', np.nan)
                        if np.isfinite(x_raw) and np.isfinite(y_raw):
                            try:
                                sky = wcs.pixel_to_world(x_raw / binning,
                                                         y_raw / binning)
                                ra_out[ap][fi]  = float(sky.ra.deg)
                                dec_out[ap][fi] = float(sky.dec.deg)
                            except Exception:
                                pass
                    if fi % 10 == 0 or fi == n_frames - 1:
                        print(f"[WCS]   frame {fi + 1:04d}/{n_frames} | applying frame-0 WCS")
                    continue

                # ── per-frame solve ────────────────────────────────────────
                stars = []
                for ap in ap_labels:
                    x = getattr(row, f'x_{ap}', np.nan) / binning
                    y = getattr(row, f'y_{ap}', np.nan) / binning
                    if np.isfinite(x) and np.isfinite(y):
                        stars.append((x, y))

                position_hint = None
                if prev_ra is not None and prev_dec is not None:
                    search_r = follow_radius if field_solved else radius
                    position_hint = astrometry.PositionHint(
                        ra_deg=prev_ra, dec_deg=prev_dec, radius_deg=search_r)

                wcs = None
                newly_solved = False
                if len(stars) >= 6:
                    solution = solver.solve(
                        stars=stars,
                        size_hint=size_hint,
                        position_hint=position_hint,
                        solution_parameters=params,
                    )
                    if solution.has_match():
                        match        = solution.best_match()
                        wcs          = match.astropy_wcs()
                        prev_ra      = match.center_ra_deg
                        prev_dec     = match.center_dec_deg
                        prev_wcs     = wcs
                        newly_solved = True
                        if not field_solved:
                            field_solved = True
                            print(f"[WCS] Frame {fi} solved: "
                                  f"RA={prev_ra:.5f}  Dec={prev_dec:.5f} "
                                  f"→ follow radius={follow_radius}°")
                    elif fi == 0:
                        hint_str = (f"RA={prev_ra:.4f} Dec={prev_dec:.4f}"
                                    if prev_ra is not None else "no hint")
                        print(f"[WCS] Frame 0 solve failed ({hint_str}). "
                              f"Scale hint [{sl:.3f}, {sh:.3f}] arcsec/px.")
                        if verbose:
                            logging.getLogger().setLevel(logging.WARNING)
                        return

                if wcs is None:
                    wcs = prev_wcs

                for ap in ap_labels:
                    x_raw = getattr(row, f'x_{ap}', np.nan)
                    y_raw = getattr(row, f'y_{ap}', np.nan)
                    if wcs is not None and np.isfinite(x_raw) and np.isfinite(y_raw):
                        try:
                            sky = wcs.pixel_to_world(x_raw / binning,
                                                     y_raw / binning)
                            ra_out[ap][fi]  = float(sky.ra.deg)
                            dec_out[ap][fi] = float(sky.dec.deg)
                        except Exception:
                            pass

                status = "OK" if newly_solved else "fallback"
                print(f"[WCS]   frame {fi + 1:04d}/{n_frames} | "
                      f"RA={prev_ra:.5f} Dec={prev_dec:.5f} | {status}")

        if verbose:
            logging.getLogger().setLevel(logging.WARNING)

        self.star_radec = {}
        radec_cols = {}
        for ap in ap_labels:
            radec_cols[f'ra_{ap}']  = ra_out[ap]
            radec_cols[f'dec_{ap}'] = dec_out[ap]
            med_ra  = np.nanmedian(ra_out[ap])
            med_dec = np.nanmedian(dec_out[ap])
            if np.isfinite(med_ra):
                self.star_radec[ap] = co.SkyCoord(
                    ra=med_ra * u.deg, dec=med_dec * u.deg, frame='icrs')
                print(f"  Star {ap}: RA={med_ra:.5f}  Dec={med_dec:.5f}")
        self.open_df = pd.concat(
            [self.open_df, pd.DataFrame(radec_cols, index=self.open_df.index)], axis=1
        )

        # Update src_pos to the first solved aperture so barycenter_times()
        # uses the accurate field position instead of the hardcoded default
        first_ap = ap_labels[0] if ap_labels else None
        if first_ap and first_ap in self.star_radec:
            self.src_pos = self.star_radec[first_ap]
            print(f"[WCS] src_pos updated to Star {first_ap}: "
                  f"RA={self.src_pos.ra.deg:.5f}  Dec={self.src_pos.dec.deg:.5f}")

        print("[WCS] Done. ra_N / dec_N written to open_df.")
    
    def find_most_variable_star(self):
        mag_cols = [c for c in self.data.columns if c.startswith('instru_')]
        if not mag_cols or 'exposure_corr' not in self.data.columns:
            return {}
        corr = self.data['exposure_corr'].values
        variability = {col: np.nanstd(self.data[col].values - corr) for col in mag_cols}
        return dict(sorted(variability.items(), key=lambda x: x[1], reverse=True))


    # ------------------------------------------------------------------
    def run(self, target_rms=0.02, numstars=12, ignor_stars=None,
            ra_center=None, dec_center=None, radius=5.0, follow_radius=0.3,
            scale_low=None, scale_high=None,
            astrometry_cache='astrometry_cache', binning=1, single_frame=False):

        print("=" * 60)
        print(f"[1/7] Reading log file: {self.file}")
        self.read_log_file(logfile=self.file, instrument="hcam")
        n_frames_raw = len(self.open_df)
        n_stars_raw  = len([c for c in self.open_df.columns if c.startswith('counts_')])
        print(f"      → {n_frames_raw} frames, {n_stars_raw} apertures loaded")

        if self.solvwcs:
            mode_str = "single-frame → apply-all" if single_frame else f"per-frame, follow={follow_radius}°"
            print(f"[2/7] Solving WCS (radius={radius}°, mode={mode_str})...")
            self.solve_wcs(ra_center=ra_center, dec_center=dec_center,
                           radius=radius, follow_radius=follow_radius,
                           scale_low=scale_low, scale_high=scale_high,
                           astrometry_cache=astrometry_cache,
                           binning=binning, single_frame=single_frame)
        else:
            print(f"[2/7] WCS solve skipped (solvwcs=False)")

        print(f"[3/7] Computing BJD & instrumental mags...")
        self.barycenter_times()
        print(f"      → BJD using RA={self.src_pos.ra.deg:.5f} Dec={self.src_pos.dec.deg:.5f}")
        self.get_instrumental_mags()
        self.data = self.open_df.copy()
        print(f"      → Instrumental mags computed for {n_stars_raw} stars")

        print(f"[4/7] Initial star filtering...")
        n_before = len([c for c in self.data.columns if c.startswith('instru_')])

        bad_sky_stars = self.filter_by_sky(tolerance=0.03)
        for s_num in set(bad_sky_stars):
            self.data.drop(columns=[f'instru_{s_num}', f'instrue_{s_num}'],
                           errors='ignore', inplace=True)
            print(f"      [Sky]  Drop Star {s_num}")

        low_snr_stars = self.filter_by_snr(min_snr=10.0)
        for s_num in set(low_snr_stars):
            self.data.drop(columns=[f'instru_{s_num}', f'instrue_{s_num}'],
                           errors='ignore', inplace=True)
            print(f"      [SNR]  Drop Star {s_num}")

        n_after = len([c for c in self.data.columns if c.startswith('instru_')])
        print(f"      → {n_before - n_after} stars removed, {n_after} remaining")

        print(f"[5/7] Sigma-clipping outlier exposures (σ={self.sigma_clip})...")
        n_exp_before = len(self.data)
        err_cols_init = [c for c in self.data.columns if c.startswith('instrue_')]
        if err_cols_init:
            row_err = self.data[err_cols_init].median(axis=1)
            err_med = row_err.median()
            err_std = row_err.std()
            bad_rows = self.data.index[row_err > err_med + self.sigma_clip * err_std]
            if len(bad_rows) > 0:
                self.data.drop(index=bad_rows, inplace=True)
                print(f"      → Removed {len(bad_rows)} outlier exposure(s) "
                      f"(threshold={err_med + self.sigma_clip * err_std:.4f} mag-err)")
            else:
                print(f"      → No outlier exposures found")
        print(f"      → {len(self.data)}/{n_exp_before} exposures kept")

        self.df_keep = self.data.copy()

        if ignor_stars is not None:
            print(f"      [User] Ignoring stars: {ignor_stars}")
            for t_id in ignor_stars:
                t_id = str(t_id)
                self.data.drop(columns=[f'instru_{t_id}', f'instrue_{t_id}'],
                               errors='ignore', inplace=True)

        print(f"[6/7] Ensemble cleaning loop (target RMS={target_rms}, min stars={numstars})...")
        self.history = []
        cleaning  = True
        iteration = 1

        while cleaning:
            mag_cols = [c for c in self.data.columns if c.startswith('instru_')]
            if len(mag_cols) < 2:
                print("      → Less than 2 stars remaining — stopping.")
                break

            residuals = self.solve_ensemble()
            current_rms = float(np.sqrt(np.nanmean(residuals ** 2)))
            print(f"      [iter {iteration:03d}] {len(mag_cols)} stars | RMS={current_rms:.5f}")

            if current_rms <= target_rms:
                print(f"      → RMS ({current_rms:.5f}) ≤ target ({target_rms}) — converged!")
                break

            if len(mag_cols) <= numstars:
                print(f"      → {len(mag_cols)} stars ≤ min ({numstars}) — stopping.")
                break

            worst_col = list(self.find_most_variable_star().keys())[0]
            star_id   = worst_col.split('_')[-1]
            err_col   = worst_col.replace('instru_', 'instrue_')
            self.data.drop(columns=[worst_col, err_col], errors='ignore', inplace=True)
            print(f"             → Reject Star {star_id} (most variable)")

            iteration += 1
            self.history.append({
                "iteration": iteration,
                "rms": current_rms,
                "reject": star_id,
            })
            if iteration > 100:
                print("      → Max iterations (100) reached.")
                break

        self.surviving_stars = [c.split('_')[-1] for c in self.data.columns
                                if c.startswith('instru_')]
        print(f"      → Final: {len(self.surviving_stars)} reference stars | "
              f"RMS={float(np.sqrt(np.nanmean(residuals**2))):.5f}")
        print(f"      → Surviving stars: {', '.join(self.surviving_stars)}")

        print(f"[7/7] Saving results & plots...")
        self.diag_plot()
        self.plot_rms_history(target_rms=target_rms)
        self.save_results()
        print("=" * 60)
        print("Done.")



    def diag_plot(self):
        mag_cols = [c for c in self.df_keep.columns if c.startswith('instru_')]
        n_stars = len(mag_cols)
        ncols = 4
        nrows = int(np.ceil(n_stars / ncols))

        # map star_id → (iteration, rms_at_rejection)
        reject_info = {str(h['reject']): (h['iteration'], h['rms']) for h in self.history}

        print(f"Analyzing & Plotting {n_stars} stars for variability...")
        fig, axes = plt.subplots(nrows, ncols, figsize=(20, nrows * 3.5), sharex=True)
        if n_stars == 1:
            axes = [axes]
        axes = axes.flatten()

        x = np.arange(len(self.df_keep))

        for ax, col in zip(axes, mag_cols):
            star_id = col.split('_')[-1]
            err_col = col.replace('instru_', 'instrue_')
            info = reject_info.get(star_id)
            is_rejected = info is not None
            color = 'crimson' if is_rejected else 'steelblue'

            y   = self.df_keep[col].values - self.data['exposure_corr'].reindex(self.df_keep.index).values
            err = self.df_keep[err_col].values if err_col in self.df_keep.columns else np.zeros_like(y)
            y_centered = y - np.nanmedian(y)

            ax.errorbar(x, y_centered, yerr=err,
                        fmt='o', ms=2, lw=0, elinewidth=0.8,
                        color=color, ecolor=color, alpha=0.6)
            ax.axhline(0, color='gray', lw=0.5, ls='--')
            ax.set_facecolor('#fff0f0' if is_rejected else 'white')

            if is_rejected:
                iter_num, rms_at = info
                label = f"Star {star_id} [REJECTED iter={iter_num}, rms={rms_at:.4f}]"
            else:
                label = f"Star {star_id}"
            ax.set_title(label, fontsize=8,
                         fontweight='bold',
                         color='crimson' if is_rejected else 'black')
            ax.set_ylabel('Δmag', fontsize=8)

        for ax in axes[n_stars:]:
            fig.delaxes(ax)

        fig.suptitle('Comparison Star Diagnostics', fontsize=14, y=1.01)
        plt.tight_layout()

        if hasattr(self, 'figs_dir'):
            out = os.path.join(self.figs_dir, 'diag_plot.png')
            plt.savefig(out, dpi=150, bbox_inches='tight')
            print(f"Saved: {out}")
        plt.show()

