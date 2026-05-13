"""
Multi-epoch ensemble calibration.

Reads calibrated_lightcurves.txt files produced by the per-night pipeline,
cross-matches stars across all epochs using RA/DEC, then re-runs ensemble
calibration on the combined dataset to produce a consistent long-baseline
light curve.

Usage example
-------------
from multi_epoch import MultiEpochEnsemble

epoch_dirs = [
    '/MSSP/sittipong/reduce2/2023-11-06/hcam_reduction_2023_11_06_i/results/run001',
    '/MSSP/sittipong/reduce2/2023-11-07/hcam_reduction_2023_11_07_i/results/run001',
]

me = MultiEpochEnsemble(
    epoch_dirs     = epoch_dirs,
    output_dir     = '/MSSP/sittipong/multi_epoch_results',
    tolerance_arcsec = 2.0,
    target_star_ra   = 322.4377,
    target_star_dec  = -4.4853,
    diagnostics      = True,
)
me.run_multi_epoch(target_rms=0.02, numstars=12)
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import coordinates as co, units as u
from scipy import stats

from ensemble import Ensemble


class MultiEpochEnsemble(Ensemble):
    """Cross-epoch ensemble calibration."""

    LC_FILE      = 'calibrated_lightcurves.txt'

    def __init__(self, epoch_dirs, output_dir,
                 tolerance_arcsec=2.0,
                 target_star_ra=None, target_star_dec=None,
                 diagnostics=False):
        """
        Parameters
        ----------
        epoch_dirs : list[str]
            Each path must point to a run directory that contains
            calibrated_lightcurves.txt (with RA_N / DEC_N columns).
        output_dir : str
            Where to write multi-epoch results.
        tolerance_arcsec : float
            Sky-matching radius used when cross-identifying stars across epochs.
        target_star_ra, target_star_dec : float
            ICRS coordinates of the science target (deg).
            Used to identify which global star ID is the target.
        """
        super().__init__(file=None, save_path=output_dir,
                         solvwcs=False, diagnostics=diagnostics)

        self.epoch_dirs       = epoch_dirs
        self.tolerance_arcsec = tolerance_arcsec
        self.target_ra        = target_star_ra
        self.target_dec       = target_star_dec

        self.epochs           = []   # list of per-epoch dicts
        self.global_catalog   = {}   # {global_id (int): SkyCoord}
        self.epoch_mappings   = []   # [{local_star_id: global_id}, ...]
        self.target_global_id = None

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

    def _find_lc_file(self, dirpath):
        lc = os.path.join(dirpath, self.LC_FILE)
        if not os.path.exists(lc):
            print(f"  [skip] {self.LC_FILE} not found in {dirpath}")
            return None
        return lc

    def _load_lightcurve(self, lc_path):
        """
        Parse calibrated_lightcurves.txt.

        Returns (DataFrame, ref_star_list).
        Columns include MJD, BJD, Exptim, secz, exposure_corr,
        mag_N, emag_N, RA_N, DEC_N for each aperture N.
        """
        ref_stars = []
        with open(lc_path) as fh:
            for line in fh:
                if line.startswith('# Ensemble reference stars:'):
                    ref_stars = [s.strip() for s in line.split(':')[1].split(',')]
                    break
        df = pd.read_csv(lc_path, comment='#')
        return df, ref_stars

    def _extract_star_catalog(self, df):
        """
        Build {local_star_id (str): SkyCoord} from RA_N / DEC_N columns.
        Uses the per-frame median to get a single position per star.
        """
        ra_cols = [c for c in df.columns if c.startswith('RA_')]
        catalog = {}
        for ra_col in ra_cols:
            star_id = ra_col[3:]          # strip 'RA_'
            dec_col = f'DEC_{star_id}'
            if dec_col not in df.columns:
                continue
            med_ra  = df[ra_col].median()
            med_dec = df[dec_col].median()
            if np.isfinite(med_ra) and np.isfinite(med_dec):
                catalog[star_id] = co.SkyCoord(
                    ra=med_ra * u.deg, dec=med_dec * u.deg, frame='icrs')
        return catalog

    def load_epochs(self):
        """Read all epoch files and extract per-star positions."""
        print(f"[MultiEpoch] Loading {len(self.epoch_dirs)} epoch directories...")
        for dirpath in self.epoch_dirs:
            lc_path = self._find_lc_file(dirpath)
            if lc_path is None:
                continue
            lc_df, ref_stars = self._load_lightcurve(lc_path)
            star_cat = self._extract_star_catalog(lc_df)
            if not star_cat:
                print(f"  [skip] No RA_N/DEC_N columns in {lc_path} — "
                      "re-run the per-night pipeline first.")
                continue
            label = '/'.join(lc_path.split(os.sep)[-4:-1])
            self.epochs.append({
                'label':     label,
                'lc_path':   lc_path,
                'lc_df':     lc_df,
                'star_cat':  star_cat,
                'ref_stars': ref_stars,
            })
            print(f"  Loaded: {label} | {len(lc_df)} frames | "
                  f"{len(star_cat)} stars with positions")
        print(f"[MultiEpoch] {len(self.epochs)} epochs loaded.")

    # ------------------------------------------------------------------
    # Cross-matching
    # ------------------------------------------------------------------

    def build_global_catalog(self):
        """
        Assign a unique global integer ID to every distinct sky position
        found across all epochs using astropy match_to_catalog_sky for
        robust cross-matching even when per-night WCS has systematic offsets.
        """
        print(f"[MultiEpoch] Cross-matching stars "
              f"(tolerance={self.tolerance_arcsec}\")...")
        next_gid = 0

        for ep_idx, epoch in enumerate(self.epochs):
            mapping  = {}
            local_ids   = list(epoch['star_cat'].keys())
            local_coords = co.SkyCoord(
                [epoch['star_cat'][lid].ra  for lid in local_ids],
                [epoch['star_cat'][lid].dec for lid in local_ids],
                frame='icrs')

            if not self.global_catalog:
                # Seed catalog from first epoch
                for lid, sky in zip(local_ids, local_coords):
                    self.global_catalog[next_gid] = sky
                    mapping[lid] = next_gid
                    next_gid += 1
            else:
                gids       = list(self.global_catalog.keys())
                cat_coords = co.SkyCoord(
                    [self.global_catalog[g].ra  for g in gids],
                    [self.global_catalog[g].dec for g in gids],
                    frame='icrs')

                # Match all stars at once (O(N log N) k-d tree)
                idx, sep2d, _ = local_coords.match_to_catalog_sky(cat_coords)

                for i, (lid, sky) in enumerate(zip(local_ids, local_coords)):
                    if sep2d[i].arcsec <= self.tolerance_arcsec:
                        matched_gid = gids[idx[i]]
                        # Running average of position
                        old = self.global_catalog[matched_gid]
                        self.global_catalog[matched_gid] = co.SkyCoord(
                            ra =(old.ra.deg  + sky.ra.deg)  / 2 * u.deg,
                            dec=(old.dec.deg + sky.dec.deg) / 2 * u.deg,
                            frame='icrs')
                        mapping[lid] = matched_gid
                    else:
                        self.global_catalog[next_gid] = sky
                        mapping[lid] = next_gid
                        next_gid += 1

            self.epoch_mappings.append(mapping)
            n_matched_existing = (
                int(np.sum(sep2d.arcsec <= self.tolerance_arcsec))
                if ep_idx > 0 else 0)
            n_new_stars = len(local_ids) - n_matched_existing
            print(f"  Epoch {ep_idx} ({epoch['label']}): "
                  f"{len(local_ids)} stars → "
                  f"matched existing: {n_matched_existing}  "
                  f"new: {n_new_stars}")

        print(f"[MultiEpoch] Global catalog: {len(self.global_catalog)} unique stars")

        # Identify target star
        if self.target_ra is not None and self.target_dec is not None:
            target = co.SkyCoord(ra=self.target_ra  * u.deg,
                                 dec=self.target_dec * u.deg, frame='icrs')
            gids       = list(self.global_catalog.keys())
            cat_coords = co.SkyCoord(
                [self.global_catalog[g].ra  for g in gids],
                [self.global_catalog[g].dec for g in gids],
                frame='icrs')
            sep  = target.separation(cat_coords).arcsec
            best = int(np.argmin(sep))
            self.target_global_id = gids[best]
            print(f"  Target star → global ID {self.target_global_id} "
                  f"(sep={sep[best]:.2f}\")")

    # ------------------------------------------------------------------
    # Merge
    # ------------------------------------------------------------------

    def build_merged_dataframe(self):
        """
        Combine all epochs into one wide DataFrame.

        Recovers per-frame instrumental magnitudes from the saved
        calibrated values:  instru_N = mag_N + exposure_corr
        (this reverses the per-night ensemble correction so the
        multi-epoch solver starts from scratch).

        Output columns:
            MJD, BJD, Exptim, secz, epoch_label, epoch_idx,
            instru_{gid}, instrue_{gid}  for every global star gid.
        """
        print("[MultiEpoch] Building merged DataFrame...")
        all_frames = []

        for ep_idx, (epoch, mapping) in enumerate(
                zip(self.epochs, self.epoch_mappings)):
            lc_df    = epoch['lc_df'].copy()
            exp_corr = (lc_df['exposure_corr'].values
                        if 'exposure_corr' in lc_df.columns
                        else np.zeros(len(lc_df)))

            frame = {}
            for col in ('MJD', 'BJD', 'Exptim', 'secz'):
                if col in lc_df.columns:
                    frame[col] = lc_df[col].values
            frame['epoch_label'] = epoch['label']
            frame['epoch_idx']   = ep_idx

            for mag_col in (c for c in lc_df.columns if c.startswith('mag_')):
                local_id = mag_col[4:]          # strip 'mag_'
                gid      = mapping.get(local_id)
                if gid is None:
                    continue
                emag_col = f'emag_{local_id}'
                frame[f'instru_{gid}']  = lc_df[mag_col].values + exp_corr
                if emag_col in lc_df.columns:
                    frame[f'instrue_{gid}'] = lc_df[emag_col].values

            all_frames.append(pd.DataFrame(frame))

        merged = pd.concat(all_frames, ignore_index=True)

        time_col = 'BJD' if 'BJD' in merged.columns else 'MJD'
        if time_col in merged.columns:
            merged.sort_values(time_col, inplace=True)
            merged.reset_index(drop=True, inplace=True)

        # Ensure every global star has a column (NaN where not detected)
        for gid in self.global_catalog:
            for prefix in ('instru_', 'instrue_'):
                col = f'{prefix}{gid}'
                if col not in merged.columns:
                    merged[col] = np.nan

        self.open_df = merged
        self.data    = merged.copy()
        self.df_keep = merged.copy()

        n_inst = len([c for c in merged.columns if c.startswith('instru_')])
        print(f"[MultiEpoch] Merged: {len(merged)} frames | "
              f"{n_inst} global stars")
        return merged

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def save_global_catalog(self):
        rows = [{'global_id': gid,
                 'ra_deg':    coord.ra.deg,
                 'dec_deg':   coord.dec.deg}
                for gid, coord in self.global_catalog.items()]
        path = os.path.join(self.save_path, 'global_star_catalog.txt')
        pd.DataFrame(rows).to_csv(path, index=False)
        print(f"Saved: {path}")

    def plot_star_map(self):
        """
        Sky map (RA vs Dec) of every star in the global catalog.

        - Target star      → gold star marker (large)
        - Surviving refs   → filled green circles
        - Other (rejected) → open grey circles
        - Each point is labelled with its global ID.
        RA axis is inverted (astronomical convention).
        """
        if not self.global_catalog:
            print("[StarMap] Global catalog empty — nothing to plot.")
            return

        gids = list(self.global_catalog.keys())
        ras  = np.array([self.global_catalog[g].ra.deg  for g in gids])
        decs = np.array([self.global_catalog[g].dec.deg for g in gids])

        target_gid = self.target_global_id
        surviving  = set(str(s) for s in getattr(self, 'surviving_stars', []))

        fig, ax = plt.subplots(figsize=(10, 9))

        # Categorise each star
        for gid, ra, dec in zip(gids, ras, decs):
            gid_str = str(gid)
            is_target = (target_gid is not None and gid_str == str(target_gid))
            is_ref    = gid_str in surviving and not is_target

            if is_target:
                ax.plot(ra, dec, marker='*', ms=22,
                        color='gold', markeredgecolor='black',
                        markeredgewidth=1.2, zorder=5,
                        label='Target' if 'Target' not in ax.get_legend_handles_labels()[1] else None)
                ax.annotate(f' {gid_str}', (ra, dec),
                            fontsize=10, fontweight='bold',
                            color='darkgoldenrod', zorder=6)
            elif is_ref:
                ax.plot(ra, dec, marker='o', ms=9,
                        color='seagreen', markeredgecolor='black',
                        markeredgewidth=0.6, zorder=4,
                        label=('Surviving reference'
                               if 'Surviving reference'
                               not in ax.get_legend_handles_labels()[1]
                               else None))
                ax.annotate(f' {gid_str}', (ra, dec),
                            fontsize=7, color='darkgreen', zorder=5)
            else:
                ax.plot(ra, dec, marker='o', ms=6,
                        markerfacecolor='none', markeredgecolor='grey',
                        markeredgewidth=0.8, alpha=0.7, zorder=2,
                        label=('Other / rejected'
                               if 'Other / rejected'
                               not in ax.get_legend_handles_labels()[1]
                               else None))
                ax.annotate(f' {gid_str}', (ra, dec),
                            fontsize=6, color='grey', alpha=0.7, zorder=3)

        # Field centre cross-hair
        cen_ra, cen_dec = float(np.mean(ras)), float(np.mean(decs))
        ax.plot(cen_ra, cen_dec, marker='+', ms=18,
                color='crimson', mew=1.5, zorder=1,
                label=f'Field centre  ({cen_ra:.4f}, {cen_dec:.4f})')

        ax.invert_xaxis()
        ax.set_xlabel('RA (deg)')
        ax.set_ylabel('Dec (deg)')
        ax.set_aspect('equal', adjustable='datalim')
        ax.grid(True, ls=':', alpha=0.5)

        n_total = len(gids)
        n_ref   = len(surviving)
        n_tgt   = 1 if target_gid is not None else 0
        ax.set_title(
            f'Global Star Map  |  {n_total} stars  '
            f'(target: {n_tgt}, surviving refs: {n_ref})',
            fontweight='bold', fontsize=12)
        ax.legend(loc='best', fontsize=9, framealpha=0.9)

        plt.tight_layout()
        out = os.path.join(self.save_path, 'star_map.png')
        plt.savefig(out, dpi=300, bbox_inches='tight')
        print(f"Saved: {out}")
        plt.show()

    def plot_epoch_overlay(self):
        """
        Overlay sky map showing each epoch's star positions in a
        different colour.  Useful for visualising WCS systematic
        offsets between nights and confirming the cross-matching.

        - Each epoch is one colour (legend = epoch label)
        - A star detected in multiple epochs appears multiple times
          (one point per epoch) at slightly offset positions.
        - The target star is overdrawn as a gold star marker.
        """
        if not self.epochs:
            print("[EpochOverlay] No epochs loaded — nothing to plot.")
            return

        fig, ax = plt.subplots(figsize=(11, 10))
        ep_colors = plt.cm.tab10(np.linspace(0, 1, max(len(self.epochs), 2)))

        all_ras, all_decs = [], []

        for ep_idx, (epoch, color) in enumerate(
                zip(self.epochs, ep_colors)):
            ras  = np.array([sk.ra.deg
                             for sk in epoch['star_cat'].values()])
            decs = np.array([sk.dec.deg
                             for sk in epoch['star_cat'].values()])
            local_ids = list(epoch['star_cat'].keys())
            mapping   = (self.epoch_mappings[ep_idx]
                         if ep_idx < len(self.epoch_mappings) else {})

            all_ras.extend(ras)
            all_decs.extend(decs)

            ax.plot(ras, decs, marker='o', ms=6,
                    color=color, markeredgecolor='black',
                    markeredgewidth=0.4, alpha=0.75, ls='none',
                    label=f'Epoch {ep_idx}: {epoch["label"]}  ({len(ras)} stars)',
                    zorder=3 + ep_idx)

            # Annotate each point with its global ID (small, same colour)
            for lid, ra, dec in zip(local_ids, ras, decs):
                gid = mapping.get(lid, '?')
                ax.annotate(f' {gid}', (ra, dec),
                            fontsize=5.5, color=color, alpha=0.8,
                            zorder=4 + ep_idx)

        # Target overlay (gold star)
        if self.target_global_id is not None:
            t_coord = self.global_catalog.get(self.target_global_id)
            if t_coord is not None:
                ax.plot(t_coord.ra.deg, t_coord.dec.deg,
                        marker='*', ms=24,
                        color='gold', markeredgecolor='black',
                        markeredgewidth=1.4, zorder=10,
                        label=f'Target (gid {self.target_global_id})')

        # Field centre
        cen_ra  = float(np.mean(all_ras))
        cen_dec = float(np.mean(all_decs))
        ax.plot(cen_ra, cen_dec, marker='+', ms=18,
                color='crimson', mew=1.5, zorder=1,
                label=f'Field centre  ({cen_ra:.4f}, {cen_dec:.4f})')

        ax.invert_xaxis()
        ax.set_xlabel('RA (deg)')
        ax.set_ylabel('Dec (deg)')
        ax.set_aspect('equal', adjustable='datalim')
        ax.grid(True, ls=':', alpha=0.5)
        ax.set_title(
            f'Per-Epoch Star Position Overlay  '
            f'(match tolerance = {self.tolerance_arcsec}″)',
            fontweight='bold', fontsize=12)
        ax.legend(loc='best', fontsize=8, framealpha=0.9)

        plt.tight_layout()
        out = os.path.join(self.save_path, 'epoch_overlay.png')
        plt.savefig(out, dpi=300, bbox_inches='tight')
        print(f"Saved: {out}")
        plt.show()

    def plot_epoch_maps(self):
        """
        Grid of sky maps — one subplot per epoch.

        Each panel shows:
          - All global-catalog positions as faint grey background
          - That epoch's stars in colour, labelled with global ID
          - The target star as a gold ★ if detected this epoch
        All panels share the same RA/Dec extent for direct comparison.
        """
        if not self.epochs:
            print("[EpochMaps] No epochs loaded — nothing to plot.")
            return

        cat_ra  = np.array([c.ra.deg  for c in self.global_catalog.values()])
        cat_dec = np.array([c.dec.deg for c in self.global_catalog.values()])
        ra_min, ra_max   = float(np.min(cat_ra)),  float(np.max(cat_ra))
        dec_min, dec_max = float(np.min(cat_dec)), float(np.max(cat_dec))
        ra_pad  = 0.05 * (ra_max - ra_min  + 1e-6)
        dec_pad = 0.05 * (dec_max - dec_min + 1e-6)

        n_ep   = len(self.epochs)
        ncols  = min(3, n_ep)
        nrows  = int(np.ceil(n_ep / ncols))
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(6 * ncols, 5.5 * nrows),
                                 sharex=True, sharey=True)
        axes = np.atleast_1d(axes).flatten()

        ep_colors = plt.cm.tab10(np.linspace(0, 1, max(n_ep, 2)))

        for ax, (ep_idx, epoch, color) in zip(
                axes, zip(range(n_ep), self.epochs, ep_colors)):

            ax.plot(cat_ra, cat_dec, marker='o', ms=4, ls='none',
                    markerfacecolor='none', markeredgecolor='lightgrey',
                    markeredgewidth=0.6, alpha=0.6, zorder=1)

            ras  = np.array([sk.ra.deg
                             for sk in epoch['star_cat'].values()])
            decs = np.array([sk.dec.deg
                             for sk in epoch['star_cat'].values()])
            local_ids = list(epoch['star_cat'].keys())
            mapping   = (self.epoch_mappings[ep_idx]
                         if ep_idx < len(self.epoch_mappings) else {})

            ax.plot(ras, decs, marker='o', ms=7, ls='none',
                    color=color, markeredgecolor='black',
                    markeredgewidth=0.4, alpha=0.85, zorder=3)

            for lid, ra, dec in zip(local_ids, ras, decs):
                gid = mapping.get(lid, '?')
                ax.annotate(f' {gid}', (ra, dec),
                            fontsize=6, color='black', alpha=0.9, zorder=4)

            if self.target_global_id is not None:
                tgt_local = [lid for lid, gid in mapping.items()
                             if str(gid) == str(self.target_global_id)]
                if tgt_local:
                    tc = epoch['star_cat'][tgt_local[0]]
                    ax.plot(tc.ra.deg, tc.dec.deg,
                            marker='*', ms=20,
                            color='gold', markeredgecolor='black',
                            markeredgewidth=1.2, zorder=6)

            ax.set_xlim(ra_max + ra_pad, ra_min - ra_pad)   # RA inverted
            ax.set_ylim(dec_min - dec_pad, dec_max + dec_pad)
            ax.set_aspect('equal', adjustable='box')
            ax.grid(True, ls=':', alpha=0.4)
            ax.set_title(
                f'Epoch {ep_idx}: {epoch["label"]}\n'
                f'{len(ras)} stars detected',
                fontsize=10, fontweight='bold', color=color)
            ax.set_xlabel('RA (deg)')
            ax.set_ylabel('Dec (deg)')

        for ax in axes[n_ep:]:
            fig.delaxes(ax)

        fig.suptitle(
            f'Per-Epoch Sky Maps  '
            f'(match tolerance = {self.tolerance_arcsec}″,  '
            f'gold ★ = target,  grey ○ = global catalog)',
            fontsize=13, fontweight='bold', y=1.01)
        plt.tight_layout()

        out = os.path.join(self.save_path, 'epoch_maps.png')
        plt.savefig(out, dpi=250, bbox_inches='tight')
        print(f"Saved: {out}")
        plt.show()

    def save_results_multiepoch(self):
        if 'exposure_corr' not in self.data.columns:
            print("Warning: run solve_ensemble() before saving.")
            return

        exp_corr  = self.data['exposure_corr'].reindex(self.df_keep.index)
        time_cols = [c for c in
                     ('MJD', 'BJD', 'Exptim', 'secz',
                      'epoch_label', 'epoch_idx')
                     if c in self.df_keep.columns]
        df_out = self.df_keep[time_cols].copy()

        new_cols = {'exposure_corr': exp_corr}
        for col in (c for c in self.df_keep.columns if c.startswith('instru_')):
            gid     = col.split('_', 1)[1]
            err_col = f'instrue_{gid}'
            new_cols[f'mag_{gid}']  = self.df_keep[col] - exp_corr
            if err_col in self.df_keep.columns:
                new_cols[f'emag_{gid}'] = self.df_keep[err_col]
            # carry RA/DEC from global catalog
            coord = self.global_catalog.get(int(gid) if gid.isdigit() else gid)
            if coord is not None:
                new_cols[f'RA_{gid}']  = coord.ra.deg
                new_cols[f'DEC_{gid}'] = coord.dec.deg

        df_out = pd.concat(
            [df_out, pd.DataFrame(new_cols, index=self.df_keep.index)], axis=1)

        path = os.path.join(self.save_path, 'multi_epoch_lightcurves.txt')
        with open(path, 'w') as fh:
            if hasattr(self, 'surviving_stars'):
                fh.write("# Multi-epoch ensemble reference stars (global IDs): "
                         + ", ".join(self.surviving_stars) + "\n")
            if self.target_global_id is not None:
                fh.write(f"# Target star global ID: {self.target_global_id}\n")
            df_out.to_csv(fh, index=False)
        print(f"Saved: {path}")

    @staticmethod
    def _fold_phase(t, period, t0):
        """Return phase in [0, 1) given a period and reference epoch."""
        return ((t - t0) / period) % 1.0

    def plot_multi_epoch_lightcurve(self, target_global_id=None,
                                    period=None, t0=None):
        """
        Two-panel plot:
          Top    — target star light curve coloured by epoch, with RA/DEC in title
          Bottom — all surviving reference stars (grey, median-centred)

        Parameters
        ----------
        period : float | None
            Orbital/rotational period in days.  When given the x-axis is
            plotted as phase = ((t - t0) / period) % 1  instead of time.
        t0 : float | None
            Reference epoch (same time system as the light curve, e.g. BJD).
            Defaults to the first data point when period is given.
        """
        if 'exposure_corr' not in self.data.columns:
            print("[MultiEpoch] No exposure_corr — ensemble not yet run.")
            return

        target_gid = (target_global_id if target_global_id is not None
                      else self.target_global_id)

        exp_corr = self.data['exposure_corr'].reindex(self.df_keep.index)
        time_col = 'BJD' if 'BJD' in self.df_keep.columns else 'MJD'
        t        = self.df_keep[time_col]

        # Phase-fold if requested
        if period is not None:
            t0_val  = t0 if t0 is not None else float(t.min())
            x       = self._fold_phase(t, period, t0_val)
            xlabel  = f'Phase  (P = {period:.6f} d,  T₀ = {t0_val:.4f})'
        else:
            x      = t
            xlabel = time_col

        epochs     = (self.df_keep['epoch_label'].unique()
                      if 'epoch_label' in self.df_keep.columns else ['all'])
        ep_colors  = plt.cm.tab20(np.linspace(0, 1, max(len(epochs), 1)))

        surviving  = getattr(self, 'surviving_stars', [])
        ref_colors = plt.cm.turbo(np.linspace(0, 1, max(len(surviving), 1)))

        has_target = (target_gid is not None and
                      f'instru_{target_gid}' in self.df_keep.columns)

        n_panels = 2 if has_target else 1
        fig, axes = plt.subplots(
            n_panels, 1,
            figsize=(16, 5 * n_panels),
            sharex=True,
            gridspec_kw={'height_ratios': [1.5, 1] if n_panels == 2 else [1]},
        )
        axes = np.atleast_1d(axes)
        ax_target = axes[0] if has_target else None
        ax_refs   = axes[1] if has_target else axes[0]

        # ── Top panel: target star coloured by epoch ──────────────────
        if has_target:
            instru_col  = f'instru_{target_gid}'
            instrue_col = f'instrue_{target_gid}'
            cal_mag     = self.df_keep[instru_col] - exp_corr

            # Get target RA/DEC from global catalog for the title
            t_int  = int(target_gid) if str(target_gid).isdigit() else target_gid
            t_coord = self.global_catalog.get(t_int)
            if t_coord is not None:
                ra_str  = f"RA={t_coord.ra.deg:.5f}°"
                dec_str = f"Dec={t_coord.dec.deg:.5f}°"
            else:
                ra_str  = f"RA={self.target_ra}°"  if self.target_ra  is not None else "RA=?"
                dec_str = f"Dec={self.target_dec}°" if self.target_dec is not None else "Dec=?"

            for ep, color in zip(epochs, ep_colors):
                if 'epoch_label' in self.df_keep.columns:
                    mask = self.df_keep['epoch_label'] == ep
                else:
                    mask = np.ones(len(self.df_keep), dtype=bool)
                x_ep   = x[mask]
                mag_ep = cal_mag[mask]
                err_ep = (self.df_keep.loc[mask, instrue_col]
                          if instrue_col in self.df_keep.columns else None)
                ax_target.errorbar(
                    x_ep, mag_ep, yerr=err_ep,
                    fmt='o', ms=3, color=color, ecolor=color,
                    alpha=0.7, elinewidth=0.8, capsize=0, label=str(ep))

            ax_target.invert_yaxis()
            ax_target.set_ylabel('Calibrated Mag')
            ax_target.set_title(
                f'Target Star  |  {ra_str}  {dec_str}  '
                f'(global ID {target_gid})',
                fontweight='bold', fontsize=12)
            ax_target.legend(loc='upper left', bbox_to_anchor=(1.01, 1),
                             fontsize=7, ncol=2, title='Epoch')
            ax_target.grid(True, ls=':', alpha=0.5)
        else:
            print("[MultiEpoch] No target star — plotting reference stars only.")

        # ── Bottom panel: surviving reference stars (median-centred) ──
        for sid, color in zip(surviving, ref_colors):
            col    = f'instru_{sid}'
            errcol = f'instrue_{sid}'
            if col not in self.df_keep.columns:
                continue
            s_int  = int(sid) if str(sid).isdigit() else sid
            s_coord = self.global_catalog.get(s_int)
            if s_coord is not None:
                lbl = (f'Ref {sid}  '
                       f'({s_coord.ra.deg:.4f}, {s_coord.dec.deg:.4f})')
            else:
                lbl = f'Ref {sid}'

            cal  = self.df_keep[col] - exp_corr
            cent = cal - np.nanmedian(cal)
            err  = (self.df_keep[errcol]
                    if errcol in self.df_keep.columns else None)
            ax_refs.errorbar(
                x, cent, yerr=err,
                fmt='o', ms=1.5, color=color, ecolor=color,
                alpha=0.45, elinewidth=0.5, capsize=0, label=lbl)

        ax_refs.axhline(0, color='black', lw=0.8, ls='--', alpha=0.5)
        ax_refs.invert_yaxis()
        ax_refs.set_xlabel(xlabel)
        ax_refs.set_ylabel('Δ Mag (median-centred)')
        ax_refs.set_title('Surviving Reference Stars', fontweight='bold')
        ax_refs.legend(loc='upper left', bbox_to_anchor=(1.01, 1),
                       fontsize=6, ncol=1)
        ax_refs.grid(True, ls=':', alpha=0.5)

        plt.tight_layout()
        suffix = f'_P{period:.6f}d' if period is not None else ''
        fname = (f'multi_epoch_lc_star_{target_gid}{suffix}.png' if has_target
                 else f'multi_epoch_reference_stars{suffix}.png')
        out = os.path.join(self.save_path, fname)
        plt.savefig(out, dpi=300, bbox_inches='tight')
        print(f"Saved: {out}")
        plt.show()

    # ------------------------------------------------------------------
    # SDSS photometric calibration
    # ------------------------------------------------------------------

    def calibrate_with_sdss(self, filter_band='i',
                             search_radius_arcsec=3.0,
                             data_release=18,
                             period=None, t0=None):
        """
        Query SDSS with ONE bulk request centred on the field, then
        cross-match locally to surviving reference stars.  Fit a linear
        regression (sdss_mag = slope × ens_mag + intercept) and apply
        it to all stars.

        Parameters
        ----------
        filter_band : str
            SDSS band: 'u', 'g', 'r', 'i', or 'z'.
        search_radius_arcsec : float
            Matching radius per star when cross-matching locally (default 3").
        data_release : int
            SDSS data release (default 18).
        """
        try:
            from astroquery.sdss import SDSS
        except ImportError:
            print("[SDSS] astroquery not installed: pip install astroquery")
            return

        if not hasattr(self, 'surviving_stars') or not self.surviving_stars:
            print("[SDSS] No surviving stars — run ensemble first.")
            return

        if 'exposure_corr' not in self.data.columns:
            print("[SDSS] exposure_corr missing — run ensemble first.")
            return

        exp_corr     = self.data['exposure_corr'].reindex(self.df_keep.index)
        query_fields = ['objid', 'ra', 'dec',
                        filter_band, f'err_{filter_band}']

        # ── ONE bulk query centred on the field ───────────────────────
        # Compute field centre and radius from global catalog positions
        all_ra  = np.array([c.ra.deg  for c in self.global_catalog.values()])
        all_dec = np.array([c.dec.deg for c in self.global_catalog.values()])
        cen_ra  = float(np.mean(all_ra))
        cen_dec = float(np.mean(all_dec))
        field_center = co.SkyCoord(ra=cen_ra * u.deg,
                                   dec=cen_dec * u.deg, frame='icrs')
        # Radius = largest separation from centre + matching buffer
        # SDSS caps single requests at 3 arcmin — tile if field is larger
        star_coords  = co.SkyCoord(ra=all_ra * u.deg,
                                   dec=all_dec * u.deg, frame='icrs')
        max_sep_arcmin = float(np.max(
            field_center.separation(star_coords).arcmin))
        query_radius_arcmin = max_sep_arcmin + search_radius_arcsec / 60.0
        SDSS_MAX_ARCMIN = 3.0

        import time

        def _query_with_retry(pos, radius_arcsec):
            """Single SDSS query with 3 retries and backoff."""
            for attempt in range(3):
                try:
                    return SDSS.query_region(
                        pos,
                        radius=radius_arcsec * u.arcsec,
                        photoobj_fields=query_fields,
                        data_release=data_release,
                    )
                except Exception as exc:
                    print(f"    attempt {attempt + 1}/3 failed: {exc}")
                    if attempt < 2:
                        time.sleep(5 * (attempt + 1))
            return None

        print(f"\n[SDSS] Field radius = {query_radius_arcmin:.2f} arcmin  "
              f"(SDSS cap = {SDSS_MAX_ARCMIN} arcmin)")

        from astropy.table import vstack as table_vstack

        if query_radius_arcmin <= SDSS_MAX_ARCMIN:
            # Field fits in one request
            print(f"[SDSS] Single bulk query: centre "
                  f"RA={cen_ra:.4f}° Dec={cen_dec:.4f}°  "
                  f"radius={query_radius_arcmin*60:.1f}\"  "
                  f"DR{data_release} {filter_band}-band...")
            sdss_table = _query_with_retry(
                field_center, query_radius_arcmin * 60.0)
        else:
            # Tile: query around each surviving star individually.
            # Use a wider QUERY radius than the MATCH radius so SDSS
            # objects just outside `search_radius_arcsec` still appear
            # as candidates — the final cross-match still enforces
            # `search_radius_arcsec`.
            per_star_query_arcsec = max(30.0, 4 * search_radius_arcsec)
            print(f"[SDSS] Field too large — querying each star "
                  f"(query r={per_star_query_arcsec:.0f}\", "
                  f"match r={search_radius_arcsec}\") ...")
            tables = []
            for gid_str in self.surviving_stars:
                gid   = int(gid_str) if gid_str.isdigit() else gid_str
                coord = self.global_catalog.get(gid)
                if coord is None:
                    continue
                pos = co.SkyCoord(ra=coord.ra.deg * u.deg,
                                  dec=coord.dec.deg * u.deg, frame='icrs')
                tbl = _query_with_retry(pos, per_star_query_arcsec)
                if tbl is not None:
                    tables.append(tbl)
                time.sleep(0.3)   # gentle rate limiting
            sdss_table = table_vstack(tables) if tables else None

        if sdss_table is None:
            print("[SDSS] All queries failed. Skipping SDSS calibration.")
            return

        # Remove duplicate SDSS objects (objid)
        if 'objid' in sdss_table.colnames:
            _, keep = np.unique(sdss_table['objid'], return_index=True)
            sdss_table = sdss_table[keep]

        print(f"  → {len(sdss_table)} unique SDSS objects retrieved.")

        sdss_coords = co.SkyCoord(
            ra=np.array(sdss_table['ra'])  * u.deg,
            dec=np.array(sdss_table['dec']) * u.deg,
            frame='icrs')
        sdss_mag_arr = np.array(sdss_table[filter_band],  dtype=float)
        sdss_err_arr = np.array(sdss_table[f'err_{filter_band}'], dtype=float)

        # ── Local cross-match: each surviving star → nearest SDSS object ─
        ens_mags     = []
        sdss_mags    = []
        sdss_errs    = []
        matched_gids = []

        for gid_str in self.surviving_stars:
            gid   = int(gid_str) if gid_str.isdigit() else gid_str
            coord = self.global_catalog.get(gid)
            if coord is None:
                continue

            instru_col = f'instru_{gid_str}'
            if instru_col not in self.df_keep.columns:
                continue

            mean_ens = float(np.nanmedian(
                self.df_keep[instru_col] - exp_corr))
            if not np.isfinite(mean_ens):
                continue

            seps = coord.separation(sdss_coords).arcsec
            best = int(np.argmin(seps))
            if seps[best] > search_radius_arcsec:
                print(f"  Star {gid_str}: no SDSS match within "
                      f"{search_radius_arcsec}\" (closest {seps[best]:.2f}\")")
                continue

            smag = float(sdss_mag_arr[best])
            serr = float(sdss_err_arr[best])

            if not (10.0 < smag < 25.0):
                print(f"  Star {gid_str}: SDSS {filter_band}={smag:.3f} "
                      "out of range — skipped")
                continue

            ens_mags.append(mean_ens)
            sdss_mags.append(smag)
            sdss_errs.append(serr)
            matched_gids.append(gid_str)
            print(f"  Star {gid_str}: ens={mean_ens:.4f}  "
                  f"SDSS {filter_band}={smag:.4f}±{serr:.4f}  "
                  f"sep={seps[best]:.2f}\"")

        n = len(matched_gids)
        if n < 2:
            print(f"[SDSS] Only {n} match(es) — need ≥2 for regression. "
                  "Skipping calibration.")
            return

        # Linear regression: sdss_mag = slope * ens_mag + intercept
        ens_arr  = np.array(ens_mags)
        sdss_arr = np.array(sdss_mags)
        slope, intercept, r_value, _, std_err = stats.linregress(
            ens_arr, sdss_arr)

        print(f"\n[SDSS] Linear fit ({n} stars):")
        print(f"       {filter_band}_sdss = {slope:.4f} × ens_mag "
              f"+ {intercept:.4f}")
        print(f"       R² = {r_value**2:.4f}  |  std_err = {std_err:.4f}")

        self.sdss_fit = {
            'slope':     slope,
            'intercept': intercept,
            'r2':        r_value ** 2,
            'n_stars':   n,
            'filter':    filter_band,
        }

        # Apply fit to every star in df_keep (build as dict then concat once)
        new_sdss_cols = {}
        for col in (c for c in self.df_keep.columns
                    if c.startswith('instru_')):
            gid_str = col.split('_', 1)[1]
            err_col = f'instrue_{gid_str}'
            ens_cal = self.df_keep[col] - exp_corr
            new_sdss_cols[f'sdss_mag_{gid_str}'] = slope * ens_cal + intercept
            if err_col in self.df_keep.columns:
                new_sdss_cols[f'sdss_emag_{gid_str}'] = (
                    abs(slope) * self.df_keep[err_col])
        self.df_keep = pd.concat(
            [self.df_keep, pd.DataFrame(new_sdss_cols, index=self.df_keep.index)],
            axis=1)

        # Regression diagnostic plot
        self._plot_sdss_regression(
            ens_arr, sdss_arr, sdss_errs,
            matched_gids, slope, intercept, r_value ** 2, filter_band)

        # Save SDSS-calibrated output
        self._save_sdss_calibrated(filter_band)

        # Light-curve plot using SDSS-calibrated mags
        self._plot_sdss_lightcurves(filter_band, period=period, t0=t0)
        print("[SDSS] Done.\n")

    # ------------------------------------------------------------------
    def _plot_sdss_regression(self, ens_mags, sdss_mags, sdss_errs,
                               labels, slope, intercept, r2, filter_band):
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.errorbar(ens_mags, sdss_mags, yerr=sdss_errs,
                    fmt='o', ms=6, color='steelblue', ecolor='steelblue',
                    elinewidth=1.2, capsize=4, label='Reference stars')

        x_fit = np.linspace(min(ens_mags) - 0.1, max(ens_mags) + 0.1, 100)
        ax.plot(x_fit, slope * x_fit + intercept,
                color='crimson', lw=2,
                label=f'Fit: {slope:.3f}×x + {intercept:.3f}\nR²={r2:.4f}')

        for x, y, lbl in zip(ens_mags, sdss_mags, labels):
            ax.annotate(f' {lbl}', (x, y), fontsize=7, color='gray')

        ax.set_xlabel('Ensemble calibrated mag (instrumental)')
        ax.set_ylabel(f'SDSS {filter_band} (mag)')
        ax.set_title('SDSS Photometric Calibration — Linear Regression',
                     fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, ls=':', alpha=0.5)
        plt.tight_layout()

        out = os.path.join(self.save_path,
                           f'sdss_regression_{filter_band}.png')
        plt.savefig(out, dpi=300, bbox_inches='tight')
        print(f"Saved: {out}")
        plt.show()

    # ------------------------------------------------------------------
    def _save_sdss_calibrated(self, filter_band):
        time_cols = [c for c in
                     ('MJD', 'BJD', 'Exptim', 'secz',
                      'epoch_label', 'epoch_idx')
                     if c in self.df_keep.columns]

        new_cols = {
            'exposure_corr': self.data['exposure_corr'].reindex(
                self.df_keep.index)
        }
        for col in (c for c in self.df_keep.columns
                    if c.startswith('sdss_mag_')):
            gid_str  = col[9:]               # strip 'sdss_mag_'
            emag_col = f'sdss_emag_{gid_str}'
            new_cols[col] = self.df_keep[col]
            if emag_col in self.df_keep.columns:
                new_cols[emag_col] = self.df_keep[emag_col]
            coord = self.global_catalog.get(
                int(gid_str) if gid_str.isdigit() else gid_str)
            if coord is not None:
                new_cols[f'RA_{gid_str}']  = coord.ra.deg
                new_cols[f'DEC_{gid_str}'] = coord.dec.deg

        df_out = pd.concat(
            [self.df_keep[time_cols],
             pd.DataFrame(new_cols, index=self.df_keep.index)],
            axis=1)

        path = os.path.join(self.save_path,
                            'sdss_calibrated_lightcurves.txt')
        fit = self.sdss_fit
        with open(path, 'w') as fh:
            fh.write(f"# SDSS {filter_band}-band calibration: "
                     f"sdss_mag = {fit['slope']:.6f} * ens_mag "
                     f"+ {fit['intercept']:.6f}  "
                     f"R2={fit['r2']:.4f}  n={fit['n_stars']}\n")
            if self.target_global_id is not None:
                fh.write(f"# Target star global ID: "
                         f"{self.target_global_id}\n")
            df_out.to_csv(fh, index=False)
        print(f"Saved: {path}")

    # ------------------------------------------------------------------
    def _plot_sdss_lightcurves(self, filter_band, period=None, t0=None):
        """
        Two-panel plot using SDSS-calibrated magnitudes:
          Top    — target star  sdss_mag_{gid}  coloured by epoch, RA/DEC in title
          Bottom — all surviving reference stars sdss_mag_{sid} (median-centred)

        Parameters
        ----------
        period : float | None
            Fold period in days.  None = plot vs time.
        t0 : float | None
            Reference epoch for phase fold.  Defaults to first data point.
        """
        time_col   = 'BJD' if 'BJD' in self.df_keep.columns else 'MJD'
        t          = self.df_keep[time_col]

        if period is not None:
            t0_val = t0 if t0 is not None else float(t.min())
            x      = self._fold_phase(t, period, t0_val)
            xlabel = f'Phase  (P = {period:.6f} d,  T₀ = {t0_val:.4f})'
        else:
            x      = t
            xlabel = time_col

        epochs     = (self.df_keep['epoch_label'].unique()
                      if 'epoch_label' in self.df_keep.columns else ['all'])
        ep_colors  = plt.cm.tab20(np.linspace(0, 1, max(len(epochs), 1)))
        surviving  = getattr(self, 'surviving_stars', [])
        ref_colors = plt.cm.turbo(np.linspace(0, 1, max(len(surviving), 1)))

        target_gid    = self.target_global_id
        target_col    = f'sdss_mag_{target_gid}' if target_gid is not None else None
        has_target    = (target_col is not None and
                         target_col in self.df_keep.columns)

        n_panels = 2 if has_target else 1
        fig, axes = plt.subplots(
            n_panels, 1,
            figsize=(16, 5 * n_panels),
            sharex=True,
            gridspec_kw={'height_ratios': [1.5, 1] if n_panels == 2 else [1]},
        )
        axes    = np.atleast_1d(axes)
        ax_tgt  = axes[0] if has_target else None
        ax_refs = axes[1] if has_target else axes[0]

        # ── Top: target star ─────────────────────────────────────────
        if has_target:
            errcol = f'sdss_emag_{target_gid}'
            t_int  = int(target_gid) if str(target_gid).isdigit() else target_gid
            coord  = self.global_catalog.get(t_int)
            if coord is not None:
                ra_str  = f"RA={coord.ra.deg:.5f}°"
                dec_str = f"Dec={coord.dec.deg:.5f}°"
            else:
                ra_str  = f"RA={self.target_ra}°"  if self.target_ra  is not None else "RA=?"
                dec_str = f"Dec={self.target_dec}°" if self.target_dec is not None else "Dec=?"

            for ep, color in zip(epochs, ep_colors):
                mask = (self.df_keep['epoch_label'] == ep
                        if 'epoch_label' in self.df_keep.columns
                        else np.ones(len(self.df_keep), dtype=bool))
                err = (self.df_keep.loc[mask, errcol]
                       if errcol in self.df_keep.columns else None)
                ax_tgt.errorbar(
                    x[mask], self.df_keep.loc[mask, target_col],
                    yerr=err, fmt='o', ms=3,
                    color=color, ecolor=color,
                    alpha=0.7, elinewidth=0.8, capsize=0, label=str(ep))

            ax_tgt.invert_yaxis()
            ax_tgt.set_ylabel(f'SDSS {filter_band} (mag)')
            ax_tgt.set_title(
                f'Target Star  |  {ra_str}  {dec_str}  '
                f'(global ID {target_gid})  —  SDSS {filter_band}-band',
                fontweight='bold', fontsize=12)
            ax_tgt.legend(loc='upper left', bbox_to_anchor=(1.01, 1),
                          fontsize=7, ncol=2, title='Epoch')
            ax_tgt.grid(True, ls=':', alpha=0.5)

        # ── Bottom: reference stars ────────────────────────────────────
        for sid, color in zip(surviving, ref_colors):
            scol   = f'sdss_mag_{sid}'
            errcol = f'sdss_emag_{sid}'
            if scol not in self.df_keep.columns:
                continue
            s_int  = int(sid) if str(sid).isdigit() else sid
            coord  = self.global_catalog.get(s_int)
            lbl    = (f'Ref {sid}  ({coord.ra.deg:.4f}, {coord.dec.deg:.4f})'
                      if coord is not None else f'Ref {sid}')

            mag  = self.df_keep[scol]
            cent = mag - np.nanmedian(mag)
            err  = (self.df_keep[errcol]
                    if errcol in self.df_keep.columns else None)
            ax_refs.errorbar(
                x, cent, yerr=err,
                fmt='o', ms=1.5, color=color, ecolor=color,
                alpha=0.45, elinewidth=0.5, capsize=0, label=lbl)

        ax_refs.axhline(0, color='black', lw=0.8, ls='--', alpha=0.5)
        ax_refs.invert_yaxis()
        ax_refs.set_xlabel(xlabel)
        ax_refs.set_ylabel(f'Δ SDSS {filter_band} (median-centred)')
        ax_refs.set_title(
            f'Surviving Reference Stars — SDSS {filter_band}-band',
            fontweight='bold')
        ax_refs.legend(loc='upper left', bbox_to_anchor=(1.01, 1),
                       fontsize=6, ncol=1)
        ax_refs.grid(True, ls=':', alpha=0.5)

        plt.tight_layout()
        suffix = f'_P{period:.6f}d' if period is not None else ''
        out = os.path.join(self.save_path,
                           f'sdss_lightcurves_{filter_band}{suffix}.png')
        plt.savefig(out, dpi=300, bbox_inches='tight')
        print(f"Saved: {out}")
        plt.show()

        # One subplot per star
        self._plot_sdss_all_stars_grid(filter_band, period=period, t0=t0)

    # ------------------------------------------------------------------
    def _plot_sdss_all_stars_grid(self, filter_band, period=None, t0=None):
        """
        One subplot per star (target + all surviving reference stars).
        Target panel has a gold border and is placed first.
        Each panel shows the SDSS-calibrated LC coloured by epoch,
        with RA/DEC in the subplot title.

        Parameters
        ----------
        period : float | None
            Fold period in days.  None = plot vs time.
        t0 : float | None
            Reference epoch for phase fold.  Defaults to first data point.
        """
        time_col  = 'BJD' if 'BJD' in self.df_keep.columns else 'MJD'
        t         = self.df_keep[time_col]

        if period is not None:
            t0_val = t0 if t0 is not None else float(t.min())
            x      = self._fold_phase(t, period, t0_val)
            xlabel = f'Phase  (P = {period:.6f} d)'
        else:
            x      = t
            xlabel = time_col

        epochs    = (self.df_keep['epoch_label'].unique()
                     if 'epoch_label' in self.df_keep.columns else ['all'])
        ep_colors = plt.cm.tab20(np.linspace(0, 1, max(len(epochs), 1)))

        target_gid = self.target_global_id
        surviving  = getattr(self, 'surviving_stars', [])

        # Build ordered star list: target first, then reference stars
        all_gids = []
        if target_gid is not None:
            all_gids.append(str(target_gid))
        for sid in surviving:
            if str(sid) != str(target_gid):
                all_gids.append(str(sid))

        # Keep only stars that have sdss_mag_ columns
        all_gids = [g for g in all_gids
                    if f'sdss_mag_{g}' in self.df_keep.columns]

        if not all_gids:
            print("[SDSS] No sdss_mag_ columns found — skipping grid plot.")
            return

        ncols  = 4
        nrows  = int(np.ceil(len(all_gids) / ncols))
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(6 * ncols, 4 * nrows),
                                 sharex=True)
        axes = np.atleast_1d(axes).flatten()

        for ax, gid_str in zip(axes, all_gids):
            is_target = (str(gid_str) == str(target_gid))
            scol      = f'sdss_mag_{gid_str}'
            errcol    = f'sdss_emag_{gid_str}'

            # RA/DEC label from global catalog
            g_int  = int(gid_str) if gid_str.isdigit() else gid_str
            coord  = self.global_catalog.get(g_int)
            if coord is not None:
                pos_str = (f"RA={coord.ra.deg:.4f}°\n"
                           f"Dec={coord.dec.deg:.4f}°")
            else:
                pos_str = ""

            # Plot each epoch
            for ep, color in zip(epochs, ep_colors):
                mask = (self.df_keep['epoch_label'] == ep
                        if 'epoch_label' in self.df_keep.columns
                        else np.ones(len(self.df_keep), dtype=bool))
                err  = (self.df_keep.loc[mask, errcol]
                        if errcol in self.df_keep.columns else None)
                ax.errorbar(
                    x[mask],
                    self.df_keep.loc[mask, scol],
                    yerr=err,
                    fmt='o', ms=2,
                    color=color, ecolor=color,
                    alpha=0.6, elinewidth=0.7, capsize=0,
                    label=str(ep))

            ax.invert_yaxis()
            ax.grid(True, ls=':', alpha=0.5)
            ax.set_ylabel(f'SDSS {filter_band}', fontsize=8)

            # Highlight target panel
            if is_target:
                for spine in ax.spines.values():
                    spine.set_edgecolor('gold')
                    spine.set_linewidth(2.5)
                ax.set_facecolor('#fffbea')
                title = (f'TARGET  ID={gid_str}\n{pos_str}')
                ax.set_title(title, fontsize=8, fontweight='bold',
                             color='goldenrod')
            else:
                title = f'Ref {gid_str}\n{pos_str}'
                ax.set_title(title, fontsize=8)

        # Add epoch colour legend on the last used axis
        last_ax = axes[len(all_gids) - 1]
        handles = [plt.Line2D([0], [0], marker='o', color='w',
                              markerfacecolor=ep_colors[i], markersize=6,
                              label=str(ep))
                   for i, ep in enumerate(epochs)]
        last_ax.legend(handles=handles, loc='upper left',
                       bbox_to_anchor=(1.02, 1),
                       fontsize=6, title='Epoch', title_fontsize=7)

        # Remove empty axes
        for ax in axes[len(all_gids):]:
            fig.delaxes(ax)

        # Shared x-label on bottom row
        for ax in axes[(nrows - 1) * ncols: len(all_gids)]:
            ax.set_xlabel(xlabel, fontsize=8)

        phase_info = (f'  |  Phase-folded  P={period:.6f} d'
                      if period is not None else '')
        fig.suptitle(
            f'All Stars — SDSS {filter_band}-band Calibrated Light Curves'
            f'{phase_info}\n(gold border = target star)',
            fontsize=14, fontweight='bold', y=1.01)
        plt.tight_layout()

        suffix = f'_P{period:.6f}d' if period is not None else ''
        out = os.path.join(self.save_path,
                           f'sdss_all_stars_{filter_band}{suffix}.png')
        plt.savefig(out, dpi=200, bbox_inches='tight')
        print(f"Saved: {out}")
        plt.show()

    # ------------------------------------------------------------------
    # Main pipeline
    # ------------------------------------------------------------------

    def _filter_min_epochs(self, min_epochs):
        """
        Drop from self.data any reference star whose instru_ column is
        non-NaN in fewer than min_epochs distinct epoch_idx values.
        Stars dropped here are kept in self.df_keep so their light curves
        are still saved in the output file.
        """
        if 'epoch_idx' not in self.data.columns:
            print(f"  [min_epochs] No epoch_idx column — skipping filter.")
            return

        n_epochs_total = self.data['epoch_idx'].nunique()
        mag_cols = [c for c in self.data.columns if c.startswith('instru_')]
        drop_ids = []

        for col in mag_cols:
            gid_str   = col.split('_', 1)[1]
            # Count how many distinct epochs have at least one valid frame
            n_seen = (self.data.groupby('epoch_idx')[col]
                      .apply(lambda s: s.notna().any()).sum())
            if n_seen < min_epochs:
                drop_ids.append(gid_str)

        if drop_ids:
            for gid_str in drop_ids:
                self.data.drop(
                    columns=[f'instru_{gid_str}', f'instrue_{gid_str}'],
                    errors='ignore', inplace=True)
            print(f"  [min_epochs≥{min_epochs}] Dropped {len(drop_ids)} stars "
                  f"seen in <{min_epochs}/{n_epochs_total} epochs: "
                  f"{drop_ids}")
        else:
            print(f"  [min_epochs≥{min_epochs}] All stars seen in ≥{min_epochs} epochs — none dropped.")

    # ------------------------------------------------------------------
    def run_multi_epoch(self, target_rms=0.02, numstars=12,
                        min_epochs=None,
                        ignor_global_ids=None,
                        sdss_filter=None,
                        sdss_search_radius_arcsec=3.0,
                        sdss_data_release=18,
                        period=None,
                        t0=None):
        """
        Full multi-epoch ensemble pipeline.

        Parameters
        ----------
        target_rms : float
            Convergence criterion (mag RMS).
        numstars : int
            Minimum number of reference stars to keep.
        min_epochs : int | None
            Only use reference stars detected in at least this many
            distinct epochs.  None = use all stars (default).
            Set to len(epoch_dirs) to require stars in every epoch.
        ignor_global_ids : list[int|str] | None
            Global star IDs to exclude from the reference ensemble
            (e.g. the science target).  If None and target_global_id is
            known, it is excluded automatically.
        sdss_filter : str | None
            If given ('g', 'r', 'i', 'z'), query SDSS for surviving
            reference stars and fit a linear regression to convert
            ensemble mags → SDSS-calibrated mags.  None skips this step.
        sdss_search_radius_arcsec : float
            SDSS matching radius per star (default 3").
        sdss_data_release : int
            SDSS data release (default 18).
        period : float | None
            If given, all light-curve plots use phase = ((t - t0) / period) % 1
            on the x-axis instead of time.
        t0 : float | None
            Reference epoch for phase folding (same units as BJD/MJD).
            Defaults to the first data point when period is given.
        """
        print("=" * 60)
        print("[MultiEpoch] 1/5  Loading epoch files...")
        self.load_epochs()
        if not self.epochs:
            print("[MultiEpoch] No valid epochs found. Aborting.")
            return

        print("[MultiEpoch] 2/5  Cross-matching stars by RA/DEC...")
        self.build_global_catalog()

        print("[MultiEpoch] 3/5  Building merged DataFrame...")
        self.build_merged_dataframe()

        # Determine which global IDs to exclude from reference set
        exclude = set(str(x) for x in (ignor_global_ids or []))
        if self.target_global_id is not None and not ignor_global_ids:
            exclude.add(str(self.target_global_id))
        if exclude:
            print(f"[MultiEpoch]       Excluding global IDs from reference: {exclude}")
            for gid_str in exclude:
                self.data.drop(
                    columns=[f'instru_{gid_str}', f'instrue_{gid_str}'],
                    errors='ignore', inplace=True)

        # Drop stars not seen in enough epochs
        if min_epochs is not None:
            self._filter_min_epochs(min_epochs)

        print("[MultiEpoch] 4/5  Ensemble calibration loop "
              f"(target RMS={target_rms}, min stars={numstars})...")
        self.history  = []
        iteration     = 1
        residuals     = np.array([np.nan])

        while True:
            mag_cols = [c for c in self.data.columns if c.startswith('instru_')]
            if len(mag_cols) < 2:
                print("  → Fewer than 2 stars remaining — stopping.")
                break

            residuals    = self.solve_ensemble()
            current_rms  = float(np.sqrt(np.nanmean(residuals ** 2)))
            print(f"  [iter {iteration:03d}] {len(mag_cols)} stars | "
                  f"RMS={current_rms:.5f}")

            if current_rms <= target_rms:
                print(f"  → Converged (RMS={current_rms:.5f} ≤ {target_rms})")
                break
            if len(mag_cols) <= numstars:
                print(f"  → Minimum star count ({numstars}) reached.")
                break
            if iteration > 100:
                print("  → Max iterations (100) reached.")
                break

            worst_col = list(self.find_most_variable_star().keys())[0]
            star_id   = worst_col.split('_')[-1]
            err_col   = f'instrue_{star_id}'
            self.data.drop(columns=[worst_col, err_col],
                           errors='ignore', inplace=True)
            print(f"      → Reject global star {star_id}")
            self.history.append({'iteration': iteration,
                                 'rms': current_rms, 'reject': star_id})
            iteration += 1

        self.surviving_stars = [c.split('_')[-1]
                                 for c in self.data.columns
                                 if c.startswith('instru_')]
        final_rms = float(np.sqrt(np.nanmean(residuals**2)))
        rms_str = f"{final_rms:.5f}" if np.isfinite(final_rms) else "n/a"
        print(f"  → Final: {len(self.surviving_stars)} reference stars | RMS={rms_str}")

        n_steps = "6" if sdss_filter else "5"
        print(f"[MultiEpoch] 5/{n_steps}  Saving results...")
        self.save_global_catalog()
        self.save_results_multiepoch()
        self.plot_star_map()
        self.plot_epoch_overlay()
        self.plot_epoch_maps()
        self.plot_rms_history(target_rms=target_rms,
                              save_folder=self.save_path)
        self.plot_multi_epoch_lightcurve(period=period, t0=t0)

        if sdss_filter:
            print(f"[MultiEpoch] 6/{n_steps}  SDSS photometric calibration "
                  f"({sdss_filter}-band)...")
            self.calibrate_with_sdss(
                filter_band=sdss_filter,
                search_radius_arcsec=sdss_search_radius_arcsec,
                data_release=sdss_data_release,
                period=period,
                t0=t0,
            )

        print("=" * 60)
        print("[MultiEpoch] Done.")


# ---------------------------------------------------------------------------
# Convenience: auto-discover run directories from a pattern
# ---------------------------------------------------------------------------

def find_epoch_dirs(base_pattern):
    """
    Return a sorted list of run directories matching a glob pattern.

    Example
    -------
    dirs = find_epoch_dirs(
        '/MSSP/sittipong/reduce2/*/hcam_reduction_*_i/results/run*')
    """
    return sorted(glob.glob(base_pattern))
