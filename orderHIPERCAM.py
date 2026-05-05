"""
auto_hcam_v3
"""
import os
import math
import glob
import json
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.visualization import ZScaleInterval
from photutils.detection import DAOStarFinder

import hipercam as hcam
from hipercam import MCCD
import hipercam.scripts as scripts

BLUE = '\033[94m'
RED = '\033[91m'
RESET = '\033[0m'


class Hipercam_setup:
    """
    Parameters
    ----------
    source : str
        Data source:
          'hs' : HiPERCAM server
          'hl' : local HiPERCAM FITS file
          'us' : ULTRACAM server
          'ul' : local ULTRACAM .xml/.dat files
          'hf' : list of HiPERCAM hcm FITS-format files
    """

    def __init__(self, source, input_run_ul=None, ccd='1',
                 save_dir=None, raw_dir=None, diagnostics=False):
        self.bias = {'file': [], 'master': [], 'lis': []}
        self.dark = {'file': [], 'master': [], 'lis': []}
        self.flat = {'file': [], 'master': [], 'lis': []}
        self.data = {'file': [], 'lis': []}

        self.bias_dir = None
        self.dark_dir = None
        self.flat_dir = None

        self.raw_dir = raw_dir
        self.base_dir = save_dir
        self.source = source

        if self.base_dir:
            os.makedirs(self.base_dir, exist_ok=True)

        if source == 'ul':
            if input_run_ul is not None:
                print('Work with source : ul')
                self.input_run_ul = input_run_ul
                self.source_ul()
            else:
                print('Source : ul must provide --input_run_ul--')
                return

        self.reduction()

        self.diagnostics = diagnostics
        if self.diagnostics:
            self.figs_dir = os.path.join(self.base_dir, 'figs')
            os.makedirs(self.figs_dir, exist_ok=True)
            self.diagnostic()

    # ------------------------------------------------------------------
    def diagnostic(self):
        """Plot raw frames, master frames, and median trend for all calibrations."""
        print(f"{BLUE}--- Running Diagnostics: Calibration Analysis ---{RESET}")
        
        categories = [
            ('bias', self.bias, "Bias"),
            ('dark', self.dark, "Dark"),
            ('flat', self.flat, "Flat"),
            ('data', self.data, "Data") 
        ]
        
        for key, data_dict, label in categories:
            # print(key, data_dict,)
            if data_dict.get('file'):
                print(f" Plotting raw {label} frames & trend...")
                self.plot_run(data_dict['file'],
                              save_name=f"fig_{key}.png")
                              

            if data_dict.get('master'):
                print(f" Plotting Master {label}...")
                self.plot_run(data_dict['master'],
                              save_name=f"master_{key}.png",
                              show_median=True)
                              
        print(f"{BLUE}--- Diagnostics Complete! Files saved in: {self.figs_dir} ---{RESET}\n")
        
    # ------------------------------------------------------------------
    def plot_run(self, files, save_name, target_ccd='1', df=None,
                 aps_count=0, show_median=False):
        """
        Plot a grid of HCM images and an optional median-level trend.

        Parameters
        ----------
        files       : list of full paths to .hcm files
        save_name   : output filename for the grid image
        target_ccd  : CCD label to display (e.g. '1', '2', '3')
        df          : optional DataFrame with aperture coordinates
        aps_count   : number of apertures to overlay from df
        show_median : display sigma-clipped median/min/max in titles
        """
        if not files:
            return

        frame_numbers = []
        median_levels = []

        ncols = 5
        nrows = max(1, math.ceil(len(files) / ncols))

        fig_grid, axs = plt.subplots(
            nrows, ncols,
            figsize=(6 * ncols, 6 * nrows),
            gridspec_kw={'hspace': 0.25, 'wspace': 0.05},
            facecolor='white',
        )
        axs_flat = axs.flatten() if hasattr(axs, 'flatten') else [axs]
        zscale = ZScaleInterval()

        for i, fname in enumerate(files):
            ax = axs_flat[i]
            try:
                mccd = MCCD.read(fname)
                data = mccd[str(target_ccd)]['1'].data

                clipped = sigma_clip(data, sigma=3.0)
                median_val = np.ma.median(clipped)
                min_val = np.ma.min(clipped)
                max_val = np.ma.max(clipped)

                frame_numbers.append(i + 1)
                median_levels.append(median_val)

                vmin, vmax = zscale.get_limits(data)
                ax.imshow(data, origin='lower', cmap='gray_r',
                          vmin=vmin, vmax=vmax)

                title = f"[#{i+1}] {os.path.basename(fname)}"
                if show_median:
                    title += (f"\nMed: {median_val:.5f} | "
                              f"Min: {min_val:.5f} | Max: {max_val:.5f}")
                ax.set_title(title, fontsize=8)
                ax.axis('off')

                if df is not None and i < len(df):
                    row = df.iloc[[i]]
                    if not row.empty:
                        for ap in range(1, aps_count + 1):
                            col_x, col_y = f'x_{ap}', f'y_{ap}'
                            if col_x in row.columns and col_y in row.columns:
                                x_, y_ = row.iloc[0][col_x], row.iloc[0][col_y]
                                if np.isfinite(x_) and np.isfinite(y_):
                                    circle = patches.Circle(
                                        (x_ - 0.5, y_ - 0.5), 10,
                                        edgecolor='lime', facecolor='none', lw=1.0)
                                    ax.add_patch(circle)
                                    ax.text(x_, y_ + 15, f"{ap}",
                                            color='red', fontsize=6, ha='center')

            except Exception as e:
                ax.set_title("Read Error", color='red', fontsize=9)
                ax.axis('off')
                print(f"{RED}Error plotting {fname}: {e}{RESET}")

        for j in range(i + 1, len(axs_flat)):
            axs_flat[j].axis('off')

        grid_path = os.path.join(self.figs_dir, save_name)
        plt.savefig(grid_path, bbox_inches='tight', dpi=100)
        plt.close(fig_grid)

        if len(frame_numbers) > 1:
            plt.figure(figsize=(12, 5))
            plt.plot(frame_numbers, median_levels, 'o-',
                     markersize=4, color='tab:blue', label='Median ADU')
            plt.title(f"Median Level Trend: "
                      f"{os.path.basename(save_name).replace('.png', '')}")
            plt.xlabel("Frame Index")
            plt.ylabel("Counts (ADU)")
            plt.grid(True, alpha=0.3)
            plt.legend()
            trend_path = os.path.join(
                self.figs_dir, save_name.replace('.png', '_trend.png'))
            plt.savefig(trend_path, bbox_inches='tight')
            plt.close()
            print(f"   -> Saved {save_name} and trend plot.")
        else:
            print(f"   -> Saved {save_name}.")

    # ------------------------------------------------------------------
    def grab(self, work_dir, run_path, prefix_out, **kwargs):
        """Grab raw frames from a single run and rename with prefix."""
        os.makedirs(work_dir, exist_ok=True)
        os.chdir(work_dir)
        run = os.path.basename(run_path)
        args = [
            "nodef",
            f"source={self.source}",
            f"run={run_path}",
            f"first={kwargs.get('f1', 1)}",
            f"last={kwargs.get('f2', 0)}",
            f"ndigit={kwargs.get('ndigit', 3)}",
            f"trim={kwargs.get('trim', False)}",
            f"bias={kwargs.get('bias', 'none')}",
            f"dark={kwargs.get('dark', 'none')}",
            f"flat={kwargs.get('flat', 'none')}",
            "fmap=none",
        ]
        try:
            print(f"{BLUE}#### Grabbing {run} from {run_path} ####{RESET}")
            scripts.grab(args)
            raw_files = sorted(glob.glob(f"{run}_*.hcm"))
            full_paths = []
            for f in raw_files:
                new_name = f"{prefix_out}_{f}"
                os.rename(f, new_name)
                full_paths.append(os.path.abspath(new_name))

            if 'bias' in prefix_out:
                self.bias['file'].extend(full_paths)
            elif 'dark' in prefix_out:
                self.dark['file'].extend(full_paths)
            elif 'flat' in prefix_out:
                self.flat['file'].extend(full_paths)
            elif 'data' in prefix_out:
                self.data['file'].extend(full_paths)

        except Exception as e:
            print(f"{RED}Error on {run}: {e}{RESET}")

    # ------------------------------------------------------------------
    def makebias(self, lis_filename, output_name, sigma=3):
        """Combine bias frames into a master bias."""
        args = [
            "nodef", "source=hf",
            f"flist={lis_filename}",
            f"sigma={sigma}",
            "plot=no",
            f"output={output_name}",
        ]
        try:
            print(f"{BLUE}--- Making Master Bias: "
                  f"{os.path.basename(output_name)} ---{RESET}")
            scripts.makebias(args)
            self.bias['master'].append(
                os.path.abspath(f'{output_name}.hcm'))
        except Exception as e:
            print(f"{RED}Error processing bias: {e}{RESET}")

    # ------------------------------------------------------------------
    def makeflat(self, lis_filename, output_name, sigma=3,
                 chosen_bias='none', chosen_dark='none'):
        """Combine flat frames into a master flat."""
        args = [
            "nodef", "source=hf", "ngroup=5",
            f"flist={lis_filename}",
            f"sigma={sigma}",
            "lower=10000", "upper=35000",
            "plot=no",
            f"bias={chosen_bias}",
            f"dark={chosen_dark}",
            f"output={output_name}",
        ]
        try:
            full_master_path = os.path.abspath(f'{output_name}.hcm')
            if os.path.exists(full_master_path):
                os.remove(full_master_path)
                print(f"Removed existing: {os.path.basename(full_master_path)}")
            print(f"{BLUE}--- Making Master Flat: "
                  f"{os.path.basename(output_name)} ---{RESET}")
            scripts.makeflat(args)
            self.flat['master'].append(full_master_path)
        except Exception as e:
            print(f"{RED}Error processing flat: {e}{RESET}")

    # ------------------------------------------------------------------
    def source_ul(self):
        """Grab and calibrate bias and flat frames from ULTRACAM local files."""
        # --- Bias ---
        self.bias_dir = os.path.join(self.base_dir, 'bias')
        os.makedirs(self.bias_dir, exist_ok=True)

        bias_subset = {k: v for k, v in self.input_run_ul.items()
                       if k.startswith('bias_')}
        for key, value in bias_subset.items():
            if value is None:
                continue
            rawdir = value.get('rawdir') or self.raw_dir
            for run, f1, f2 in value.get('runs', []):
                self.grab(work_dir=self.bias_dir,
                          run_path=os.path.join(rawdir, run),
                          prefix_out=key, f1=f1, f2=f2)

            relevant = [f for f in self.bias['file']
                        if f"{key}_" in os.path.basename(f)]
            if relevant:
                lis_path = os.path.join(self.bias_dir, f"{key}.lis")
                with open(lis_path, 'w') as f:
                    f.writelines(p + '\n' for p in relevant)
                full_lis = os.path.abspath(lis_path)
                self.bias['lis'].append(full_lis)
                self.makebias(
                    lis_filename=full_lis,
                    output_name=os.path.join(self.bias_dir, f"master_{key}"))
            print(f"Processing {key} finished...\n")

        # --- Flat ---
        self.flat_dir = os.path.join(self.base_dir, 'flat')
        os.makedirs(self.flat_dir, exist_ok=True)

        flat_subset = {k: v for k, v in self.input_run_ul.items()
                       if k.startswith('flat_')}
        for key, value in flat_subset.items():
            if value is None:
                continue

            # Match bias to this flat: flat_data -> bias_flat_data,
            # then fallback to any bias_flat, then bias_data/bias_all
            chosen_bias = 'none'
            bias_key = key.replace('flat_', 'bias_flat_', 1)
            for b in self.bias['master']:
                if f'master_{bias_key}' in os.path.basename(b):
                    chosen_bias = b
                    break
            if chosen_bias == 'none':
                for b in self.bias['master']:
                    if 'master_bias_flat' in os.path.basename(b):
                        chosen_bias = b
                        break
            if chosen_bias == 'none':
                for b in self.bias['master']:
                    if ('master_bias_all' in os.path.basename(b) or
                            'master_bias_data' in os.path.basename(b)):
                        chosen_bias = b
                        break

            rawdir = value.get('rawdir') or self.raw_dir
            for run, f1, f2 in value.get('runs', []):
                self.grab(work_dir=self.flat_dir,
                          run_path=os.path.join(rawdir, run),
                          prefix_out=key, f1=f1, f2=f2,
                          bias=chosen_bias)

            relevant = [f for f in self.flat['file']
                        if f"{key}_" in os.path.basename(f)]
            if relevant:
                lis_path = os.path.join(self.flat_dir, f"{key}.lis")
                with open(lis_path, 'w') as f:
                    f.writelines(p + '\n' for p in relevant)
                full_lis = os.path.abspath(lis_path)
                self.flat['lis'].append(full_lis)
                self.makeflat(
                    lis_filename=full_lis,
                    output_name=os.path.join(self.flat_dir, f"master_{key}"))
            print(f"Processing {key} finished...\n")

    # ------------------------------------------------------------------
    def reduction(self):
        """Grab science data frames with bias and flat applied."""
        if self.source == 'ul':
    
            self.data_dir = os.path.join(self.base_dir, 'data')
            os.makedirs(self.data_dir, exist_ok=True)
            data_subset = {k: v for k, v in self.input_run_ul.items() if k in ['data', 'science', 'target'] or k.startswith('data')}
            
            if not data_subset:
                print(f"{RED}Warning: No 'data' key found in your input_run dictionary! Skipping reduction.{RESET}")
                return
                
    
            # data_subset = {k: v for k, v in self.input_run_ul.items()
            #                if k.startswith('data')}
            
            for key, value in data_subset.items():
                if value is None:
                    continue
                rawdir = value.get('rawdir') or self.raw_dir
    
                chosen_bias = 'none'
                for b in self.bias['master']:
                    if 'master_bias_data' in os.path.basename(b):
                        chosen_bias = b
                        break
    
                chosen_flat = 'none'
                for f in self.flat['master']:
                    if 'master_flat_data' in os.path.basename(f):
                        chosen_flat = f
                        break
    
                print(f" Use Bias: {os.path.basename(chosen_bias) if chosen_bias != 'none' else 'none'}")
                print(f" Use Flat: {os.path.basename(chosen_flat) if chosen_flat != 'none' else 'none'}")
    
                for run, f1, f2 in value.get('runs', []):
                    self.grab(work_dir=self.data_dir,
                              run_path=os.path.join(rawdir, run),
                              prefix_out=key, f1=f1, f2=f2,
                              bias=chosen_bias, flat=chosen_flat)
    
                relevant = [f for f in self.data['file']
                            if f"{key}_" in os.path.basename(f)]
                if relevant:
                    lis_path = os.path.join(self.data_dir, f"{key}.lis")
                    with open(lis_path, 'w') as f:
                        f.writelines(p + '\n' for p in relevant)
                    self.data['lis'].append(os.path.abspath(lis_path))


# ======================================================================

class Photometry:
    def __init__(self, base_dir, lis, diagnostics=True, single_run=True):
        self.diagnostics = diagnostics
        self.base_dir = base_dir
        self.lis = lis
        self.single_run = single_run
        self.data = []
        self.entries_to_process = []

    def solvwcs(self, log_file=None, ccd_label='1', win_label='1',
                scale_low=0.3, scale_high=1.5,
                ra_center=None, dec_center=None, radius=5.0,
                output_csv='wcs_radec.csv', frame_idx=3,
                astrometry_cache='astrometry_cache'):
        """
        Solve WCS from the .ape star catalog produced by setaper(), then
        transform every aperture (x, y) in the reduce .log to (RA, DEC)
        for all frames.

        Requires the ``astrometry`` Python package (pip install astrometry).

        Parameters
        ----------
        log_file : str, optional
            Path to reduce .log file. Defaults to ``base_dir/all.log``.
        ccd_label : str
            CCD label used in setaper() and reduce().
        win_label : str
            Window label (e.g. '1').
        scale_low, scale_high : float
            Plate-scale bounds in arcsec/pixel.
        ra_center, dec_center : float, optional
            Field-centre hint in degrees (speeds up solve considerably).
        radius : float
            Search radius in degrees when a centre hint is given.
        output_csv : str
            Output CSV filename saved inside ``base_dir``.
        frame_idx : int
            Index into the HCM file list used as the reference frame.
            Should match the frame used by setaper() so the .ape file
            is consistent.
        astrometry_cache : str
            Directory where astrometry.net index files are cached.

        Returns
        -------
        pandas.DataFrame or None
            Columns: frame, mjd, aperture, x_log, y_log, x_pix, y_pix,
            ra_deg, dec_deg.
        """
        import astrometry
        import pandas as pd

        # ── 1. Collect HCM file list ──────────────────────────────────
        hcm_files = []
        for lis_path in self.lis:
            if os.path.exists(lis_path):
                with open(lis_path) as fh:
                    hcm_files.extend(l.strip() for l in fh if l.strip())
        if not hcm_files:
            print(f"{RED}[solvwcs] No HCM files found in list files.{RESET}")
            return None

        ref_idx = min(frame_idx, len(hcm_files) - 1)
        ref_hcm = hcm_files[ref_idx]
        print(f"{BLUE}[solvwcs] Reference frame [{ref_idx}]: "
              f"{os.path.basename(ref_hcm)}{RESET}")

        # ── 2. Get binning from reference window ──────────────────────
        try:
            mccd = hcam.MCCD.read(ref_hcm)
            window = mccd[ccd_label][win_label]
            xbin = int(getattr(window, 'xbin', 1))
            ybin = int(getattr(window, 'ybin', 1))
        except Exception as e:
            print(f"{RED}[solvwcs] Cannot read reference frame: {e}{RESET}")
            return None

        # ── 3. Build star list from .ape file ─────────────────────────
        # .ape x/y  =  DAOStarFinder (0-indexed window px) × xbin
        # astrometry.Solver.solve(stars=...) expects 0-indexed pixel coords
        ape_path = (self.data[0]['out_ape']
                    if self.data else
                    os.path.splitext(ref_hcm)[0] + '.ape')
        if not os.path.exists(ape_path):
            print(f"{RED}[solvwcs] .ape file not found: {ape_path}. "
                  f"Run setaper() first.{RESET}")
            return None

        try:
            with open(ape_path) as fh:
                ape_json = json.load(fh)
            # layout: ["hipercam.MccdAper",
            #          [ccd, ["hipercam.CcdAper", [id, {...}], ...]]]
            ccd_block = ape_json[1][1]   # starts with "hipercam.CcdAper"
            stars = []
            for entry in ccd_block[1:]:  # skip the type-string element
                ap = entry[1]
                stars.append((ap['x'] / xbin, ap['y'] / ybin))
            print(f"[solvwcs] {len(stars)} stars loaded from .ape")
        except Exception as e:
            print(f"{RED}[solvwcs] Failed to parse .ape file: {e}{RESET}")
            return None

        # ── 4. Build position hint ────────────────────────────────────
        # Priority: manual ra_center/dec_center → HCM header → None
        if ra_center is None or dec_center is None:
            try:
                hdr = mccd.header
                # common HiPERCAM/ULTRASPEC header keys for telescope pointing
                for ra_key in ('RA', 'RADEG', 'TELRA', 'RA_TEL'):
                    if ra_key in hdr:
                        ra_val = hdr[ra_key]
                        # convert HH:MM:SS string to degrees if needed
                        if isinstance(ra_val, str) and ':' in ra_val:
                            from astropy.coordinates import Angle
                            import astropy.units as u
                            ra_val = Angle(ra_val, unit=u.hourangle).deg
                        ra_center = float(ra_val)
                        break
                for dec_key in ('DEC', 'DECDEG', 'TELDEC', 'DEC_TEL'):
                    if dec_key in hdr:
                        dec_val = hdr[dec_key]
                        if isinstance(dec_val, str) and ':' in dec_val:
                            from astropy.coordinates import Angle
                            import astropy.units as u
                            dec_val = Angle(dec_val, unit=u.deg).deg
                        dec_center = float(dec_val)
                        break
                if ra_center is not None and dec_center is not None:
                    print(f"[solvwcs] Position hint from HCM header: "
                          f"RA={ra_center:.4f}°  Dec={dec_center:.4f}°")
                else:
                    print("[solvwcs] No position hint found in header "
                          "— blind solve (slow).")
            except Exception:
                print("[solvwcs] Could not read position from header "
                      "— blind solve (slow).")

        size_hint = astrometry.SizeHint(
            lower_arcsec_per_pixel=scale_low,
            upper_arcsec_per_pixel=scale_high,
        )
        position_hint = None
        if ra_center is not None and dec_center is not None:
            position_hint = astrometry.PositionHint(
                ra_deg=ra_center,
                dec_deg=dec_center,
                radius_deg=radius,
            )

        print(f"{BLUE}[solvwcs] Running astrometry solver ...{RESET}")
        with astrometry.Solver(
            astrometry.series_5200.index_files(
                cache_directory=astrometry_cache,
                scales={4, 5, 6},
            )
            + astrometry.series_4200.index_files(
                cache_directory=astrometry_cache,
                scales={6, 7, 12},
            )
        ) as solver:
            solution = solver.solve(
                stars=stars,
                size_hint=size_hint,
                position_hint=position_hint,
                solution_parameters=astrometry.SolutionParameters(),
            )

        if not solution.has_match():
            print(f"{RED}[solvwcs] No WCS solution found. "
                  f"Try adjusting scale_low/scale_high or providing "
                  f"ra_center/dec_center hint.{RESET}")
            return None

        match = solution.best_match()
        wcs = match.wcs   # astropy WCS object, 0-indexed pixel convention
        print(f"[solvwcs] Solved:  RA={match.center_ra_deg:.5f}°  "
              f"Dec={match.center_dec_deg:.5f}°")

        # ── 5. Parse reduce log file ──────────────────────────────────
        if log_file is None:
            log_file = os.path.join(self.base_dir, 'all.log')
        if not os.path.exists(log_file):
            print(f"{RED}[solvwcs] Log file not found: {log_file}{RESET}")
            return None

        ccd_log = self._read_hlog(log_file, ccd_label)
        if ccd_log is None:
            return None

        # ── 6. Apply WCS → RA/DEC for every aperture, every frame ─────
        # log x/y = 0-indexed window pixel × xbin
        #   x_pix (0-indexed) = log_x / xbin
        # astropy WCS.pixel_to_world() takes 0-indexed pixels
        records = []
        for ap_label in sorted(ccd_log.keys(),
                                key=lambda k: int(k) if str(k).isdigit() else k):
            ts = ccd_log[ap_label]
            x_log = np.asarray(ts['x'], dtype=float)
            y_log = np.asarray(ts['y'], dtype=float)
            t_arr = np.asarray(ts['t'], dtype=float)
            x_pix = x_log / xbin
            y_pix = y_log / ybin

            try:
                sky     = wcs.pixel_to_world(x_pix, y_pix)
                ra_arr  = np.atleast_1d(sky.ra.deg)
                dec_arr = np.atleast_1d(sky.dec.deg)
            except Exception as exc:
                print(f"[solvwcs] WCS transform failed (ap {ap_label}): {exc}")
                ra_arr  = np.full(len(x_log), np.nan)
                dec_arr = np.full(len(x_log), np.nan)

            for i in range(len(x_log)):
                records.append({
                    'frame':    i + 1,
                    'mjd':      t_arr[i],
                    'aperture': int(ap_label) if str(ap_label).isdigit()
                                else ap_label,
                    'x_log':   x_log[i],
                    'y_log':   y_log[i],
                    'x_pix':   x_pix[i],
                    'y_pix':   y_pix[i],
                    'ra_deg':  ra_arr[i],
                    'dec_deg': dec_arr[i],
                })

        if not records:
            print(f"{RED}[solvwcs] No data extracted from log file.{RESET}")
            return None

        df = pd.DataFrame(records)
        out_path = os.path.join(self.base_dir, output_csv)
        df.to_csv(out_path, index=False)
        n_ap = df['aperture'].nunique()
        n_fr = df['frame'].nunique()
        print(f"{BLUE}[solvwcs] {n_ap} apertures × {n_fr} frames "
              f"→ {out_path}{RESET}")
        return df

    # ------------------------------------------------------------------
    def _read_hlog(self, log_file, ccd_label):
        """
        Parse a HiPERCAM reduce log file.

        Returns
        -------
        dict  {ap_label: {'t': ndarray, 'x': ndarray, 'y': ndarray}}
        or None on failure.
        """
        # ── attempt 1: hipercam.hlog.Hlog ────────────────────────────
        try:
            import importlib
            hlog_mod = importlib.import_module('hipercam.hlog')
            log_data = hlog_mod.Hlog.read(log_file)
            ccd = log_data[ccd_label]
            result = {
                ap: {'t': np.asarray(ts.t),
                     'x': np.asarray(ts.x),
                     'y': np.asarray(ts.y)}
                for ap, ts in ccd.items()
            }
            print(f"[solvwcs] Log parsed via hipercam.hlog "
                  f"({len(result)} apertures)")
            return result
        except Exception:
            pass

        # ── attempt 2: pandas fallback (ASCII log) ────────────────────
        try:
            import pandas as pd
            header_cols = None
            with open(log_file) as fh:
                for line in fh:
                    stripped = line.lstrip('#').strip()
                    if (line.startswith('#') and
                            ('MJD' in stripped or 'mjd' in stripped)):
                        header_cols = stripped.split()
                        break

            if header_cols is None:
                print(f"{RED}[solvwcs] Cannot identify column names "
                      f"in {log_file}{RESET}")
                return None

            df = pd.read_csv(log_file, comment='#', sep=r'\s+',
                             names=header_cols, engine='python')
            df.columns = [c.lower() for c in df.columns]

            mjd_col = next((c for c in df.columns
                            if c in ('mjd', 'tmid', 't')), df.columns[1])
            ap_ids = sorted(
                {c.split('_')[1] for c in df.columns
                 if c.startswith('x_') and c.split('_')[1].isdigit()},
                key=int)

            result = {}
            for ap in ap_ids:
                xc, yc = f'x_{ap}', f'y_{ap}'
                if xc in df.columns and yc in df.columns:
                    result[ap] = {
                        't': df[mjd_col].to_numpy(dtype=float),
                        'x': df[xc].to_numpy(dtype=float),
                        'y': df[yc].to_numpy(dtype=float),
                    }
            print(f"[solvwcs] Log parsed via pandas fallback "
                  f"({len(result)} apertures)")
            return result if result else None
        except Exception as exc:
            print(f"{RED}[solvwcs] Log parse failed: {exc}{RESET}")
            return None
        
    def setaper(self, ccd_label='1', win_label='1', SIGMA_THRESHOLD=1.5,
                output_plot="detection_labeled.png", SKIP_BRIGHTEST=10,
                MARGIN_LEFT=15, MARGIN_RIGHT=15, MARGIN_BOTTOM=5, MARGIN_TOP=27,
                R_TARG=None, R_SKY1=16, R_SKY2=24, frame=5, diagnostics=False,
                solve_wcs=False, **wcs_kwargs):
        """
        Detect stars and write a HiPERCAM aperture (.ape) file.

        In single_run mode one representative frame (index=frame) is used.
        In multi-run mode every file in the list gets its own .ape file.

        Parameters
        ----------
        solve_wcs : bool
            If True, call solvwcs() automatically after the .ape file is
            written.  Pass any solvwcs keyword arguments via **wcs_kwargs
            (e.g. ra_center=123.4, dec_center=45.6, scale_low=0.35).
        """
        hcm_files = []
        for lis_path in self.lis:
            if os.path.exists(lis_path):
                with open(lis_path, 'r') as f:
                    hcm_files.extend(
                        line.strip() for line in f if line.strip())
        # print(f"{RED}  {hcm_files}.{RESET}" )
        if len(hcm_files) == 0:
            print(f"{RED}Error: No HCM files found in list files.{RESET}")
            return

        if self.single_run:
            frame_idx = 3 #min(frame, len(hcm_files) - 1)
            print(frame_idx)
            self._detect_and_write_ape(
                hcm_files[frame_idx], ccd_label, win_label,
                SIGMA_THRESHOLD, SKIP_BRIGHTEST,
                MARGIN_LEFT, MARGIN_RIGHT, MARGIN_BOTTOM, MARGIN_TOP,
                R_TARG, R_SKY1, R_SKY2, diagnostics)
        else:
            for target_file in hcm_files:
                self._detect_and_write_ape(
                    target_file, ccd_label, win_label,
                    SIGMA_THRESHOLD, SKIP_BRIGHTEST,
                    MARGIN_LEFT, MARGIN_RIGHT, MARGIN_BOTTOM, MARGIN_TOP,
                    R_TARG, R_SKY1, R_SKY2, diagnostics)

        if solve_wcs:
            self.solvwcs(ccd_label=ccd_label, win_label=win_label,
                         **wcs_kwargs)

    def _detect_and_write_ape(self, target_file, ccd_label, win_label,
                               SIGMA_THRESHOLD, SKIP_BRIGHTEST,
                               MARGIN_LEFT, MARGIN_RIGHT, MARGIN_BOTTOM, MARGIN_TOP,
                               R_TARG, R_SKY1, R_SKY2, diagnostics):
        """Detect stars in one file, write a .ape file, and optionally plot the result."""
        if not os.path.exists(target_file):
            print(f"{RED}Error: HCM file '{target_file}' not found. Skipping.{RESET}")
            return

        try:
            mccd = hcam.MCCD.read(target_file)
            window = mccd[ccd_label][win_label]
            data = window.data
            binning = getattr(window, 'xbin', 1)
            self.binning = binning
            h, w = data.shape
        except Exception as e:
            print(f"{RED}Error reading '{target_file}': {e}. Skipping.{RESET}")
            return

        _, median, std = sigma_clipped_stats(data, sigma=3)
        daofind = DAOStarFinder(fwhm=4.0, threshold=SIGMA_THRESHOLD * std)
        sources = daofind(data - median)
        
        if not sources:
            print(f"{RED}No sources detected in {os.path.basename(target_file)}.{RESET}")
            return
        
        # 1. Apply spatial mask to filter out edge stars
        mask = ((sources['xcentroid'] > MARGIN_LEFT) &
                (sources['xcentroid'] < w - MARGIN_RIGHT) &
                (sources['ycentroid'] > MARGIN_BOTTOM) &
                (sources['ycentroid'] < h - MARGIN_TOP))
        sources = sources[mask]
        
        # Check if any sources survived the mask
        if len(sources) == 0:
            print(f"{RED}No sources left in {os.path.basename(target_file)} after applying margins.{RESET}")
            return

        # 2. Sort by flux (brightest first)
        sources.sort('flux')
        sources.reverse()
        
        # 3. Calculate FWHM
        sources['fwhm'] = np.sqrt(sources['npix'] / np.pi) * 2.35
        
        # 4. Skip the brightest stars
        if len(sources) > SKIP_BRIGHTEST:
            sources = sources[SKIP_BRIGHTEST:]
            
        self.num_aps = len(sources)

        # 5. Build the Aperture list
        ccd_ape = []
        for i, source in enumerate(sources):
            x, y, fwhm = source['xcentroid'], source['ycentroid'], source['fwhm']
            r_t = (R_TARG if R_TARG is not None else fwhm) * binning
            
            ccd_ape.append([str(i + 1), {
                "Comment": "hipercam.Aperture",
                "x": float(x) * binning,
                "y": float(y) * binning,
                "rtarg": r_t,
                "rsky1": R_SKY1 * binning,
                "rsky2": R_SKY2 * binning,
                "ref": (i == 0), # Brightest remaining star becomes the reference
                "compo": False,
                "mask": [],
                "extra": [],
                "link": "" if i == 0 else "1", # Link all others to the reference
            }])

        # Write to JSON .ape file
        ape_json = ["hipercam.MccdAper",
                    [ccd_label, ["hipercam.CcdAper"] + ccd_ape]]
        out_ape = os.path.splitext(target_file)[0] + '.ape'
        
        with open(out_ape, 'w') as f:
            json.dump(ape_json, f, indent=2)
            
        self.data.append({
            'target': target_file,
            'out_ape': out_ape,
            'sources': len(sources),
        })

        # ==========================================
        #  Plotting Block (when diagnostics=True)
        # ==========================================
        if diagnostics:
            
            fig, ax = plt.subplots(figsize=(10, 10))
            
            zscale = ZScaleInterval()
            vmin, vmax = zscale.get_limits(data)
            ax.imshow(data, origin='lower', cmap='Greys', vmin=vmin, vmax=vmax)
            
            # Draw Margin Rectangle
            rect_width = w - MARGIN_LEFT - MARGIN_RIGHT
            rect_height = h - MARGIN_BOTTOM - MARGIN_TOP
            rect = patches.Rectangle((MARGIN_LEFT, MARGIN_BOTTOM), 
                                     rect_width, rect_height,
                                     linewidth=2, edgecolor='yellow', facecolor='none', linestyle='--')
            ax.add_patch(rect)
    
            # Draw Apertures
            for i, source in enumerate(sources):
                x, y = source['xcentroid'], source['ycentroid']
                ap_number = i + 1
                r_plot = (R_TARG if R_TARG is not None else source['fwhm'])
                
                # Target Aperture
                ax.add_patch(patches.Circle((x, y), r_plot, edgecolor='lime', facecolor='none', lw=1.5))
                # Sky Annulus Inner
                ax.add_patch(patches.Circle((x, y), R_SKY1, edgecolor='red', facecolor='none', lw=1, linestyle='--'))
                # Sky Annulus Outer
                ax.add_patch(patches.Circle((x, y), R_SKY2, edgecolor='red', facecolor='none', lw=1, linestyle='--'))
                ax.text(x, y + R_SKY2 + 2, str(ap_number), color='cyan', fontsize=12, fontweight='bold', ha='center', va='bottom')
            
            plt.title(f"Apertures on {os.path.basename(target_file)}", fontsize=14)
            plt.axis('off')
            plt.tight_layout()
            
            # Save Figure
            base_name = os.path.basename(target_file).replace('.hcm', '_aper.png')
            figs_dir = os.path.join(self.base_dir, 'figs')
            os.makedirs(figs_dir, exist_ok=True)
            save_path = os.path.join(figs_dir, base_name)
            
            plt.savefig(save_path, bbox_inches='tight', dpi=150)
            plt.close(fig)
            print(f"Aperture Plot saved: {save_path}")
            
    # # ------------------------------------------------------------------
    def genred(self, instru='ultraspec-tnt', **kwargs):
        """Create a HiPERCAM reduction control file (.red) via genred."""
        if not self.data:
            print(f"{RED}Error: No apertures defined. Run setaper() first.{RESET}")
            return

        self.entries_to_process = (
            [self.data[0]] if self.single_run else self.data)

        for entry in self.entries_to_process:
            ape_file = entry.get('out_ape')
            if not ape_file:
                continue

            out_red = (os.path.join(os.path.dirname(ape_file), 'all.red')
                       if self.single_run
                       else os.path.splitext(ape_file)[0] + '.red')

            args = ["genred", ape_file, out_red,
                    "none", "none", "none", "none",
                    "0", "none", instru]
            try:
                scripts.genred(args)
                entry['out_red'] = out_red
                print(f"Successfully created: {os.path.basename(out_red)}")
            except Exception as e:
                print(f"{RED}Error running genred on {ape_file}: {e}{RESET}")

    # ------------------------------------------------------------------
    def modify_red(self, target_section=None, **kwargs):
        """Edit key=value pairs in the .red reduction file."""
        if not self.entries_to_process:
            print(f"{RED}Error: No reduction files. Run genred() first.{RESET}")
            return
        if not kwargs:
            return
        print(self.entries_to_process)
        for entry in self.entries_to_process:
            red_file = entry.get('out_red')
            
            if not red_file or not os.path.exists(red_file):
                continue

            print(f"{BLUE}#### Modifying {red_file} "
                  f"[{target_section or 'All'}] ####{RESET}")

            with open(red_file, 'r') as f:
                lines = f.readlines()

            new_lines = []
            modified_keys = set()
            current_section = None
            keys_pattern = '|'.join(re.escape(str(k)) for k in kwargs)
            pattern = rf"^\s*\b({keys_pattern})\b\s*="

            for line in lines:
                strip_line = line.strip()
                if strip_line.startswith('[') and strip_line.endswith(']'):
                    current_section = strip_line[1:-1].strip().lower()

                match = re.match(pattern, line)
                if match:
                    key = match.group(1)
                    target = target_section.lower() if target_section else None
                    if target is None or current_section == target:
                        comment = (" # " + line.split("#", 1)[1].strip()
                                   if "#" in line else "")
                        new_lines.append(f"{key} = {kwargs[key]}{comment}\n")
                        modified_keys.add(key)
                        continue
                new_lines.append(line)

            with open(red_file, 'w') as f:
                f.writelines(new_lines)

            for k in kwargs:
                status = "[Modified]" if k in modified_keys else "[Skip/Fail]"
                print(f"   {status} {k}")

    # ------------------------------------------------------------------
    def reduce(self, sep_mode=False, diagnostics=False):
        """Run the HiPERCAM reduce pipeline on the prepared .red file."""
        if not self.entries_to_process:
            print(f"{RED}Error: No reduction files. Run genred() first.{RESET}")
            return

        if self.single_run:
            if not self.lis:
                print(f"{RED}Error: No list files available for reduction.{RESET}")
                return
            red_file = self.entries_to_process[0].get('out_red')
            if not red_file or not os.path.exists(red_file):
                print(f"{RED}Error: Reduction file not found. Run genred() first.{RESET}")
                return

            os.chdir(self.base_dir)
            list_file = self.lis[0]
            log_file = "all.log"
            print(list_file)

            reduce_args = [
                "reduce", "source=hf",
                f"flist={list_file}", "trim=False",
                f"rfile={red_file}", f"log={log_file}",
                "tkeep=1", "lplot=False", "implot=False",
            ]
            try:
                print(f"{BLUE} Starting FULLY AUTOMATED Reduction...{RESET}")
                scripts.reduce(reduce_args)
                print(f" SUCCESS! Log saved to: {log_file}")
            except Exception as e:
                print(f"{RED} Reduction Failed: {e}{RESET}")
        else:
            print("Not implemented yet")
