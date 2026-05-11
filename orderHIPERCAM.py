"""
auto_hcam_v3
"""
import logging
import os
import math
import glob
import json
import re
from pathlib import Path
import pandas as pd
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

# Instrument plate scale arcsec/px unbinned — ULTRASPEC on TNT
ULTRASPEC_SCALE_ARCSEC_PX = 0.45


def read_hipercam_log(filename):
    """read .log then return DataFrame """
    if not os.path.exists(filename):
        return None, 0

    base_cols = ['CCD', 'nframe', 'MJD', 'MJDok', 'Exptim', 'mfwhm', 'mbeta']
    ap_cols_template = ['x', 'xe', 'y', 'ye', 'fwhm', 'fwhme', 'beta', 'betae', 
                        'counts', 'countse', 'sky', 'skye', 'nsky', 'nrej', 'cmax', 'flag']

    max_aperture = 0
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') and 'flag_' in line:
                parts = line.split()
                for part in parts:
                    if part.startswith('flag_'):
                        try:
                            ap = int(part.replace('flag_', ''))
                            if ap > max_aperture: max_aperture = ap
                        except: pass
            if not line.startswith('#'): break 

    all_cols = base_cols.copy()
    for i in range(1, max_aperture + 1):
        all_cols.extend([f"{col}_{i}" for col in ap_cols_template])

    df = pd.read_csv(filename, comment='#', sep=r'\s+', header=None, names=all_cols)
    return df, max_aperture

class Reduction:
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

    def __init__(self, source, save_dir=None, raw_dir=None, input_run_ul=None,
                 fix_pixel=False, sigma_low=5, sigma_high=8, diagnostics=False):
        self.bias = {'file': [], 'master': [], 'lis': []}
        self.dark = {'file': [], 'master': [], 'lis': []}
        self.flat = {'file': [], 'master': [], 'lis': []}
        self.data = {'file': [], 'lis': []}

        self.bias_dir = None
        self.dark_dir = None
        self.flat_dir = None

        self.source = source
        self.base_dir = save_dir
        self.raw_dir = raw_dir

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
        elif self.source in ['hs', 'hl', 'us', 'hf']:
            print(f"Source '{self.source}' is not supported yet!")
            return

        self.reduction()

        if fix_pixel:
            self.make_bad_pixel_mask(sigma_low=sigma_low, sigma_high=sigma_high)
            self.fix_bad_pixels(self.data['file'])

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
                if vmin >= vmax:
                    vmin, vmax = float(np.nanmin(data)), float(np.nanmax(data))
                if vmin >= vmax:
                    vmin, vmax = vmin - 1, vmax + 1
                ax.imshow(data, origin='lower', cmap='gray_r',
                          vmin=vmin, vmax=vmax)

                title = f"[#{i+1}] {os.path.basename(fname)} median: {median_val} "
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
        
        print(f" Use Bias: {kwargs.get('bias', 'none')}")
        print(f" Use Dark: {kwargs.get('dark', 'none')}")
        print(f" Use Flat: {kwargs.get('flat', 'none')}")
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
    def make_bad_pixel_mask(self, sigma_low=5, sigma_high=8,
                            margin_left=10, margin_right=10,
                            margin_bottom=10, margin_top=25):
        """
        Detect dead and hot pixels from the master flat using sigma clipping.
        Dead pixels: normalised response < 1 - sigma_low  * sigma
        Hot  pixels: normalised response > 1 + sigma_high * sigma
        Result stored in self.bad_pixel_masks[ccd_label][win_label].

        Margin parameters (pixels): if provided, also report how many bad pixels
        fall inside the science region vs. the excluded margin area.
        """
        if not self.flat['master']:
            print(f"{RED}No master flat found. Run makeflat() first.{RESET}")
            return

        flat_file = self.flat['master'][0]
        print(f"{BLUE}--- Detecting bad pixels from: {os.path.basename(flat_file)} ---{RESET}")
        mccd = MCCD.read(flat_file)
        self.bad_pixel_masks = {}

        for ccd_label, ccd in mccd.items():
            self.bad_pixel_masks[ccd_label] = {}
            for win_label, window in ccd.items():
                data = window.data.astype(float)
                h, w = data.shape
                _, median, std = sigma_clipped_stats(data, sigma=3)
                # Normalise to response map (1.0 = normal pixel)
                norm     = data / median if median != 0 else data
                norm_std = std  / median if median != 0 else std
                dead = norm < (1.0 - sigma_low  * norm_std)
                hot  = norm > (1.0 + sigma_high * norm_std)
                mask = dead | hot
                self.bad_pixel_masks[ccd_label][win_label] = mask

                n_dead = int(np.sum(dead))
                n_hot  = int(np.sum(hot))
                n_total = int(np.sum(mask))
                print(f"  CCD {ccd_label} Win {win_label}: "
                      f"{n_dead} dead, {n_hot} hot ({n_total} total) "
                      f"out of {h*w} pixels ({100*n_total/(h*w):.1f}%)")

                # Report science-region count if margins are provided
                if any([margin_left, margin_right, margin_bottom, margin_top]):
                    ml, mr = margin_left, margin_right
                    mb, mt = margin_bottom, margin_top
                    science_mask = np.zeros((h, w), dtype=bool)
                    y0 = mb
                    y1 = h - mt if h > mt else h
                    x0 = ml
                    x1 = w - mr if w > mr else w
                    science_mask[y0:y1, x0:x1] = True
                    n_sci   = int(np.sum(mask &  science_mask))
                    n_marg  = int(np.sum(mask & ~science_mask))
                    sci_pix = int(np.sum(science_mask))
                    print(f"    -> Science region ({x0}:{x1}, {y0}:{y1}) = {sci_pix} px: "
                          f"{n_sci} bad ({100*n_sci/sci_pix:.1f}%) | "
                          f"Margins: {n_marg} bad (not used for photometry)")

        print(f"{BLUE}--- Bad pixel mask ready ---{RESET}\n")

    # ------------------------------------------------------------------
    def fix_bad_pixels(self, files):
        """
        Replace each bad pixel with the median of its 8 surrounding neighbours.
        Overwrites the .hcm files in-place.
        Requires make_bad_pixel_mask() to have been called first.
        """
        if not hasattr(self, 'bad_pixel_masks'):
            print(f"{RED}No bad pixel mask. Run make_bad_pixel_mask() first.{RESET}")
            return

        from scipy.ndimage import generic_filter

        def _neighbor_median(x):
            # x is the 3x3 patch flattened; index 4 is the centre pixel
            neighbors = np.concatenate([x[:4], x[5:]])
            valid = neighbors[np.isfinite(neighbors)]
            return float(np.median(valid)) if len(valid) > 0 else x[4]

        n_files = len(files)
        print(f"{BLUE}--- Fixing bad pixels in {n_files} files ---{RESET}")
        for idx, file_path in enumerate(files, 1):
            try:
                mccd = MCCD.read(file_path)
                any_fixed = False
                for ccd_label, ccd in mccd.items():
                    masks = self.bad_pixel_masks.get(ccd_label, {})
                    for win_label, window in ccd.items():
                        mask = masks.get(win_label)
                        if mask is None or not np.any(mask):
                            continue
                        data = window.data.astype(float)
                        filled = generic_filter(data, _neighbor_median, size=3, mode='nearest')
                        data[mask] = filled[mask]
                        window.data = data
                        any_fixed = True
                if any_fixed:
                    mccd.write(file_path, overwrite=True)
                print(f"  [{idx}/{n_files}] Fixed: {os.path.basename(file_path)}")
            except Exception as e:
                print(f"{RED}  Error fixing {os.path.basename(file_path)}: {e}{RESET}")

        print(f"{BLUE}--- Bad pixel correction complete ---{RESET}\n")

    # ------------------------------------------------------------------
    def source_ul(self):
        """Grab and calibrate bias and flat frames from ULTRACAM local files."""
        # --- Bias ---
        self.bias_dir = os.path.join(self.base_dir, 'bias')
        os.makedirs(self.bias_dir, exist_ok=True)

        bias_subset = {k: v for k, v in self.input_run_ul.items()
                       if k.startswith('bias_')}
        for key, value in bias_subset.items():    #bias_data, bias_dark, bias_flat
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
        for key, value in flat_subset.items():   # 'flat_data'
            if value is None:
                continue
                
            # Match bias to this flat: flat_data -> bias_flat_data,
            # then fallback to any bias_flat, then bias_data/bias_all 
            chosen_bias = 'none'
            # Priority 1: Look for an a general 'master_bias_flat'
            if chosen_bias == 'none':
                for b in self.bias['master']:
                    if 'master_bias_flat' in os.path.basename(b):
                        chosen_bias = b
                        break

            # Priority 1: Fallback to 'master_bias_data' or 'master_bias_all'
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
    def __init__(self, base_dir, lis, diagnostics=True, bad_pixel_masks=None):
        self.diagnostics = diagnostics
        self.base_dir = base_dir
        self.lis = lis
        self.bad_pixel_masks = bad_pixel_masks  # {ccd_label: {win_label: bool_array}}
        self.entries_to_process = []
        self.data = {'lis':[], 'hcam_file': [], 'ape': [], 'red_file': [], 'log_file': [] }
        self._read_lis()

    def _read_lis(self):
        for lis in self.lis:
            # key = os.path.splitext(os.basename(lis))[0]
            self.data['lis'].append(lis)
            if os.path.exists(lis):
                hcm_files = []
                with open(lis, 'r') as f:
                    hcm_files.extend(
                        line.strip() for line in f if line.strip())
                    self.data['hcam_file'].append(hcm_files)
                    
    def setaper(self, ccd_label='1', win_label='1', SIGMA_THRESHOLD=3, SKIP_BRIGHTEST=10,
                MARGIN_LEFT=15, MARGIN_RIGHT=15, MARGIN_BOTTOM=5, MARGIN_TOP=27, ref_index=1,
                R_TARG=None, R_SKY1=16, R_SKY2=24, frame=2, SATURATION=50000, diagnostics=False):
        """
        Detect stars and write a HiPERCAM aperture (.ape) file.
        """

        for hcm_files in self.data['hcam_file']:
            self._detect_and_write_ape(
                hcm_files[frame], ccd_label, win_label, ref_index,
                SIGMA_THRESHOLD, SKIP_BRIGHTEST,
                MARGIN_LEFT, MARGIN_RIGHT, MARGIN_BOTTOM, MARGIN_TOP,
                R_TARG, R_SKY1, R_SKY2, SATURATION, diagnostics)

    def _detect_and_write_ape(self, target_file, ccd_label, win_label, ref_index,
                               SIGMA_THRESHOLD, SKIP_BRIGHTEST,
                               MARGIN_LEFT, MARGIN_RIGHT, MARGIN_BOTTOM, MARGIN_TOP,
                               R_TARG, R_SKY1, R_SKY2, SATURATION, diagnostics):
        """Detect stars in one file, write a .ape file, and optionally plot the result."""
        if not os.path.exists(target_file):
            print(f"{RED}Error: HCM file '{target_file}' not found. Skipping.{RESET}")
            return

        try:
            mccd = hcam.MCCD.read(target_file)
            window = mccd[ccd_label][win_label]
            data = window.data
            binning = getattr(window, 'xbin', 1)
            llx = getattr(window, 'llx', 0)  # 1-indexed full-CCD lower-left x
            lly = getattr(window, 'lly', 0)  # 1-indexed full-CCD lower-left y
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

        # 0. Remove saturated stars (peak above background + median >= SATURATION)
        sat_mask = (sources['peak'] + median) < SATURATION
        n_sat = np.sum(~sat_mask)
        if n_sat > 0:
            print(f" Removed {n_sat} saturated star(s) (peak >= {SATURATION} ADU)")
        sources = sources[sat_mask]
        if len(sources) == 0:
            print(f"{RED}No unsaturated sources in {os.path.basename(target_file)}.{RESET}")
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

        # ref_index is 1-based; clamp then convert to 0-based
        ref_idx    = max(1, min(ref_index, len(sources))) - 1
        ref_ap_num = str(ref_idx + 1)
        print(f" Reference aperture: {ref_ap_num} "
              f"(brightness rank {ref_idx + 1 + SKIP_BRIGHTEST} overall)")

        # 5. Build the Aperture list
        ccd_ape = []
        for i, source in enumerate(sources):
            x, y, fwhm = source['xcentroid'], source['ycentroid'], source['fwhm']
            r_t = (R_TARG if R_TARG is not None else fwhm) * binning
            is_ref = (i == ref_idx)
            ccd_ape.append([str(i + 1), {
                "Comment": "hipercam.Aperture",
                "x": float(x) * binning,
                "y": float(y) * binning,
                "rtarg": r_t,
                "rsky1": R_SKY1 * binning,
                "rsky2": R_SKY2 * binning,
                "ref": is_ref,
                "compo": False,
                "mask": [],
                "extra": [],
                "link": "" if is_ref else ref_ap_num,
            }])

        # Write to JSON .ape file
        ape_json = ["hipercam.MccdAper",
                    [ccd_label, ["hipercam.CcdAper"] + ccd_ape]]
        out_ape = os.path.splitext(target_file)[0] + '.ape'
        
        with open(out_ape, 'w') as f:
            json.dump(ape_json, f, indent=2)
            
        self.data['ape'].append(out_ape)

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
                x, y = source['xcentroid'], source['ycentroid']   # window-pixel coords
                ap_number = i + 1
                is_ref_ap = (i == ref_idx)
                r_plot = (R_TARG if R_TARG is not None else source['fwhm'])
                edge_color = 'red' if is_ref_ap else 'lime'
                lw = 2.5 if is_ref_ap else 1.5

                ax.add_patch(patches.Circle((x, y), r_plot,
                                            edgecolor=edge_color, facecolor='none', lw=lw))
                ax.add_patch(patches.Circle((x, y), R_SKY1,
                                            edgecolor='orange', facecolor='none',
                                            lw=1.0, linestyle='--', alpha=0.7))
                ax.add_patch(patches.Circle((x, y), R_SKY2,
                                            edgecolor='deepskyblue', facecolor='none',
                                            lw=1.0, linestyle='--', alpha=0.7))
                label = f"[REF]{ap_number}" if is_ref_ap else str(ap_number)
                ax.text(x, y + R_SKY2 + 2, label,
                        color='red' if is_ref_ap else 'cyan',
                        fontsize=10 if is_ref_ap else 8,
                        fontweight='bold', ha='center', va='bottom')
            plt.title(f"Apertures on {os.path.basename(target_file)}", fontsize=14)
            plt.axis('off')
            plt.tight_layout()
            
            # Save Figure
            base_name = os.path.basename(target_file).replace('.hcm', '_aper.png')
            figs_dir = os.path.join(self.base_dir)
            # os.makedirs(figs_dir, exist_ok=True)
            save_path = os.path.join(figs_dir, base_name)
            
            plt.savefig(save_path, bbox_inches='tight', dpi=150)
            plt.close(fig)
            print(f"Aperture Plot saved: {save_path}")
            
    # # ------------------------------------------------------------------
    def genred(self, instru='ultraspec-tnt', **kwargs):
        """Create a HiPERCAM reduction control file (.red) via genred."""
        if not self.data['ape']:
            print(f"{RED}Error: No apertures defined. Run setaper() first.{RESET}")
            return

        for lis, ape in zip(self.data['lis'], self.data['ape']):
            red_name = os.path.basename(lis).replace('.lis', '.red')
            out_red = os.path.join(self.base_dir, red_name)

            args = ["genred", ape, out_red,
                    "none", "none", "none", "none",
                    "0", "none", instru]
            try:
                scripts.genred(args)
                self.data['red_file'].append(out_red)
                print(f"[Genred] Successfully created: {red_name}")
            except Exception as e:
                print(f"{RED}[Genred] Error running genred on {ape}: {e}{RESET}")

    # ------------------------------------------------------------------
    def modify_red(self, target_section=None, **kwargs):
        """Edit key=value pairs in the .red reduction file."""
        if not self.data['red_file']:
            print(f"{RED}[[Modified] Error: No reduction files. Run genred() first.{RESET}")
            return
            
        if not kwargs:
            return
        # print(self.entries_to_process)
        for red_file in self.data['red_file']:
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
    def reduce(self, R_EXTARCT=None, plot_all=False, plot_zoom=False, plot_with_log=True, plot_with_zoom=True, diagnostics=False):
        """Run the HiPERCAM reduce pipeline on the prepared .red file."""
        if not self.data['red_file']:
            print(f"{RED}Error: No reduction files. Run genred() first.{RESET}")
            return

        for red_file, list_file in zip(self.data['red_file'], self.data['lis']):
            if not os.path.exists(red_file):
                print(f"{RED}Error: '{red_file}' not found. Skipping.{RESET}")
                continue

            log_file = os.path.splitext(os.path.basename(red_file))[0] + '.log'
            full_log = os.path.join(self.base_dir, log_file)

            os.chdir(self.base_dir)
            reduce_args = [
                "reduce", "source=hf",
                f"flist={list_file}", "trim=False",
                f"rfile={red_file}", f"log={log_file}",
                "tkeep=1", "lplot=False", "implot=False",
            ]
            try:
                print(f"{BLUE} Starting Reduction: {red_file}...{RESET}")
                scripts.reduce(reduce_args)
                self.data['log_file'].append(full_log)
                print(f" SUCCESS! Log saved to: {log_file}")
            except Exception as e:
                print(f"{RED} Reduction Failed: {e}{RESET}")
                continue

            if plot_with_log:
                print(f"{BLUE}--- Diagnostics: running plot_with_log ---{RESET}")
                self.plot_with_log(full_log)
            if plot_with_zoom:
                print(f"{BLUE}--- Diagnostics: running plot_with_zoom ---{RESET}")
                self.plot_with_zoom(full_log, R_EXTARCT=R_EXTARCT)

    # ------------------------------------------------------------------
    def _overlay_bad_pixels(self, ax, ccd_label, win_label):
        """Overlay bad pixels as a semi-transparent red mask on ax."""
        if self.bad_pixel_masks is None:
            return
        mask = self.bad_pixel_masks.get(ccd_label, {}).get(win_label)
        if mask is None or not np.any(mask):
            return
        overlay = np.zeros((*mask.shape, 4), dtype=float)  # RGBA
        overlay[mask] = [1.0, 0.0, 0.0, 0.55]             # red, 55% opacity
        ax.imshow(overlay, origin='lower', aspect='auto', interpolation='nearest')

    # ------------------------------------------------------------------
    def plot_with_zoom(self, log_filename, ccd_num='1', zoom_box=15, R_EXTARCT=None):
        """
        Plots a zoomed-in cutout of EVERY aperture for EACH frame.
        Includes X and Y marginal profile histograms for each star.
        """
        zoom_dir = os.path.join(self.base_dir, 'zoom')
        os.makedirs(zoom_dir, exist_ok=True)

        from mpl_toolkits.axes_grid1 import make_axes_locatable

        df_log, n_aps = read_hipercam_log(log_filename)
        if df_log is None:
            print(f"{RED}Error: Cannot read log file.{RESET}"); return

        lis_file = self.lis[0] if self.lis else None
        if not lis_file or not os.path.exists(lis_file):
            print(f"{RED}Error: List file not found!{RESET}"); return

        with open(lis_file, 'r') as f:
            hcm_files = [line.strip() for line in f.readlines() if line.strip()]

        n_plot = len(hcm_files)
        print(f"{BLUE}--- Generating zoomed plots with X/Y Profiles for {n_plot} frames ---{RESET}")

        zscale = ZScaleInterval()

        for i in range(n_plot):
            file_path = hcm_files[i]
            frame_name = os.path.basename(file_path)
            print(f'Generating zoomed plots F{i} {frame_name}')
            try:
                mccd = hcam.MCCD.read(file_path)
                win_obj = mccd[ccd_num]['1']
                data = win_obj.data
                llx = getattr(win_obj, 'llx', 1)
                lly = getattr(win_obj, 'lly', 1)
                vmin, vmax = zscale.get_limits(data)

                if i >= len(df_log):
                    print(f"{RED}Warning: No log data for frame {i}. Skipping.{RESET}")
                    continue

                frame_data = df_log.iloc[[i]]
                if frame_data.empty:
                    continue

                ncols = 5
                nrows = math.ceil(n_aps / ncols)
                fig, axs = plt.subplots(nrows, ncols, figsize=(4.5*ncols, 4.5*nrows),
                                        gridspec_kw={'wspace': 0.4, 'hspace': 0.4}, facecolor='whitesmoke')
                axs_flat = np.atleast_1d(axs).flatten()

                _bin = getattr(self, 'binning', 1)
                _E = R_EXTARCT if R_EXTARCT is not None else [1.8, 3.0, 15.0, 2.5, 15.0, 18.0, 3.0, 18.0, 20.0]
                mfwhm_unbinned = frame_data['mfwhm'].values[0]
                def _r(scale, rmin, rmax):
                    return float(np.clip(scale * mfwhm_unbinned, rmin, rmax)) / _bin if pd.notna(mfwhm_unbinned) else 5.0
                r_targ = _r(_E[0], _E[1], _E[2])
                r_sky1 = _r(_E[3], _E[4], _E[5])
                r_sky2 = _r(_E[6], _E[7], _E[8])
                fig.suptitle(
                    f"Frame {i+1}: {frame_name} | "
                    f"rtarg={r_targ*_bin:.2f} rsky1={r_sky1*_bin:.2f} rsky2={r_sky2*_bin:.2f} (unbin)",
                    fontsize=14, fontweight='bold')

                for ap in range(1, n_aps + 1):
                    ax = axs_flat[ap - 1]

                    x = (frame_data[f'x_{ap}'].values[0] - llx) / _bin
                    y = (frame_data[f'y_{ap}'].values[0] - lly) / _bin
                    flag = frame_data[f'flag_{ap}'].values[0]

                    ax.imshow(data, cmap='gray_r', vmin=vmin, vmax=vmax, origin='lower')
                    self._overlay_bad_pixels(ax, ccd_num, '1')

                    if pd.notna(x) and pd.notna(y):
                        color = 'lime' if flag == 0 else 'red'

                        ax.add_patch(patches.Circle((x, y), r_targ, edgecolor=color,
                                                    facecolor='none', lw=1.5, alpha=0.9))
                        ax.add_patch(patches.Circle((x, y), r_sky1, edgecolor='orange',
                                                    facecolor='none', lw=1.0, linestyle='--', alpha=0.8))
                        ax.add_patch(patches.Circle((x, y), r_sky2, edgecolor='red',
                                                    facecolor='none', lw=1.0, linestyle='--', alpha=0.8))

                        ax.set_xlim(x - zoom_box, x + zoom_box)
                        ax.set_ylim(y - zoom_box, y + zoom_box)
                        ax.set_title(f"Star {ap} (Flag: {flag})", fontsize=10, color=color)

                        divider = make_axes_locatable(ax)
                        ax_histx = divider.append_axes("top", size="25%", pad=0.05, sharex=ax)
                        ax_histy = divider.append_axes("right", size="25%", pad=0.05, sharey=ax)

                        x_start = max(0, int(x - zoom_box))
                        x_end = min(data.shape[1], int(x + zoom_box))
                        y_start = max(0, int(y - zoom_box))
                        y_end = min(data.shape[0], int(y + zoom_box))
                        cutout = data[y_start:y_end, x_start:x_end]

                        if cutout.size > 0:
                            local_median = np.nanmedian(cutout)
                            cutout_subbed = cutout - local_median
                            x_range = np.arange(x_start, x_end)
                            y_range = np.arange(y_start, y_end)
                            profile_x = np.nanmean(cutout_subbed, axis=0)
                            profile_y = np.nanmean(cutout_subbed, axis=1)

                            ax_histx.plot(x_range, profile_x, color='blue', lw=1.2, drawstyle='steps-mid')
                            ax_histy.plot(profile_y, y_range, color='blue', lw=1.2, drawstyle='steps-mid')
                            ax_histx.fill_between(x_range, profile_x, 0, step='mid', alpha=0.3, color='blue')
                            ax_histy.fill_betweenx(y_range, 0, profile_y, step='mid', alpha=0.3, color='blue')
                            ax_histx.axhline(0, color='black', lw=0.5, linestyle='--')
                            ax_histy.axvline(0, color='black', lw=0.5, linestyle='--')
                            ax_histx.set_yscale('symlog', linthresh=100.0)
                            ax_histy.set_xscale('symlog', linthresh=100.0)

                        ax_histx.axis('off')
                        ax_histy.axis('off')

                    else:
                        color = 'red'
                        ax.set_title(f"Star {ap} (Not Found)", fontsize=10, color='red')

                    ax.set_xticks([])
                    ax.set_yticks([])
                    for spine in ax.spines.values():
                        spine.set_color(color if pd.notna(x) else 'red')
                        spine.set_linewidth(2)

                for j in range(n_aps, len(axs_flat)):
                    axs_flat[j].axis('off')

                plt.subplots_adjust(top=0.92, bottom=0.05, left=0.05, right=0.95)
                save_file = os.path.join(zoom_dir, f'zoom_F{i+1}_{frame_name}.png')
                plt.savefig(save_file, dpi=150, bbox_inches='tight')
                plt.show()
                plt.close(fig)

            except Exception as e:
                print(f"{RED}Error plotting frame {i}: {e}{RESET}")

    # ------------------------------------------------------------------
    def plot_with_log(self, log_filename, ccd_num='1'):
        df_log, n_aps = read_hipercam_log(log_filename)
        if df_log is None:
            print(f"{RED}Error: Cannot read log file.{RESET}"); return

        csv_path = os.path.join(self.base_dir, 'all.csv')
        df_log.to_csv(csv_path)
        print(f"[Plot Check Tracking] Save {csv_path}")

        lis_file = self.lis[0] if self.lis else None
        if not lis_file or not os.path.exists(lis_file):
            print(f"{RED}Error: List file not found!{RESET}"); return

        with open(lis_file, 'r') as f:
            hcm_files = [line.strip() for line in f.readlines() if line.strip()]

        n_plot = len(hcm_files)
        print(f"{BLUE}--- Plotting {n_plot} frames ---{RESET}")

        ncols = 5
        nrows = math.ceil(n_plot / ncols)
        fig, axs = plt.subplots(nrows, ncols, figsize=(5*ncols, 5*nrows), facecolor='whitesmoke')
        axs_flat = np.atleast_1d(axs).flatten()
        zscale = ZScaleInterval()

        for i in range(n_plot):
            file_path = hcm_files[i]
            ax = axs_flat[i]

            try:
                mccd = hcam.MCCD.read(file_path)
                win_obj = mccd[ccd_num]['1']
                data = win_obj.data
                llx = getattr(win_obj, 'llx', 1)
                lly = getattr(win_obj, 'lly', 1)

                vmin, vmax = zscale.get_limits(data)
                ax.imshow(data, cmap='gray_r', vmin=vmin, vmax=vmax, origin='lower')
                self._overlay_bad_pixels(ax, ccd_num, '1')

                if i < len(df_log):
                    frame_data = df_log.iloc[[i]]
                    _bin = getattr(self, 'binning', 1)
                    if not frame_data.empty:
                        for ap in range(1, n_aps + 1):
                            x = (frame_data[f'x_{ap}'].values[0] - llx) / _bin
                            y = (frame_data[f'y_{ap}'].values[0] - lly) / _bin
                            fwhm = frame_data['mfwhm'].values[0] / _bin
                            flag = frame_data[f'flag_{ap}'].values[0]

                            color = 'lime' if flag == 0 else 'red'
                            plot_radius = fwhm
                            if pd.notna(x) and pd.notna(y):
                                ax.add_patch(patches.Circle((x, y), plot_radius * 1.8,
                                                            edgecolor=color, facecolor='none', lw=0.5, alpha=0.8))
                                ax.text(x + 1.5 * plot_radius, y, str(ap),
                                        color=color, fontsize=5, ha='center')

                ax.set_title(f"F{i+1}: {os.path.basename(file_path)}", fontsize=8)
                ax.axis('off')

            except Exception as e:
                ax.axis('off')
                print(f"Error plotting frame {i}: {e}")

        for j in range(n_plot, len(axs_flat)):
            axs_flat[j].axis('off')

        plt.tight_layout()
        save_file = os.path.join(self.base_dir, f'reduction_check_{ccd_num}.png')
        plt.savefig(save_file, dpi=150, bbox_inches='tight')
        plt.show()
        plt.close(fig)
