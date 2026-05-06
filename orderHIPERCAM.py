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

    def setaper(self, ccd_label='1', win_label='1', SIGMA_THRESHOLD=1.5,
                output_plot="detection_labeled.png", SKIP_BRIGHTEST=10,
                MARGIN_LEFT=15, MARGIN_RIGHT=15, MARGIN_BOTTOM=5, MARGIN_TOP=27,
                R_TARG=None, R_SKY1=16, R_SKY2=24, frame=5, diagnostics=False):
        """
        Detect stars and write a HiPERCAM aperture (.ape) file.

        In single_run mode one representative frame (index=frame) is used.
        In multi-run mode every file in the list gets its own .ape file.

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
            print(red_file)
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
                return

            if diagnostics:
                
                full_log = os.path.join(self.base_dir, log_file)
                print(red_file, full_log)
                print(f"{BLUE}--- Diagnostics: running plot_with_log ---{RESET}")
                self.plot_with_log(full_log)
                print(f"{BLUE}--- Diagnostics: running plot_with_zoom ---{RESET}")
                self.plot_with_zoom(full_log)
        else:
            print("Not implemented yet")

    # ------------------------------------------------------------------
    def solvwcs(self, ccd_label='1', win_label='1', log_file=None,
                output_csv='wcs_radec.csv',
                ra_center=None, dec_center=None, radius=5.0,
                scale_low=None, scale_high=None,
                astrometry_cache='astrometry_cache',
                verbose=False):
        """
        Solve WCS for every frame using aperture x/y from the reduce .log.

        Strategy
        --------
        For each frame, the x/y positions of ALL apertures are passed as
        the star catalog to astrometry.net.

          - Frame 0  : blind solve (or position-hinted if ra_center/dec_center
                       are provided).
          - Frame N>0: previous frame's solved centre used as position hint,
                       which makes each subsequent solve fast.

        A single WCS is needed for the whole observation only if the field
        does not drift. The per-frame strategy handles small tracking errors
        automatically.

        Parameters
        ----------
        ccd_label : str
            CCD label matching the one used in setaper/reduce.
        win_label : str
            Window label (e.g. '1').
        log_file : str, optional
            Path to reduce .log file. Defaults to <base_dir>/all.log.
        output_csv : str
            Output CSV filename written inside base_dir.
        ra_center, dec_center : float, optional
            Position hint in degrees for the first-frame solve.
        radius : float
            Search radius in degrees for the position hint.
        scale_low, scale_high : float, optional
            Plate-scale bounds arcsec/binned-px. Derived automatically from
            ULTRASPEC_SCALE_ARCSEC_PX × xbin if not provided.
        astrometry_cache : str
            Directory where astrometry.net index files are cached.
        verbose : bool
            Print full internal solver progress (INFO-level logging).

        Returns
        -------
        pandas.DataFrame or None
            Columns: frame, mjd, aperture, x, y, ra_deg, dec_deg
        """
        import astrometry
        import pandas as pd

        # ── 1. Get binning from a reference HCM ──────────────────────
        hcm_files = []
        for lis_path in self.lis:
            if os.path.exists(lis_path):
                with open(lis_path) as fh:
                    hcm_files.extend(l.strip() for l in fh if l.strip())
        if not hcm_files:
            print(f"{RED}[solvwcs] No HCM files found in list files.{RESET}")
            return None

        ref_hcm = hcm_files[min(3, len(hcm_files) - 1)]
        try:
            mccd   = hcam.MCCD.read(ref_hcm)
            window = mccd[ccd_label][win_label]
            xbin   = int(getattr(window, 'xbin', 1))
            ybin   = int(getattr(window, 'ybin', 1))
        except Exception as e:
            print(f"{RED}[solvwcs] Cannot read reference HCM: {e}{RESET}")
            return None

        eff_scale = ULTRASPEC_SCALE_ARCSEC_PX * xbin
        if scale_low  is None:
            scale_low  = eff_scale * 0.95
        if scale_high is None:
            scale_high = eff_scale * 1.05
        print(f"[solvwcs] Binning={xbin}x{ybin}  "
              f"scale={eff_scale:.3f}\"/px  "
              f"hint range [{scale_low:.3f}, {scale_high:.3f}]")

        # ── 2. Read log ───────────────────────────────────────────────
        if log_file is None:
            log_file = os.path.join(self.base_dir, 'all.log')
        ccd_log = self._read_hlog(log_file, ccd_label)
        if ccd_log is None:
            return None

        ap_labels = sorted(ccd_log.keys(),
                           key=lambda k: int(k) if str(k).isdigit() else k)
        n_frames = len(next(iter(ccd_log.values()))['t'])
        n_aps    = len(ap_labels)
        print(f"[solvwcs] {n_aps} apertures × {n_frames} frames")

        if n_aps < 6:
            print(f"{RED}[solvwcs] Need ≥6 apertures for plate solving, "
                  f"got {n_aps}. Run setaper() with more stars.{RESET}")
            return None

        # ── 3. Load index files as Path objects ───────────────────────
        # astrometry.Solver calls path.resolve() internally — must be Path
        def _idx(series, scales, label):
            try:
                files   = series.index_files(
                    cache_directory=astrometry_cache, scales=scales)
                on_disk = [Path(f) for f in files if Path(f).exists()]
                print(f"[solvwcs]   {label}: {len(on_disk)} on disk")
                return on_disk
            except Exception as e:
                print(f"[solvwcs]   {label}: {e}")
                return []

        print("[solvwcs] Loading index files...")
        index_files = list(dict.fromkeys(
            _idx(astrometry.series_5200, {2, 3, 4}, "series_5200 {2,3,4}")
            + _idx(astrometry.series_4200, {2, 3, 4}, "series_4200 {2,3,4}")
        ))
        if not index_files:
            print(f"{RED}[solvwcs] No index files found in "
                  f"{astrometry_cache}{RESET}")
            return None

        size_hint = astrometry.SizeHint(lower_arcsec_per_pixel=scale_low,
                                         upper_arcsec_per_pixel=scale_high)
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

        # ── 4. Per-frame solve ────────────────────────────────────────
        records  = []
        prev_ra  = ra_center
        prev_dec = dec_center
        prev_wcs = None

        print(f"{BLUE}[solvwcs] Solving {n_frames} frames...{RESET}")

        with astrometry.Solver(index_files) as solver:
            for fi in range(n_frames):

                # Star list for this frame.
                # Log x/y are in unbinned window-pixel coords (HiPERCAM stores
                # positions × xbin in the ape file and the log follows the same
                # convention).  Divide by xbin/ybin to get binned pixel coords
                # so they match the SizeHint scale (arcsec / binned-px).
                stars = []
                for ap in ap_labels:
                    x = float(ccd_log[ap]['x'][fi]) / xbin
                    y = float(ccd_log[ap]['y'][fi]) / ybin
                    if np.isfinite(x) and np.isfinite(y):
                        stars.append((x, y))

                # Position hint: None for blind frame-0, previous centre after
                position_hint = None
                if prev_ra is not None and prev_dec is not None:
                    position_hint = astrometry.PositionHint(
                        ra_deg=prev_ra, dec_deg=prev_dec, radius_deg=radius)

                wcs = None
                if len(stars) >= 6:
                    solution = solver.solve(
                        stars=stars,
                        size_hint=size_hint,
                        position_hint=position_hint,
                        solution_parameters=params,
                    )
                    if solution.has_match():
                        match    = solution.best_match()
                        wcs      = match.astropy_wcs()   # 0-indexed px convention
                        prev_ra  = match.center_ra_deg
                        prev_dec = match.center_dec_deg
                        prev_wcs = wcs
                        if fi == 0:
                            print(f"{BLUE}[solvwcs] Frame 0 solved: "
                                  f"RA={prev_ra:.5f}°  "
                                  f"Dec={prev_dec:.5f}°{RESET}")
                    elif fi == 0:
                        hint_info = (f"RA={prev_ra:.4f} Dec={prev_dec:.4f}"
                                     if prev_ra is not None else "no hint")
                        print(f"{RED}[solvwcs] Frame 0 solve failed "
                              f"({hint_info}).{RESET}")
                        print(f"{RED}           scale hint [{scale_low:.3f}, "
                              f"{scale_high:.3f}] arcsec/px — verify this "
                              f"matches the true plate scale.{RESET}")
                        if verbose:
                            logging.getLogger().setLevel(logging.WARNING)
                        return None

                # Fall back to previous WCS if this frame failed
                if wcs is None:
                    wcs = prev_wcs

                # Apply WCS to every aperture in this frame.
                # Divide by xbin/ybin for the same reason as above.
                for ap in ap_labels:
                    x_raw = float(ccd_log[ap]['x'][fi])
                    y_raw = float(ccd_log[ap]['y'][fi])
                    mjd   = float(ccd_log[ap]['t'][fi])
                    x_bin = x_raw / xbin
                    y_bin = y_raw / ybin
                    ra = dec = np.nan
                    if wcs is not None and np.isfinite(x_bin) and np.isfinite(y_bin):
                        try:
                            sky = wcs.pixel_to_world(x_bin, y_bin)
                            ra  = float(sky.ra.deg)
                            dec = float(sky.dec.deg)
                        except Exception:
                            pass
                    records.append({
                        'frame':    fi + 1,
                        'mjd':      mjd,
                        'aperture': int(ap) if str(ap).isdigit() else ap,
                        'x':        x_raw,   # keep original log coords in CSV
                        'y':        y_raw,
                        'ra_deg':   ra,
                        'dec_deg':  dec,
                    })

                if (fi + 1) % 100 == 0 or fi == n_frames - 1:
                    print(f"[solvwcs]   {fi + 1}/{n_frames} frames done")

        if verbose:
            logging.getLogger().setLevel(logging.WARNING)

        if not records:
            print(f"{RED}[solvwcs] No records produced.{RESET}")
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
        Parse a HiPERCAM reduce log.

        Returns
        -------
        dict {ap_label: {'t': ndarray, 'x': ndarray, 'y': ndarray}} or None.
        """
        # Attempt 1: hipercam.hlog API
        try:
            import importlib
            hlog_mod = importlib.import_module('hipercam.hlog')
            # Constructor signature varies by hipercam version
            for loader in [lambda: hlog_mod.Hlog(log_file),
                           lambda: hlog_mod.Hlog.read(log_file)]:
                try:
                    log_data = loader()
                    break
                except Exception:
                    continue
            else:
                raise RuntimeError("no valid Hlog loader")
            ccd_data = log_data[str(ccd_label)]
            result   = {
                ap: {'t': np.asarray(ts.t),
                     'x': np.asarray(ts.x),
                     'y': np.asarray(ts.y)}
                for ap, ts in ccd_data.items()
            }
            print(f"[solvwcs] Log parsed via hipercam.hlog "
                  f"({len(result)} apertures)")
            return result
        except Exception as e:
            print(f"[solvwcs] hipercam.hlog failed ({e}), trying pandas...")

        # Attempt 2: pandas ASCII fallback
        # Log header format:  # N = CCD nframe MJD MJDok Exptim ... x_1 xe_1 y_1 ...
        # Column names start from 'nframe' — skip the leading "N = CCD" tokens.
        try:
            import pandas as pd
            header_cols = None
            with open(log_file) as fh:
                for line in fh:
                    if not line.startswith('#'):
                        break
                    tokens = line.lstrip('#').strip().split()
                    # True header line has both 'nframe' and 'MJD' among many tokens
                    if 'nframe' in tokens and 'MJD' in tokens:
                        nframe_idx = tokens.index('nframe')
                        header_cols = tokens[nframe_idx:]   # strip "N = CCD" prefix
                        break

            if header_cols is None:
                print(f"{RED}[solvwcs] Cannot find column header in "
                      f"{log_file}{RESET}")
                return None

            df = pd.read_csv(log_file, comment='#', sep=r'\s+',
                             names=header_cols, engine='python')
            df.columns = [c.lower() for c in df.columns]

            mjd_col = next((c for c in df.columns
                            if c in ('mjd', 'mjd_mid', 'tmid', 't')),
                           df.columns[0])

            # Columns are x_1, y_1, x_2, y_2, ...
            ap_ids = sorted(
                {c.split('_')[1] for c in df.columns
                 if c.startswith('x_')
                 and len(c.split('_')) == 2
                 and c.split('_')[1].isdigit()},
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

            if result:
                print(f"[solvwcs] Log parsed via pandas fallback "
                      f"({len(result)} apertures)")
                return result

            print(f"{RED}[solvwcs] No x_N/y_N columns found in log.{RESET}")
            return None

        except Exception as exc:
            print(f"{RED}[solvwcs] Log parse failed: {exc}{RESET}")
            return None

    # ------------------------------------------------------------------
    def plot_with_zoom(self, log_filename, ccd_num='1', zoom_box=15):
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
                data = mccd[ccd_num]['1'].data
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

                fig.suptitle(f"Frame {i+1}: {frame_name}", fontsize=16, fontweight='bold')
        
                for ap in range(1, n_aps + 1):
                    ax = axs_flat[ap - 1]
                    
                    x = frame_data[f'x_{ap}'].values[0] / self.binning
                    y = frame_data[f'y_{ap}'].values[0] / self.binning
                    fwhm = frame_data[f'fwhm_{ap}'].values[0] / self.binning
                    flag = frame_data[f'flag_{ap}'].values[0]

                    ax.imshow(data, cmap='gray_r', vmin=vmin, vmax=vmax, origin='lower')
                    
                    if pd.notna(x) and pd.notna(y):
                        color = 'lime' if flag == 0 else 'red'
                        plot_radius = fwhm if (pd.notna(fwhm) and fwhm > 0) else 5.0
                        
                        circle = patches.Circle((x, y), plot_radius * 1.8, edgecolor=color, 
                                                facecolor='none', lw=1.5, alpha=0.8)
                        ax.add_patch(circle)
                        
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
                        # ==========================================
                        # ==========================================

                    else:
                        ax.set_title(f"Star {ap} (Not Found)", fontsize=10, color='red')
                
                    ax.set_xticks([])
                    ax.set_yticks([])
                    for spine in ax.spines.values():
                        spine.set_color(color if pd.notna(x) else 'red')
                        spine.set_linewidth(2)


                for j in range(n_aps, len(axs_flat)):
                    axs_flat[j].axis('off')
                
                # Use subplots_adjust to manually set margins instead of tight_layout
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
        df_log.to_csv(os.path.join(self.bias_dir, 'all.csv'))
        print(f'[Plot Check Tracking] Save {os.path.join(self.bias_dir, 'all.csv')}')
        
        if df_log is None:
            print(f"{RED}Error: Cannot read log file.{RESET}"); return

        lis_file = self.lis[0] if self.lis else None
        if not lis_file or not os.path.exists(lis_file):
            print(f"{RED}Error: List file not found!{RESET}"); return

        with open(lis_file, 'r') as f:
            hcm_files = [line.strip() for line in f.readlines() if line.strip()]

        ntot = len(hcm_files)
        n_plot = ntot
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
                data = mccd[ccd_num]['1'].data
                vmin, vmax = zscale.get_limits(data)
                ax.imshow(data, cmap='gray_r', vmin=vmin, vmax=vmax, origin='lower')
                
                # Use iloc to get data for the current frame
                if i < len(df_log):
                    frame_data = df_log.iloc[[i]]

                    if not frame_data.empty:
                        for ap in range(1, n_aps + 1):
                            x = frame_data[f'x_{ap}'].values[0]/self.binning
                            y = frame_data[f'y_{ap}'].values[0]/self.binning
                            # fwhm = frame_data[f'fwhm_{ap}'].values[0]/self.binning
                            fwhm = frame_data[f'mfwhm'].values[0]/self.binning
                            flag = frame_data[f'flag_{ap}'].values[0]
                       
                            color = 'lime' if flag == 0 else 'red'
                            plot_radius = fwhm  #if (pd.notna(fwhm) and fwhm > 0) else 5.0
                            
                            # info_text = f"Ap{ap}\n({x}, {y})\nF:{fwhm}"
                            info_text = f"{ap}"
                            if pd.notna(x) and pd.notna(y):
                                circle = patches.Circle((x, y), plot_radius*1.8, edgecolor=color, 
                                                       facecolor='none', lw=.5, alpha=0.8)
                                ax.add_patch(circle)
                                ax.text(x +1.5*plot_radius, y , info_text, color=color, fontsize=5, ha='center')

                ax.set_title(f"F{i+1}: {os.path.basename(file_path)}", fontsize=8)
                ax.axis('off')

            except Exception as e:
                ax.axis('off')
                print(f"Error plotting frame {i}: {e}")

        for j in range(n_plot, len(axs_flat)): axs_flat[j].axis('off')
        
        plt.tight_layout()
        save_file = os.path.join(self.base_dir, f'reduction_check_{ccd_num}.png')
        plt.savefig(save_file, dpi=150, bbox_inches='tight')
        plt.show()
        plt.close(fig)