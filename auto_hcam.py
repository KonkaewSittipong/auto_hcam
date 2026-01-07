import os
import sys
import math
import shutil
import glob
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from datetime import datetime
from typing import Optional

# Astropy / Photutils
from astropy.visualization import ZScaleInterval
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder

# HiPERCAM specific imports
import hipercam as hcam
from hipercam import MCCD 
import hipercam.scripts as scripts
import hipercam.spooler as spooler
import hipercam.mpl as hmpl
import hipercam.ucam as ucam

# =============================================================================
# 1. HELPER FUNCTIONS (Plotting & File Reading)
# =============================================================================

def plot_run(run_name, run_type="Data", ccd_to_plot='1', save_dir=None, frames_range=None):
    """Plots a ZScale contact sheet (5 columns) AND a median level trend with range markers."""
    print(f"\n... Processing {run_type} Run: {run_name} (CCD {ccd_to_plot}) ...")
    zscale = ZScaleInterval()
    frame_numbers, median_levels = [], []
    
    try:
        rdata = ucam.Rdata(run_name, nframe=1, server=False)
        ntot = rdata.ntotal()
    except:
        ntot = 10

    # Part A: Grid Plot
    ncols = 5
    nrows = math.ceil(ntot / ncols)
    fig_grid, axs = plt.subplots(nrows, ncols, figsize=(6*ncols, 6*nrows), 
                                 gridspec_kw={'hspace': 0.15, 'wspace': 0.05})
    axs_flat = axs.flatten() if ntot > 1 else [axs]
    
    try:
        with spooler.data_source('ul', run_name) as spool:
            for n, mccd in enumerate(spool):
                if n >= len(axs_flat): break
                ax = axs_flat[n]
                frame_num = mccd.head.get('NFRAME', n+1)
                
                if ccd_to_plot in mccd:
                    ax.cla()
                    pixel_data = [wind.data.flatten() for wind in mccd[ccd_to_plot].values()]
                    all_pixels = np.concatenate(pixel_data)
                    
                    median_val = np.nanmedian(all_pixels)
                    frame_numbers.append(frame_num)
                    median_levels.append(median_val)
                    
                    vmin, vmax = zscale.get_limits(all_pixels)
                    hmpl.pCcd(ax, mccd[ccd_to_plot], iset='d', dlo=vmin, dhi=vmax)
                    ax.set_title(f"Fr {frame_num}\nMed: {median_val:.1f}", fontsize=8)
                    ax.set_aspect('equal')
                    ax.set_xticks([]); ax.set_yticks([])

            for i in range(n + 1, len(axs_flat)):
                axs_flat[i].axis('off')

            if save_dir:
                grid_path = os.path.join(save_dir, f"{run_name}_grid.png")
                fig_grid.savefig(grid_path, bbox_inches='tight')
                print(f"Saved: {grid_path}")
            
            plt.show() 

        # Part B: Median Trend (Modified with Vertical Lines)
        if len(frame_numbers) > 0:
            fig_trend = plt.figure(figsize=(15, 5))
            plt.plot(frame_numbers, median_levels, 'o-', markersize=4, color='tab:blue', label='Median Level')
            
            # --- Draw Vertical Lines for Frame Selection ---
            if frames_range:
                start_f = frames_range[0]
                end_f = frames_range[1]
                plt.axvline(start_f, color='red', linestyle='--', linewidth=2, label=f'Start: {start_f}')
                
                if end_f > 0:
                    plt.axvline(end_f, color='red', linestyle='--', linewidth=2, label=f'End: {end_f}')
                else:
                    plt.axvline(frame_numbers[-1], color='orange', linestyle=':', linewidth=2, label='End of Data')

                plt.legend()
            # ---------------------------------------------------

            plt.title(f"Median Levels: {run_name} ({run_type})")
            plt.xlabel("Frame Number"); plt.ylabel("Counts (ADU)")
            plt.grid(True, alpha=0.3)
            
            if save_dir:
                trend_path = os.path.join(save_dir, f"{run_name}_median_trend.png")
                plt.savefig(trend_path, bbox_inches='tight')
                print(f"Saved: {trend_path}")
            
            plt.show() 

    except Exception as e:
        print(f"  Plot Error: {e}")
        plt.close('all')

def read_hipercam_log(filename):
    """Reads HiPERCAM .log file into a DataFrame (Standard Read)."""
    print(f"\n--- Reading Data Log: {filename} ---")
    if not os.path.exists(filename):
        print(f"❌ Error: {filename} not found.")
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
    
    print(f"Detected {max_aperture} apertures.")

    all_cols = base_cols.copy()
    for i in range(1, max_aperture + 1):
        for col in ap_cols_template:
            all_cols.append(f"{col}_{i}")

    try:
        df = pd.read_csv(filename, comment='#', sep=r'\s+', header=None, names=all_cols)
        df['nframe'] = pd.to_numeric(df['nframe'], errors='coerce').fillna(-1).astype(int)
        return df, max_aperture
    except Exception as e:
        print(f"Error reading log: {e}")
        return None, 0

def read_hipercam_log_to_dataframe_ensemble(filename):
    """
    Reads a HiPERCAM .log file specially for the Ensemble process.
    Returns: clean_df, start_col, raw_df
    """
    print(f"--- Reading Log for Ensemble: {filename} ---")
    
    max_ap = 0
    with open(filename, 'r') as f:
        for line in f:
            if 'flag_' in line:
                parts = line.split()
                for p in parts:
                    if 'flag_' in p:
                        try:
                            ap_num = int(p.replace('flag_', ''))
                            max_ap = max(max_ap, ap_num)
                        except: pass
            if not line.startswith('#'): break

    print(f"Detected {max_ap} apertures.")

    try:
        raw_df = pd.read_csv(filename, comment='#', sep=r'\s+', header=None)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return None, None, None

    clean_df = pd.DataFrame()
    clean_df['MJD'] = raw_df.iloc[:, 2] 
    clean_df['Exptim'] = raw_df.iloc[:, 4]
    
    global_cols = 7
    block_size = 16 
    
    for i in range(max_ap):
        star_id = i + 1
        start_idx = global_cols + (i * block_size)
        
        counts_idx = start_idx + 8
        countse_idx = start_idx + 9
        
        if countse_idx < len(raw_df.columns):
            clean_df[f'flux_{star_id}'] = raw_df.iloc[:, counts_idx]
            clean_df[f'err_{star_id}'] = raw_df.iloc[:, countse_idx]
        else:
            print(f"Warning: Columns missing for star {star_id}")

    return clean_df, 2, raw_df

def plot_grid_white_background(df, aps_count, list_file="all.lis", save_path="."):
    """Plots final results with circles around detected stars."""
    if not os.path.exists(list_file):
        print(f"Error: {list_file} not found.")
        return

    with open(list_file) as f:
        files = [l.strip() for l in f if l.strip()]
    
    print(f"Plotting {len(files)} files on WHITE background...")

    ncols = 5
    nrows = math.ceil(len(files) / ncols)
    
    fig_grid, axs = plt.subplots(nrows, ncols, figsize=(6*ncols, 6*nrows), 
                                 gridspec_kw={'hspace': 0.15, 'wspace': 0.05}, facecolor='white')
    
    if len(files) > 1:
        axs_flat = axs.flatten()
    else:
        axs_flat = [axs]
    
    zscale = ZScaleInterval()
    target_ccd = 1 

    for i, fname in enumerate(files):
        ax = axs_flat[i]
        try:
            mccd = MCCD.read(fname)
            data = mccd[str(target_ccd)]['1'].data 
            
            vmin, vmax = zscale.get_limits(data)
            ax.imshow(data, origin='lower', cmap='gray_r', vmin=vmin, vmax=vmax)
            
            if df is not None:
                if i < len(df):
                    row = df.iloc[[i]]
                    if not row.empty:
                        for ap in range(1, aps_count + 1):
                            col_x, col_y, col_fwhm = f'x_{ap}', f'y_{ap}', f'fwhm_{ap}'
                            
                            if col_x in row.columns and col_y in row.columns:
                                x_ = row.iloc[0][col_x]
                                y_ = row.iloc[0][col_y]
                                fwhm_ = row.iloc[0][col_fwhm]
                                
                                h, w = data.shape
                                if np.isfinite(x_) and np.isfinite(y_):
                                    x_plot, y_plot = x_ / 2, y_ / 2

                                    if 0 <= x_plot <= w and 0 <= y_plot <= h:
                                        base_radius = fwhm_ * 2 if (np.isfinite(fwhm_) and fwhm_ > 0) else 10
                                        
                                        targ = patches.Circle((x_plot-0.5, y_plot-0.5), base_radius/2, 
                                                              edgecolor='lime', facecolor='none', lw=1.0, alpha=0.8)
                                        ax.add_patch(targ)
                                        ax.text(x_plot, y_plot + base_radius + 5, f"{ap}", 
                                                color='red', fontsize=5, ha='center', va='bottom')
                else:
                    ax.text(0.5, 0.9, "No Log Data", color='orange', ha='center', transform=ax.transAxes, fontsize=8, weight='bold')

            ax.set_title(f"#{i+1}: {os.path.basename(fname)}", fontsize=7, color='black')
            ax.axis('off')
            
        except Exception as e:
            print(f"Error plotting {fname}: {e}")
            ax.axis('off')
            
    for j in range(len(files), len(axs_flat)):
        axs_flat[j].axis('off')

    plt.tight_layout()
    output_png = os.path.join(save_path, "grid_all_stars_simple.png")
    plt.savefig(output_png, bbox_inches='tight', dpi=150)
    print(f"✅ Saved to {output_png}")
    plt.show() 

# =============================================================================
# 2. THE UNIFIED PIPELINE CLASS
# =============================================================================

class CompleteUltraspecPipeline:
    def __init__(self, config, dir_save='g'):
        self.config = config
        self.dir_save = os.path.abspath(dir_save)
        self.dir_figs = os.path.join(self.dir_save, "figs")
        self.dir_results = os.path.join(self.dir_save, "results")
        
        self.calibs = {'bias': 'none', 'dark': 'none', 'flat': 'none'}
        
        # Ensemble Class Attributes
        self.data: Optional[pd.DataFrame] = None
        self.raw_original: Optional[pd.DataFrame] = None
        self.first_column: Optional[int] = None
        self.number_stars: Optional[int] = None
        
        # Create directories
        if not os.path.exists(self.dir_save):
            os.makedirs(self.dir_save)
            print(f"Created directory: {self.dir_save}")
        
        if not os.path.exists(self.dir_figs):
            os.makedirs(self.dir_figs)
            print(f"Created directory: {self.dir_figs}")
            
        if not os.path.exists(self.dir_results):
            os.makedirs(self.dir_results)

    def get_frame_args(self, section):
        frames = self.config[section].get('frames', [1, 0])
        return str(frames[0]), str(frames[1])

    def move_generated_hcm(self, run_name):
        pattern = f"{run_name}_*.hcm"
        if run_name in ["bias", "flat"]: pattern = f"{run_name}.hcm"
        found_files = glob.glob(pattern)
        for f in found_files:
            dest = os.path.join(self.dir_save, os.path.basename(f))
            if os.path.exists(dest): os.remove(dest)
            shutil.move(f, self.dir_save)

    # --- STAGE 1-2: CALIBRATION ---
    def run_calibration(self, plot_calib = True):
        print("\n" + "="*40 + "\nSTAGE 1: Master Calibration\n" + "="*40)
        # Bias
        if self.config['bias']['runs']:
            run = self.config['bias']['runs'][0]
            f1, f2 = self.get_frame_args('bias')
            args = ["nodef", "source=ul", f"run={run}", f"first={f1}", f"last={f2}",
                    "sigma=3.0", "plot=no", "output=bias.hcm", "clobber=True"]
            scripts.makebias(args)
            self.move_generated_hcm("bias")
            self.calibs['bias'] = os.path.join(self.dir_save, "bias.hcm")
        # Flat
        if self.config['flat']['runs']:
            run = self.config['flat']['runs'][0]
            f1, f2 = self.get_frame_args('flat')
            args = ["nodef", "source=ul", f"run={run}", f"first={f1}", f"last={f2}", 
                    f"bias={self.calibs['bias']}", "dark=none", "ngroup=5", "ccd=0",
                    "lower=10000", "upper=35000", "plot=no", "output=flat.hcm", "clobber=True"]
            scripts.makeflat(args)
            self.move_generated_hcm("flat")
            self.calibs['flat'] = os.path.join(self.dir_save, "flat.hcm")
        
        if plot_calib == True : self.plot_calibrations()
        
    def plot_calibrations(self, ccd_to_plot='1'):
        print("\n" + "="*40 + "\nSTAGE 1b: Plotting Calibrations\n" + "="*40)
        
        # --- Extract frame info for lines ---
        bias_frames = self.config['bias'].get('frames')
        flat_frames = self.config['flat'].get('frames')

        if self.config['bias']['runs']:
            plot_run(self.config['bias']['runs'][0], "Bias Run", ccd_to_plot, self.dir_figs, frames_range=bias_frames)
        if self.config['flat']['runs']:
            plot_run(self.config['flat']['runs'][0], "Flat Run", ccd_to_plot, self.dir_figs, frames_range=flat_frames)

        zscale = ZScaleInterval()
        masters = [f for f in [self.calibs['bias'], self.calibs['flat']] if os.path.exists(f)]
        
        if masters:
            fig, axs = plt.subplots(1, len(masters), figsize=(14, 7))
            if len(masters) == 1: axs = [axs]
            
            for ax, fpath in zip(axs, masters):
                mccd = MCCD.read(fpath)
                pixel_data = [w.data.flatten() for w in mccd[ccd_to_plot].values()]
                all_pixels = np.concatenate(pixel_data)
                
                vmin, vmax = zscale.get_limits(all_pixels)
                hmpl.pCcd(ax, mccd[ccd_to_plot], iset='d', dlo=vmin, dhi=vmax)
                
                hplot_text = (f"CCD {ccd_to_plot}\nplot range = \n{vmin:.15f} to \n{vmax:.15f}")
                ax.text(0.05, 0.95, hplot_text, transform=ax.transAxes, fontsize=8.5, 
                        family='monospace', verticalalignment='top', color='cyan', 
                        fontweight='bold', bbox=dict(boxstyle='square', facecolor='black', alpha=0.8))
                
                ax.set_title(f"MASTER: {os.path.basename(fpath)}")
                ax.set_xticks([]); ax.set_yticks([])
            
            save_path = os.path.join(self.dir_figs, "final_masters_check.png")
            plt.savefig(save_path, bbox_inches='tight')
            print(f"Saved: {save_path}")
            plt.show() 

    def clean_and_plot_science(self):
        print("\n" + "="*40 + "\nSTAGE 2: Calibrated Extraction & Plotting\n" + "="*40)
        f1, f2 = self.get_frame_args('data')
        
        # --- Extract frame info for lines ---
        data_frames = self.config['data'].get('frames')

        for run in self.config['data']['runs']:
            args = ["nodef", f"run={run}", "source=ul", f"first={f1}", f"last={f2}",
                    "ndigit=3", "trim=False", f"bias={self.calibs['bias']}", "dark=none",
                    f"flat={self.calibs['flat']}", "fmap=none", f"output={run}"]
            try:
                scripts.grab(args)
                self.move_generated_hcm(run)
                plot_run(run, "Science Check", save_dir=self.dir_figs, frames_range=data_frames)
            except Exception as e:
                print(f"Error processing {run}: {e}")

    def plot_all_hcm(self, ccd_to_plot='1'):
        print("\n" + "="*40 + "\nSTAGE 3: Final Cleaned Grid\n" + "="*40)
        hcm_files = sorted(glob.glob(os.path.join(self.dir_save, "run*_*.hcm")))
        if not hcm_files: return
        
        zscale = ZScaleInterval()
        ncols = 5
        nrows = math.ceil(len(hcm_files) / ncols)
        fig, axs = plt.subplots(nrows, ncols, figsize=(6*ncols, 6*nrows), gridspec_kw={'hspace':0.2, 'wspace': 0.05})
        axs_flat = [axs] if len(hcm_files) == 1 else axs.flatten()

        for n, fpath in enumerate(hcm_files):
            ax = axs_flat[n]
            try:
                mccd = MCCD.read(fpath)
                pixel_data = [w.data.flatten() for w in mccd[ccd_to_plot].values()]
                all_pixels = np.concatenate(pixel_data)
                vmin, vmax = zscale.get_limits(all_pixels)
                hmpl.pCcd(ax, mccd[ccd_to_plot], iset='d', dlo=vmin, dhi=vmax)
                ax.set_title(os.path.basename(fpath), fontsize=10)
                ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
            except: ax.axis('off')

        for i in range(len(hcm_files), len(axs_flat)): axs_flat[i].axis('off')
        
        save_path = os.path.join(self.dir_figs, "ALL_CLEANED_GRID.png")
        plt.savefig(save_path, bbox_inches='tight')
        print(f"Saved: {save_path}")
        plt.show() 

    # --- STAGE 4: APERTURE SETUP ---
    def setaper(self, ccd_label='1', win_label='1', SIGMA_THRESHOLD=5,
                output_plot="detection_labeled.png", output_ape_name="ape.ape",
                MARGIN_LEFT=15, MARGIN_RIGHT=15, MARGIN_BOTTOM=5, MARGIN_TOP=25,
                R_TARG=10, R_SKY1=15, R_SKY2=20):
        
        base_path = self.dir_save
        hcm_files = glob.glob(os.path.join(base_path, "run*.hcm"))
        hcm_files.sort()
        if hcm_files:
            target_file = hcm_files[0]
            print(f"Analyzing file: {target_file}")
            
            try:
                mccd = hcam.MCCD.read(target_file)
                window_obj = mccd[ccd_label][win_label]
                data = window_obj.data
                try:
                    bining = window_obj.xbin
                    print(f"Detected Binning: {bining}")
                except AttributeError:
                    bining = 1
                height, width = data.shape
            except Exception as e:
                print(f"CRITICAL Error reading file: {e}")
                exit()
        
            print(f"Detecting sources > {SIGMA_THRESHOLD} sigma...")
            mean, median, std = sigma_clipped_stats(data, sigma=3)
            data_sub = data - median 
            
            daofind = DAOStarFinder(fwhm=4.0, threshold=SIGMA_THRESHOLD * std)
            sources = daofind(data_sub)
            if sources:
                mask = ((sources['xcentroid'] > MARGIN_LEFT) & (sources['xcentroid'] < width - MARGIN_RIGHT) & 
                        (sources['ycentroid'] > MARGIN_BOTTOM) & (sources['ycentroid'] < height - MARGIN_TOP))
                sources = sources[mask]
                sources.sort('flux')
                sources.reverse()

                print(f"Original source count: {len(sources)}")
                if len(sources) > 5:
                    print("Ignoring first 3 brightest stars...")
                    sources = sources[3:]
                else:
                    print("Warning: Less than 3 stars found. Keeping all.")

                print(f"Final usable stars: {len(sources)}")
                
                output_ape = os.path.join(base_path, output_ape_name)
                ccd_aps = []
                for i, source in enumerate(sources):
                    x, y = source['xcentroid'], source['ycentroid']
                    ap_id = str(i + 1) 
                    is_ref = (i == 0)
                    link_to = "" if is_ref else "1"
                    
                    ccd_aps.append([ap_id, {
                        "Comment": "hipercam.Aperture", 
                        "x": float(x) * bining, "y": float(y) * bining,
                        "rtarg": R_TARG * bining, "rsky1": R_SKY1 * bining, "rsky2": R_SKY2 * bining,
                        "ref": is_ref, "compo": False, "mask": [], "extra": [], "link": link_to 
                    }])

                inner_ccd_structure = ["hipercam.CcdAper"] + ccd_aps
                ape_json = ["hipercam.MccdAper", [ccd_label, inner_ccd_structure]]
                
                with open(output_ape, 'w') as f: json.dump(ape_json, f, indent=2)
                print(f"✅ Created: {output_ape}")
        
                # PLOT
                print("Generating plot...")
                fig, ax = plt.subplots(figsize=(10, 10))
                zscale = ZScaleInterval()
                vmin, vmax = zscale.get_limits(data)
                ax.imshow(data, origin='lower', cmap='Greys', vmin=vmin, vmax=vmax)
                
                rect = patches.Rectangle((MARGIN_LEFT, MARGIN_BOTTOM), width - MARGIN_LEFT - MARGIN_RIGHT, 
                                         height - MARGIN_BOTTOM - MARGIN_TOP,
                                         linewidth=2, edgecolor='yellow', facecolor='none', linestyle='--')
                ax.add_patch(rect)
        
                for i, source in enumerate(sources):
                    x, y = source['xcentroid'], source['ycentroid']
                    ap_number = i + 1
                    ax.add_patch(patches.Circle((x, y), R_TARG, edgecolor='lime', facecolor='none', lw=1.5))
                    ax.add_patch(patches.Circle((x, y), R_SKY1, edgecolor='red', facecolor='none', lw=1, linestyle='--'))
                    ax.add_patch(patches.Circle((x, y), R_SKY2, edgecolor='red', facecolor='none', lw=1, linestyle='--'))
                    ax.text(x, y + R_SKY2 + 2, str(ap_number), color='cyan', fontsize=12, fontweight='bold', ha='center', va='bottom')
                
                plt.tight_layout()
                save_path = os.path.join(self.dir_figs, output_plot)
                plt.savefig(save_path)
                print(f"✅ Plot saved: {save_path}")
                plt.show() 
            else:
                print("❌ No stars found within the margin.")
        else:
            print("No run*.hcm files found.")
 
    def make_file_list(self, list_name="all.lis"):
        print("\n" + "="*40 + "\nSTAGE 4: Generating File List\n" + "="*40)
        files = sorted(glob.glob(os.path.join(self.dir_save, "run*.hcm")))
        output_path = os.path.join(self.dir_save, list_name)
        with open(output_path, "w") as f:
            for file_path in files: f.write(os.path.basename(file_path) + "\n")
        print(f"✅ Created list file: {output_path} ({len(files)} files)")

    def genred(self):
        base_path = self.dir_save
        aperture_file = os.path.join(base_path, "ape.ape") 
        output_red = os.path.join(base_path, "reduce.red") 
        args = ["genred", aperture_file, output_red, "none", "none", "none", "none", "0", "none", "ultraspec-tnt"]
        
        print(f"Running genred...")
        try:
            scripts.genred(args)
            print(f"✅ Successfully created: {output_red}")
        except Exception as e:
            print(f"❌ Error running genred: {e}")
        
        if os.path.exists(output_red):
            print("\n--- File Created Successfully (Preview) ---")
            with open(output_red, 'r') as f:
                for _ in range(12): print(next(f), end='')

    def run_reduce_task(self, list_file="all.lis", config_file="reduce.red", log_file="all.log"):
        print("\n" + "="*40 + "\nSTAGE 5: HiPERCAM Reduction\n" + "="*40)
        original_cwd = os.getcwd()
        os.chdir(self.dir_save)
        
        reduce_args = ["reduce", "source=hf", f"flist={list_file}", "trim=False", 
                       f"rfile={config_file}", f"log={log_file}", "tkeep=1", 
                       "lplot=False", "implot=False"]

        try:
            print(f"🚀 Starting FULLY AUTOMATED Reduction...")
            scripts.reduce(reduce_args)
            print(f"🎉 SUCCESS! Log saved to: {log_file}")
        except Exception as e:
            print(f"❌ Reduction Failed: {e}")
        finally:
            os.chdir(original_cwd)

    # =========================================================================
    # ENSEMBLE PHOTOMETRY METHODS (Merged from AutoEnsembleCalibration)
    # =========================================================================

    def setup_ensemble_data(self, clean_df, raw_original, first_column: int = None):
        self.data = clean_df.copy()
        self.raw_original = raw_original.copy()
        self.first_column = int(first_column)
        n_tail = len(self.data.columns) - self.first_column
        self.number_stars = n_tail // 2
        return self.data

    def cal_instru_mag(self) -> pd.DataFrame:
        new_cols = {}
        for i in range(self.number_stars):
            flux_col = self.data.columns[self.first_column + 2 * i]
            err_col = self.data.columns[self.first_column + 2 * i + 1]
            flux = self.data[flux_col].astype(float)
            flux_err = self.data[err_col].astype(float)
            
            with np.errstate(divide="ignore", invalid="ignore"):
                mag = -2.5 * np.log10(np.where(flux > 0, flux, np.nan))
                mag_err = (2.5 / np.log(10)) * (flux_err / flux)
                mag_err = np.where(flux > 0, mag_err, np.nan)
            
            new_cols[f"instru_{i+1}"] = mag
            new_cols[f"err_instru_{i+1}"] = mag_err
            
        new_df = pd.DataFrame(new_cols, index=self.data.index)
        cols_to_drop = [c for c in new_df.columns if c in self.data.columns]
        if cols_to_drop:
            self.data = self.data.drop(columns=cols_to_drop)
        self.data = pd.concat([self.data, new_df], axis=1)
        return self.data

    def get_instrumental_dataframe(self, stars: Optional[list[int]] = None):
        if stars is None: stars = list(range(1, self.number_stars + 1))
        self.instru_cols = [f"instru_{s}" for s in stars]
        self.err_instru_cols = [f"err_instru_{s}" for s in stars]
        self.instrumental_data = self.data[self.instru_cols].copy()
        self.err_instrumental_data = self.data[self.err_instru_cols].copy()
        return self.instrumental_data, self.err_instrumental_data

    def _to_array(self, x):
        return x.to_numpy(dtype=float) if isinstance(x, (pd.DataFrame, pd.Series)) else np.asarray(x, dtype=float)

    def sum_mean_error(self, df, er=None, mode="sum", weighted=False, axis=0):
        arr = self._to_array(df)
        if mode == "mean":
            if weighted:
                err = self._to_array(er)
                with np.errstate(divide="ignore", invalid="ignore"):
                    w = 1.0 / err**2
                w[~np.isfinite(w)] = 0
                arr = np.where(np.isfinite(arr), arr, 0)
                num = np.sum(w * arr, axis=axis)
                sumw = np.sum(w, axis=axis)
                mean = np.divide(num, sumw, out=np.full_like(num, np.nan), where=sumw > 0)
                err_mean = np.where(sumw > 0, 1.0 / np.sqrt(sumw), np.nan)
                return mean, err_mean
            else:
                mean = np.nanmean(arr, axis=axis)
                if er is None: return mean, np.full_like(mean, np.nan)
                err = self._to_array(er)
                n = np.sum(np.isfinite(arr), axis=axis)
                err_mean = np.sqrt(np.nansum(err**2, axis=axis)) / np.sqrt(n)
                return mean, err_mean

    def ensemble(self, list_star=None, outdir=None, remark=None):
        self.cal_instru_mag()
        if list_star is None: list_star = list(range(1, self.number_stars + 1))
        self.list_stars = list_star
        self.get_instrumental_dataframe(list_star)

        mean_per_frame, err_per_frame = self.sum_mean_error(self.instrumental_data, er=self.err_instrumental_data, mode="mean", weighted=False, axis=1)
        ref_point, err_ref_point = self.sum_mean_error(mean_per_frame, er=err_per_frame, mode="mean", weighted=True, axis=0)
        self.sky_levels = -ref_point + mean_per_frame
        self.err_sky_levels = np.sqrt(err_per_frame**2 + err_ref_point**2)

        self.get_instrumental_dataframe(range(1, self.number_stars + 1))
        self.cor_sky_data = self.instrumental_data.subtract(self.sky_levels, axis=0)
        self.cor_sky_dataerr = np.sqrt(self.err_instrumental_data**2 + self.err_sky_levels[:, np.newaxis]**2)

        new_cor_cols = {}
        err_arr = np.asarray(self.cor_sky_dataerr) 
        for i, star in enumerate(range(1, self.number_stars + 1)):
            new_cor_cols[f"cor_instru_{star}"] = self.cor_sky_data.iloc[:, i]
            new_cor_cols[f"err_cor_instru_{star}"] = err_arr[:, i]
        
        new_df = pd.DataFrame(new_cor_cols, index=self.data.index)
        cols_to_drop = [c for c in new_df.columns if c in self.data.columns]
        if cols_to_drop:
            self.data = self.data.drop(columns=cols_to_drop)
        self.data = pd.concat([self.data, new_df], axis=1)

        if outdir:
            fig = self.plot_before_after_()
            self._auto_save(fig, outdir, remark)
        return self.data

    def plot_before_after_(self):
        list_star = self.list_stars
        instru_df, err_df = self.get_instrumental_dataframe(list_star)
        cor_df = self.cor_sky_data.loc[:, [f"instru_{s}" for s in list_star]]
        err_arr = np.asarray(self.cor_sky_dataerr)
        epochs = np.arange(len(instru_df))
        colors = plt.cm.tab10(np.linspace(0, 1, len(list_star)))

        fig = plt.figure(figsize=(12, 6))
        gs = gridspec.GridSpec(1, 2, figure=fig)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])

        for i, star in enumerate(list_star):
            col = colors[i % 10]
            raw_data = instru_df[f"instru_{star}"]
            ax1.errorbar(epochs, raw_data, yerr=err_df[f"err_instru_{star}"], fmt="o", ms=2, alpha=0.6, color=col)
            ax1.axhline(np.nanmean(raw_data), ls='--', lw=1, alpha=0.8, color=col)

            cor_data = cor_df[f"instru_{star}"]
            j = self.cor_sky_data.columns.get_loc(f"instru_{star}")
            ax2.errorbar(epochs, cor_data, yerr=err_arr[:, j], fmt="o", ms=2, alpha=0.6, label=f"Star {star}", color=col)
            ax2.axhline(np.nanmean(cor_data), ls='--', lw=1, alpha=0.8, color=col)

        ax1.set_title("Raw Data"); ax1.invert_yaxis(); ax1.grid(alpha=0.3)
        ax2.set_title("Corrected Data"); ax2.invert_yaxis(); ax2.grid(alpha=0.3); ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        return fig

    def _auto_save(self, fig, outdir, remark):
        # NOTE: Saving PNGs to self.dir_figs instead of results outdir for organization
        save_path = os.path.join(self.dir_figs, f"ensemble_plot_{remark}.png")
        fig.savefig(save_path)
        print(f"Saved Plot: {save_path}")
        plt.show() # Show the plot before continuing

    def save_hipercam_format(self, output_filename, used_stars, target_id, avg_rms, log_history):
        if self.raw_original is None:
            print("❌ Cannot save HiPERCAM format: Original raw data missing.")
            return

        print(f"--- Saving Calibrated Log: {output_filename} ---")
        final_df = self.raw_original.copy()
        global_cols = 7
        block_size = 16 

        for i in range(self.number_stars):
            star_id = i + 1
            col_name_mag = f"cor_instru_{star_id}"
            col_name_err = f"err_cor_instru_{star_id}"

            if col_name_mag in self.data.columns:
                cal_mag = self.data[col_name_mag].values
                cal_err = self.data[col_name_err].values
                new_counts = 10**(-0.4 * cal_mag)
                new_counts_err = new_counts * cal_err * (np.log(10) / 2.5)

                start_idx = global_cols + (i * block_size)
                counts_idx = start_idx + 8
                countse_idx = start_idx + 9

                final_df.iloc[:, counts_idx] = new_counts
                final_df.iloc[:, countse_idx] = new_counts_err
        
        header_str = "CCD nframe MJD MJDok Exptim mfwhm mbeta "
        for i in range(self.number_stars):
            s = i + 1
            header_str += f"x_{s} xe_{s} y_{s} ye_{s} fwhm_{s} fwhme_{s} beta_{s} betae_{s} counts_{s} countse_{s} sky_{s} skye_{s} nsky_{s} nrej_{s} cmax_{s} flag_{s} "
        
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        header_str += "\n# ========================================================\n"
        header_str += f"# ENSEMBLE PHOTOMETRY CALIBRATION LOG\n"
        header_str += f"# Generated: {now}\n"
        header_str += f"# Target ID: {target_id}\n"
        header_str += f"# Best Reference Stars: {sorted(used_stars)}\n"
        header_str += f"# Final Avg RMS of Best Stars {sorted(used_stars)}: {avg_rms:.6f}\n"
        header_str += "# ========================================================\n"
        header_str += "# OPTIMIZATION HISTORY:\n"
        for line in log_history: header_str += f"# {line}\n"
        header_str += "# ========================================================\n"
        
        with open(output_filename, 'w') as f:
            f.write(header_str)
            final_df.to_csv(f, sep=' ', index=False, header=False, float_format='%.6f', na_rep='NaN')
        print(f"✅ Saved calibrated log: {output_filename}")

    # --- MODIFIED: Optimize and Save for a Specific Target ID ---
    def optimize_references(self, target_id, min_stars=5, max_stars=10, initial_pool=None):
        if self.data is None: raise ValueError("Run setup_data() first.")

        if initial_pool is None:
            pool = list(range(1, self.number_stars + 1))
        else:
            pool = list(initial_pool)
            
        if target_id in pool: pool.remove(target_id) 

        print(f"\n" + "="*40 + f"\nSTAGE 6: Ensemble Optimization (Target {target_id})\n" + "="*40)
        
        log_buffer = [] 
        msg = f"Target: {target_id} | Initial Pool: {len(pool)} stars"
        print(msg)
        log_buffer.append(msg)

        history = []
        
        while len(pool) > min_stars:
            self.ensemble(list_star=pool, outdir=None)
            
            scatters = {}
            for star in pool:
                mags = self.data[f"cor_instru_{star}"]
                valid_mags = mags[np.isfinite(mags)]
                if len(valid_mags) > 2:
                    _, _, std = sigma_clipped_stats(valid_mags, sigma=3.0)
                    scatters[star] = std
                else:
                    scatters[star] = 99.9

            worst_star = max(scatters, key=scatters.get)
            worst_std = scatters[worst_star]
            current_rms = np.mean([val for val in scatters.values() if val < 99])
            
            history.append({'n': len(pool), 'avg_rms': current_rms, 'worst': worst_star, 'worst_rms': worst_std})

            if len(pool) <= max_stars:
                if worst_std < 1.1 * current_rms:
                    msg = f"Converged: Remaining stars are well balanced."
                    print(msg)
                    log_buffer.append(msg)
                    break
            
            msg = f"Pool: {len(pool):2d} | Avg RMS: {current_rms:.4f} | Removing Star {worst_star} (RMS: {worst_std:.4f})"
            print(msg)
            log_buffer.append(msg)
            pool.remove(worst_star)

        self.best_stars = pool
        
        # --- FINAL RUN ON BEST STARS ---
        self.ensemble(list_star=pool, outdir=self.dir_figs, remark=f"target_{target_id}_final")
        
        final_scatters = []
        for star in pool:
            mags = self.data[f"cor_instru_{star}"]
            valid_mags = mags[np.isfinite(mags)]
            if len(valid_mags) > 2:
                _, _, std = sigma_clipped_stats(valid_mags, sigma=3.0)
                final_scatters.append(std)
        
        best_stars_avg_rms = np.mean(final_scatters) if final_scatters else 0.0

        msg1 = f"OPTIMIZATION COMPLETE. Best {len(pool)} stars: {pool}"
        msg2 = f"Final Avg RMS of Best Stars {pool}: {best_stars_avg_rms:.6f}"
        
        print("\n✅ " + msg1)
        print("📊 " + msg2)
        log_buffer.append(msg1); log_buffer.append(msg2)
        
        self._plot_history(history)
        
        # --- SAVE UNIQUE FILE FOR EACH TARGET ---
        output_log = os.path.join(self.dir_results, f"calibrated_target_{target_id}.log")
        self.save_hipercam_format(output_log, used_stars=pool, target_id=target_id, avg_rms=best_stars_avg_rms, log_history=log_buffer)
        
        return pool

    def _plot_history(self, history):
        if not history: return
        n = [h['n'] for h in history]
        rms = [h['avg_rms'] for h in history]
        fig_hist = plt.figure(figsize=(6, 4))
        plt.plot(n, rms, 'o-', color='teal')
        plt.gca().invert_xaxis()
        plt.title("Optimization History")
        plt.xlabel("Number of Stars"); plt.ylabel("Ensemble RMS")
        plt.grid(True, ls='--')
        
        # Save historical plot
        save_path = os.path.join(self.dir_figs, "optimization_history.png")
        plt.savefig(save_path)
        print(f"Saved: {save_path}")
        plt.show() # Show plot
