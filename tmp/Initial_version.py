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
import re

# Astropy / Photutils
from astropy.visualization import ZScaleInterval
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder

# HiPERCAM specific imports
import hipercam as hcam
from hipercam import MCCD 
import hipercam.scripts as scripts

# ANSI Color Codes for terminal text formatting
BLUE = '\033[94m'
RED = '\033[91m'
RESET = '\033[0m'

# ==========================================
# HELPER FUNCTIONS
# ==========================================
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


class orderhcam:
    def __init__(self, save_path):
        self.save_path = os.path.abspath(save_path) 
        os.makedirs(self.save_path, exist_ok=True)

        # Paths to the .lis files generated after grabbing frames
        self.lisfilepath = {
            'bias': [None, None, None], # [0]:all, [1]:dark, [2]:flat
            'dark': [None, None],       # [0]:all, [1]:flat (dark-flat)
            'flat': [None],
            'data': [None],
        }

        # Absolute paths to the processed Master files (.hcm)
        self.masterfile = {
            'bias': [None, None, None],
            'dark': [None, None],
            'flat': [None],
        }

        # Subdirectory names for organizational clarity
        self.temps_path = {
            'bias': ['bias_all', 'bias_dark', 'bias_flat'],
            'dark': ['dark_all', 'dark_flat'],
            'flat': ['flat_all'],
            'data': ['data']
        }

    def grab(self, run_path, run, mode, **kwargs):
        """Extracts raw frames and converts them to .hcm format."""
        full_run_path = os.path.join(run_path, run)

        f1     = kwargs.get('f1', 1)
        f2     = kwargs.get('f2', 0)
        ndigit = kwargs.get('ndigit', 3)
        bias   = kwargs.get('bias', 'none')
        flat   = kwargs.get('flat', 'none')
        dark   = kwargs.get('dark', 'none')
        run_out = kwargs.get('run_out', run)
             
        args = [
            "nodef", 
            "source=ul", 
            f"run={full_run_path}", 
            f"first={f1}", 
            f"last={f2}",
            f"ndigit={ndigit}", 
            "trim=False", 
            f"bias={bias}", 
            f"flat={flat}", 
            f"dark={dark}", 
            "fmap=none", 
            f"output={run_out}"
        ]

        old_files = glob.glob(f"{run_out}*[0-9]*.hcm")
        for f in old_files:
            try:
                os.remove(f)
            except OSError:
                pass
                
        try:
            print(f"{BLUE}#### Grabbing {run} ({mode}) {full_run_path} ####{RESET}")
            scripts.grab(args)
        except Exception as e:
            print(f"{RED}Error processing {run}: {e}{RESET}")

    def grab_data(self, run_path, run, idx=0, **kwargs):
        """Extracts Science Frames and applies calibration on the fly."""
        full_run_path = os.path.join(run_path, run)
        f1 = kwargs.get('f1', 1)
        f2 = kwargs.get('f2', 0)
        run_out = kwargs.get('run_out', run)
        ndigit = kwargs.get('ndigit', 3)
        chosen_bias = self.get_bias_for('data')
        chosen_dark = self.get_dark_for('data')
        chosen_flat = self.masterfile['flat'][0] if self.masterfile['flat'][0] else 'none'

        args = [
            "nodef", "source=ul", 
            f"run={full_run_path}", 
            f"first={f1}", 
            f"last={f2}", 
            f"ndigit={ndigit}",
            "trim=False", 
            f"bias={chosen_bias}", 
            f"flat={chosen_flat}", 
            f"dark={chosen_dark}", 
            "fmap=none", 
            f"output={run_out}"
        ]
        
        old_files = glob.glob(f"{run_out}*[0-9]*.hcm")
        for f in old_files:
            try:
                os.remove(f)
            except OSError:
                pass
                
        try:
            print(f"{BLUE}#### Grabbing Data {run} {full_run_path} ####{RESET}")
            print(f"{BLUE}Calibration -> Bias: {chosen_bias}, Dark: {chosen_dark}, Flat: {chosen_flat}{RESET}")
            scripts.grab(args)
        except Exception as e:
            print(f"{RED}Error processing data {run}: {e}{RESET}")
            
    def get_bias_for(self, mode):
        """Selects the best master bias based on priority."""
        bias_map = {'all': 0, 'dark': 1, 'flat': 2}
        target_idx = bias_map.get(mode)

        if target_idx is not None and self.masterfile['bias'][target_idx] is not None:
            return self.masterfile['bias'][target_idx]
        if self.masterfile['bias'][0] is not None:
            return self.masterfile['bias'][0]
        return 'none'
    
    def get_dark_for(self, mode):
        """Selects the best master dark based on priority."""
        dark_map = {'all': 0, 'flat': 1}
        target_idx = dark_map.get(mode)

        if target_idx is not None and self.masterfile['dark'][target_idx] is not None:
            return self.masterfile['dark'][target_idx]
        if self.masterfile['dark'][0] is not None:
            return self.masterfile['dark'][0]
        return 'none'        

    def bias(self, idx=0, **kwargs):
        """Creates a Master Bias."""
        biaslis = kwargs.get('bias_lis', self.lisfilepath['bias'][idx])
        run_out = kwargs.get('run_out', f'mbias_{idx}')
        sigma = kwargs.get('sigma', 3)
        args = [
            "nodef", "source=hf", f"flist={biaslis}",
            f"sigma={sigma}", "plot=no", f"output={run_out}"
        ]
        try:
            print(f"{BLUE}--- Making Master Bias [{idx}] ---{RESET}")
            scripts.makebias(args)
            full_master_path = os.path.abspath(f'{run_out}.hcm')
            self.masterfile['bias'][idx] = full_master_path
        except Exception as e:
            print(f"{RED}Error processing bias: {e}{RESET}")

    def dark(self, idx=0, **kwargs):
        """Creates a Master Dark with bias subtraction."""
        darklis = kwargs.get('dark_lis', self.lisfilepath['dark'][idx])
        run_out = kwargs.get('run_out', f'mdark_{idx}')
        chosen_bias = self.get_bias_for('dark')

        args = [
            "nodef", "source=hf", f"flist={darklis}",
            f"bias={chosen_bias}", "plot=no", f"output={run_out}"
        ]
        try:
            print(f"{BLUE}--- Processing Master Dark [{idx}] ---{RESET}")
            print(f"{BLUE}Using Master Bias: {chosen_bias}{RESET}")
            scripts.makedark(args)
            full_master_path = os.path.abspath(f'{run_out}.hcm')
            self.masterfile['dark'][idx] = full_master_path
        except Exception as e:
            print(f"{RED}Error processing dark: {e}{RESET}")

    def flat(self, idx=0, **kwargs):
        """Creates a Master Flat with bias and dark subtraction."""
        flatlis = kwargs.get('flat_lis', self.lisfilepath['flat'][idx])
        run_out = kwargs.get('run_out', f'mflat_{idx}')
        chosen_bias = self.get_bias_for('flat')
        chosen_dark = self.get_dark_for('flat')
        
        args = [
            "nodef", "source=hf", f"flist={flatlis}",
            f"bias={chosen_bias}", f"dark={chosen_dark}",
            "ngroup=5", "ccd=0",
            f"lower={kwargs.get('lower', '3000 3000 3000')}",
            f"upper={kwargs.get('upper', '49000 29000 27000')}",
            f"output={run_out}"
        ]
        try:
            print(f"{BLUE}--- Processing Master Flat [{idx}] ---{RESET}")
            print(f"{BLUE}Using Master Bias: {chosen_bias}, Dark: {chosen_dark}{RESET}")
            scripts.makeflat(args)
            full_master_path = os.path.abspath(f'{run_out}.hcm')
            self.masterfile['flat'][idx] = full_master_path
        except Exception as e:
            print(f"{RED}Error processing flat: {e}{RESET}")    

    def generate_lis_file(self, run_list, mode, index):
        """Generates the .lis file needed for master creation."""
        list_filename = f'{mode}_{index}.lis'
        full_lis_path = os.path.abspath(list_filename)
        
        found_files = []
        for r in run_list:
            found_files.extend(sorted(glob.glob(f"{r}*[0-9]*.hcm")))
            
        found_files = list(set(found_files))
        
        if found_files:
            with open(list_filename, "w") as f:
                for filename in sorted(found_files):
                    f.write(os.path.abspath(filename) + "\n")
            self.lisfilepath[mode][index] = full_lis_path
            print(f"{BLUE}Created list file: {full_lis_path} with {len(found_files)} files.{RESET}")
        else:
            print(f"{RED}Warning: No .hcm files found to create list for {mode}_{index}{RESET}")

    def plot_lis(self, mode, index=0, ccd_num='1'):
        """Plots the .hcm files listed in a specific .lis file."""
        
        lis_file = self.lisfilepath[mode][index]
        if not lis_file or not os.path.exists(lis_file):
            print(f"{RED}Error: No .lis file found for mode '{mode}' at index {index}{RESET}")
            return

        with open(lis_file, 'r') as f:
            hcm_files = [line.strip() for line in f.readlines() if line.strip()]

        ntot = len(hcm_files)
        if ntot == 0:
            print(f"{RED}Error: The .lis file {lis_file} is empty.{RESET}")
            return

        print(f"{BLUE}--- Plotting {ntot} files from {os.path.basename(lis_file)} ---{RESET}")
        ncols = 5
        nrows = math.ceil(ntot / ncols)
        
        fig, axs = plt.subplots(nrows, ncols, figsize=(6*ncols, 6*nrows), 
                                gridspec_kw={'hspace': 0.05, 'wspace': 0.05}, facecolor='whitesmoke')
        
        if ntot == 1:
            axs_flat = [axs]
        else:
            axs_flat = np.atleast_1d(axs).flatten()
            
        zscale = ZScaleInterval()

        for i, file_path in enumerate(hcm_files):
            ax = axs_flat[i]
            try:
                mccd = hcam.MCCD.read(file_path)         
                data = mccd[ccd_num]['1'].data
                vmin, vmax = zscale.get_limits(data)

                ax.imshow(data, cmap='gray_r', vmin=vmin, vmax=vmax, origin='lower')
                file_name = os.path.basename(file_path)
                ax.set_title(file_name, fontsize=10)
                
                # Remove ticks and labels but keep the frame
                ax.set_xticks([])
                ax.set_yticks([])

                # Draw a black border around the image
                for spine in ax.spines.values():
                    spine.set_visible(True)
                    spine.set_color('black')
                    spine.set_linewidth(1.5)
                
            except Exception as e:
                ax.set_title("Error", color='red', fontsize=10)
                ax.axis('off')
                print(f"{RED}Error plotting {file_path}: {e}{RESET}")
                
        for j in range(ntot, len(axs_flat)):
            axs_flat[j].axis('off')
        x = os.path.basename(lis_file)
        figs_dir = os.path.join(self.save_path, 'figs')
        os.makedirs(figs_dir, exist_ok=True) # Create the folder if it doesn't exist
        
        save_file = os.path.join(figs_dir, f'{x}.png')
        plt.savefig(save_file)
        plt.show()
        
        
    def plot_median(self, mode, index=0, ccd_num='1'):
        """Calculates and plots the median pixel value for each frame in a .lis file."""
        lis_file = self.lisfilepath[mode][index]
        if not lis_file or not os.path.exists(lis_file):
            print(f"{RED}Error: No .lis file found for mode '{mode}' at index {index}{RESET}")
            return

        with open(lis_file, 'r') as f:
            hcm_files = [line.strip() for line in f.readlines() if line.strip()]

        if not hcm_files:
            print(f"{RED}Error: The .lis file {lis_file} is empty.{RESET}")
            return

        print(f"{BLUE}--- Calculating Medians for {len(hcm_files)} files from {os.path.basename(lis_file)} ---{RESET}")
        
        frame_numbers = []
        median_levels = []

        for i, file_path in enumerate(hcm_files):
            try:
                mccd = hcam.MCCD.read(file_path)
                data = mccd[ccd_num]['1'].data
                
                # Calculate the median, ignoring NaNs or dead pixels
                median_val = np.nanmedian(data)
                
                frame_numbers.append(i + 1)
                median_levels.append(median_val)
            except Exception as e:
                print(f"{RED}Error reading {file_path}: {e}{RESET}")

        # --- Plotting the results ---
        plt.figure(figsize=(10, 5), facecolor='white')
        plt.plot(frame_numbers, median_levels, marker='o', linestyle='-', color='b', markersize=5)
        
        plt.title(f"Median Pixel Level vs. Frame Number ({mode.upper()})", fontsize=14, fontweight='bold')
        plt.xlabel("Frame Sequence", fontsize=12)
        plt.ylabel("Median Pixel Value (ADU)", fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        x = os.path.basename(lis_file)
        x = os.path.join(self.save_path, 'figs', f'median_{x}.png')
        plt.savefig(x)
        plt.show()


    
    def run(self, input_run, default_raw_dir='', plot_show=False):
        initial_path = os.getcwd()
        os.chdir(self.save_path)

        # Helper to extract start and end frames for both 1D and 2D lists
        def get_start_end(frames_cfg, index):
            if isinstance(frames_cfg[0], list):
                return frames_cfg[index][0], frames_cfg[index][1]
            return frames_cfg[0], frames_cfg[1]

        # ------------------- BIAS run -------------------
        bias_index_map = {'all': 0, 'dark': 1, 'flat': 2}
        bias_configs = input_run.get('bias', {})
        for b_type, b_cfg in bias_configs.items():
            if b_cfg is not None:
                idx = bias_index_map[b_type]
                # If raw_dir is not in the config, use the default_raw_dir instead
                current_raw_dir = b_cfg.get('raw_dir', default_raw_dir) 
                try:
                    target_dir = os.path.join(self.save_path, self.temps_path['bias'][idx])
                    os.makedirs(target_dir, exist_ok=True)
                    os.chdir(target_dir)       
                    
                    runs_processed = []
                    for i, r in enumerate(b_cfg['runs']):
                        f_start, f_end = get_start_end(b_cfg['frames'], i)
                        self.grab(run_path=current_raw_dir, run=r, mode='bias', run_out=r, f1=f_start, f2=f_end)
                        runs_processed.append(r)
                        
                    self.generate_lis_file(run_list=runs_processed, mode='bias', index=idx)                        
                    self.bias(idx=idx)
                except Exception as e:
                    print(f'{RED}Error processing bias {b_type}: {e}{RESET}')
                finally:
                    os.chdir(self.save_path) 

        # ------------------- DARK run -------------------
        dark_index_map = {'all': 0, 'flat': 1}
        dark_configs = input_run.get('dark', {})
        for d_type, d_cfg in dark_configs.items():
            if d_cfg is not None:
                idx = dark_index_map.get(d_type, 0)
                current_raw_dir = d_cfg.get('raw_dir', default_raw_dir)
                try:
                    target_dir = os.path.join(self.save_path, self.temps_path['dark'][idx])
                    os.makedirs(target_dir, exist_ok=True)
                    os.chdir(target_dir)       
                    
                    runs_processed = []
                    for i, r in enumerate(d_cfg['runs']):
                        f_start, f_end = get_start_end(d_cfg['frames'], i)
                        self.grab(run_path=current_raw_dir, run=r, mode='dark', run_out=r, f1=f_start, f2=f_end)
                        runs_processed.append(r)
                        
                    self.generate_lis_file(run_list=runs_processed, mode='dark', index=idx)                        
                    self.dark(idx=idx)
                except Exception as e:
                    print(f'{RED}Error processing dark {d_type}: {e}{RESET}')
                finally:
                    os.chdir(self.save_path) 

        # ------------------- FLAT run -------------------
        flat_index_map = {'all': 0}
        flat_configs = input_run.get('flat', {})
        for f_type, f_cfg in flat_configs.items():
            if f_cfg is not None:
                idx = flat_index_map.get(f_type, 0)
                current_raw_dir = f_cfg.get('raw_dir', default_raw_dir)
                try:
                    target_dir = os.path.join(self.save_path, self.temps_path['flat'][idx])
                    os.makedirs(target_dir, exist_ok=True)
                    os.chdir(target_dir)       
                    
                    runs_processed = []
                    for i, r in enumerate(f_cfg['runs']):
                        f_start, f_end = get_start_end(f_cfg['frames'], i)
                        self.grab(run_path=current_raw_dir, run=r, mode='flat', run_out=r, f1=f_start, f2=f_end)
                        runs_processed.append(r)
                        
                    self.generate_lis_file(run_list=runs_processed, mode='flat', index=idx)                        
                    self.flat(idx=idx)
                except Exception as e:
                    print(f'{RED}Error processing flat {f_type}: {e}{RESET}')
                finally:
                    os.chdir(self.save_path)

        # ------------------- DATA run -------------------
        data_configs = input_run.get('data', {})
        for d_type, d_cfg in data_configs.items():
            if d_cfg is not None:
                idx = 0 
                current_raw_dir = d_cfg.get('raw_dir', default_raw_dir)
                try:
                    target_dir = os.path.join(self.save_path, self.temps_path['data'][idx])
                    os.makedirs(target_dir, exist_ok=True)
                    os.chdir(target_dir)       
                    
                    runs_processed = []
                    for i, r in enumerate(d_cfg['runs']):
                        f_start, f_end = get_start_end(d_cfg['frames'], i)
                        self.grab_data(run_path=current_raw_dir, run=r, idx=idx, run_out=r, f1=f_start, f2=f_end)
                        runs_processed.append(r)
                        
                    self.generate_lis_file(run_list=runs_processed, mode='data', index=idx)
                except Exception as e:
                    print(f'{RED}Error processing data {d_type}: {e}{RESET}')
                finally:
                    os.chdir(self.save_path)
        
        if plot_show is True:
            self.plot_all_available()
            # self.check_bad_pixels(ccd_num='1') 
        os.chdir(initial_path)     
    def plot_all_available(self, ccd_num='1'):
        """Plots all available .lis files and median data automatically."""
        print(f"\n{BLUE}=== Initiating Automated Diagnostic Plots for CCD {ccd_num} ==={RESET}")
        for mode, path_list in self.lisfilepath.items():
            for index, lis_file in enumerate(path_list):
                if lis_file is not None and os.path.exists(lis_file):
                    self.plot_lis(mode=mode, index=index, ccd_num=ccd_num)
                    self.plot_median(mode=mode, index=index, ccd_num=ccd_num)

    def check_bad_pixels(self, ccd_num='1', hot_sigma=3.0, dead_fraction=0.8,
                         MARGIN_LEFT=15, MARGIN_RIGHT=15, MARGIN_BOTTOM=5, MARGIN_TOP=27):
        """
        ตรวจสอบหา Hot/Dead Pixels เฉพาะในพื้นที่ Safe Zone (ตัดขอบภาพทิ้งตาม Margin)
        """
        print(f"\n{BLUE}=== Checking for Bad Pixels (CCD: {ccd_num} | Safe Zone Only) ==={RESET}")
        
        master_dark = self.masterfile['dark'][0] or self.masterfile['dark'][1]
        master_flat = self.masterfile['flat'][0]
        
        hot_pixels, dead_pixels = None, None
        num_hot, num_dead = 0, 0
        h, w = 0, 0 # เก็บขนาดภาพ
        
        # --- หา HOT PIXELS (ถ้ามีไฟล์ Dark) ---
        if master_dark and os.path.exists(master_dark):
            try:
                mccd_dark = hcam.MCCD.read(master_dark)
                data_dark = mccd_dark[ccd_num]['1'].data
                h, w = data_dark.shape
                
                # สร้าง Mask สำหรับ Safe Zone
                valid_mask = np.zeros((h, w), dtype=bool)
                valid_mask[MARGIN_BOTTOM : h - MARGIN_TOP, MARGIN_LEFT : w - MARGIN_RIGHT] = True
                
                # คำนวณสถิติเฉพาะข้อมูลใน Safe Zone
                valid_data_dark = data_dark[valid_mask]
                _, med_dark, std_dark = sigma_clipped_stats(valid_data_dark, sigma=3.0)
                
                hot_threshold = med_dark + (hot_sigma * std_dark)
 
                hot_pixels = (data_dark > hot_threshold) & valid_mask
                num_hot = np.sum(hot_pixels)
                print(f"-> Hot Threshold: > {hot_threshold:.2f} ADU | Found {RED}{num_hot}{RESET} Hot Pixels")
            except Exception as e: print(f"{RED}Error reading Dark: {e}{RESET}")
        else:
            print(f"{RED}Master Dark not found. Skipping hot pixel check.{RESET}")

        # --- หา DEAD PIXELS (ถ้ามีไฟล์ Flat) ---
        if master_flat and os.path.exists(master_flat):
            try:
                mccd_flat = hcam.MCCD.read(master_flat)
                data_flat = mccd_flat[ccd_num]['1'].data
                h, w = data_flat.shape
                
                # สร้าง Mask สำหรับ Safe Zone
                valid_mask = np.zeros((h, w), dtype=bool)
                valid_mask[MARGIN_BOTTOM : h - MARGIN_TOP, MARGIN_LEFT : w - MARGIN_RIGHT] = True
                
                # คำนวณสถิติเฉพาะข้อมูลใน Safe Zone
                valid_data_flat = data_flat[valid_mask]
                med_flat = np.nanmedian(valid_data_flat)
                
                dead_threshold = med_flat * dead_fraction
                # เลือกเฉพาะพิกเซลที่มืดเกินไป และ "ต้องอยู่ใน Safe Zone"
                dead_pixels = (data_flat < dead_threshold) & valid_mask
                num_dead = np.sum(dead_pixels)
                print(f"-> Dead Threshold: < {dead_threshold:.2f} ADU | Found {BLUE}{num_dead}{RESET} Dead Pixels")
            except Exception as e: print(f"{RED}Error reading Flat: {e}{RESET}")
        else:
            print(f"{RED}Master Flat not found. Skipping dead pixel check.{RESET}")


        if (hot_pixels is not None) or (dead_pixels is not None):
            fig, ax = plt.subplots(figsize=(8, 8), facecolor='whitesmoke')
            
            # ดึงขนาดภาพอ้างอิง
            ref_data = np.zeros((h, w), dtype=float) 
            
            if hot_pixels is not None: ref_data[hot_pixels] = 1.0
            if dead_pixels is not None: ref_data[dead_pixels] = -1.0
            
            cmap = plt.cm.colors.ListedColormap(['blue', 'whitesmoke', 'red'])
            ax.imshow(ref_data, origin='lower', cmap=cmap, vmin=-1, vmax=1)
            
            if h > 0 and w > 0:
                rect = patches.Rectangle((MARGIN_LEFT, MARGIN_BOTTOM), 
                                         w - MARGIN_LEFT - MARGIN_RIGHT, 
                                         h - MARGIN_BOTTOM - MARGIN_TOP,
                                         linewidth=2, edgecolor='orange', facecolor='none', linestyle='--')
                ax.add_patch(rect)
            
            ax.set_title(f"Bad Pixel Map (CCD {ccd_num} | Safe Zone)\nRed = Hot ({num_hot}), Blue = Dead ({num_dead})", 
                         fontsize=14, fontweight='bold')
            ax.set_xlabel("X (Pixels)"); ax.set_ylabel("Y (Pixels)")
            
            figs_dir = os.path.join(self.save_path, 'figs')
            os.makedirs(figs_dir, exist_ok=True)
            plt.savefig(os.path.join(figs_dir, f'bad_pixel_map_ccd{ccd_num}.png'), dpi=150, bbox_inches='tight')
            plt.show()
            plt.close(fig)



            
class Reduction:
    def __init__(self, save_path):
        self.save_path = save_path
        self.dir_figs = os.path.join(save_path,'figs')
        os.makedirs(self.dir_figs, exist_ok= True)
        # Base directory for saving files
        os.chdir(self.save_path)
        
    def setaper(self, list_file, ccd_label='1', win_label='1', SIGMA_THRESHOLD=1.5,
                output_plot="detection_labeled.png", output_ape_name="ape.ape",
                MARGIN_LEFT=15, MARGIN_RIGHT=15, MARGIN_BOTTOM=5, MARGIN_TOP=27,
                SKIP_BRIGHTEST=5, R_TARG=None, R_SKY1=16, R_SKY2=24, frame=5):
        """
        Detects stars from the first image in a .lis file and creates 
        a HiPERCAM aperture (.ape) file.
        """
        self.aperture = output_ape_name
        self.list_file= list_file
        # 1. Load and check the .lis file
        if not list_file or not os.path.exists(list_file):
            print(f"{RED}Error: List file '{list_file}' not found.{RESET}")
            return

        with open(list_file, 'r') as f:
            hcm_files = [line.strip() for line in f.readlines() if line.strip()]

        if not hcm_files:
            print(f"{RED}Error: The .lis file is empty.{RESET}")
            return
        
        # Select the first file as the reference frame
        target_file = hcm_files[frame]
        print(f"{BLUE}Analyzing reference file: {os.path.basename(target_file)}{RESET}")
        
        try:
            mccd = hcam.MCCD.read(target_file)
            window_obj = mccd[ccd_label][win_label]
            data = window_obj.data
            binning = getattr(window_obj, 'xbin', 1)
            self.binning = binning 
            h, w = data.shape
        except Exception as e:
            print(f"{RED}Critical Error reading HCM file: {e}{RESET}")
            return
        
        # 2. Source Detection logic
        print(f"Detecting sources with threshold {SIGMA_THRESHOLD} sigma...")
        _, median, std = sigma_clipped_stats(data, sigma=3)
        daofind = DAOStarFinder(fwhm=4.0, threshold=SIGMA_THRESHOLD * std)
        sources = daofind(data - median)
        
        if sources:
            # Sort by flux (brightest first)
            sources.sort('flux')
            sources = sources[::-1]

            # Filter stars based on defined margins
            mask = ((sources['xcentroid'] > MARGIN_LEFT) & 
                    (sources['xcentroid'] < w - MARGIN_RIGHT) & 
                    (sources['ycentroid'] > MARGIN_BOTTOM) & 
                    (sources['ycentroid'] < h - MARGIN_TOP))
            sources = sources[mask]
    
            # Skip the N brightest stars (useful for saturated sources)
            print(f"Skipping the first {SKIP_BRIGHTEST} brightest stars.")
            sources = sources[SKIP_BRIGHTEST:]

            sources['fwhm'] = np.sqrt(sources['npix'] / np.pi) * 2.35
            print(f"Total usable stars detected: {len(sources)}")
            self.num_aps = len(sources)  
            # 3. Construct HiPERCAM Aperture File Structure
            ccd_aps = []
            for i, source in enumerate(sources):
                x, y, fwhm = source['xcentroid'], source['ycentroid'], source['fwhm']
                ap_id = str(i + 1)
                is_ref = (i == 0) # The first star in the remaining list becomes the reference
                if R_TARG  is None : R_TARG = fwhm * 1
                # Link all comparison stars to Star 1 for group tracking
                link_to = "" if is_ref else "1"
                
                ccd_aps.append([ap_id, {
                    "Comment": "hipercam.Aperture", 
                    "x": float(x) * binning, 
                    "y": float(y) * binning,
                    "rtarg": R_TARG * binning, 
                    "rsky1": (R_TARG+4) * binning, 
                    "rsky2": (R_TARG+8) * binning,
                    "ref": is_ref, 
                    "compo": False, 
                    "mask": [], 
                    "extra": [], 
                    "link": link_to  #link_to 
                }])

            # HiPERCAM JSON schema requires this specific nested list format
            inner_ccd_structure = ["hipercam.CcdAper"] + ccd_aps
            ape_json = ["hipercam.MccdAper", [ccd_label, inner_ccd_structure]]
            
            # Save the .ape file
            output_path = os.path.join(self.save_path, output_ape_name)
            with open(output_path, 'w') as f:
                json.dump(ape_json, f, indent=2)
            print(f" Aperture file created: {output_path}")

            # 4. Generate the Verification Plot
            fig, ax = plt.subplots(figsize=(10, 10))
            zscale = ZScaleInterval()
            vmin, vmax = zscale.get_limits(data)
            ax.imshow(data, origin='lower', cmap='Greys', vmin=vmin, vmax=vmax)
            
            # Draw the safe margin box (dashed yellow line)
            rect = patches.Rectangle((MARGIN_LEFT, MARGIN_BOTTOM), 
                                     w - MARGIN_LEFT - MARGIN_RIGHT, 
                                     h - MARGIN_BOTTOM - MARGIN_TOP,
                                     linewidth=2, edgecolor='yellow', facecolor='none', linestyle='--')
            ax.add_patch(rect)
            
            # Draw circles and numbers for all detected apertures
            for i, source in enumerate(sources):
                x, y = source['xcentroid'], source['ycentroid']
                # Green circle for target aperture
                ax.add_patch(patches.Circle((x, y), R_TARG, edgecolor='lime', facecolor='none', lw=0.5))
                ax.add_patch(patches.Circle((x, y), R_SKY1, edgecolor='lime', facecolor='none', lw=0.5))
                ax.add_patch(patches.Circle((x, y), R_SKY2, edgecolor='lime', facecolor='none', lw=0.5))
                # Number label
                ax.text(x, y + R_SKY2 + 2, str(i + 1), color='cyan', 
                        fontsize=12, fontweight='bold', ha='center')
            
            plt.tight_layout()
            plt.savefig(os.path.join(self.dir_figs, output_plot))
            print(f" Verification plot saved: {output_plot}")
            plt.show()
        else:
            print(f"{RED}No stars found within defined margins.{RESET}")
    
    
    
    def genred(self, aperture="ape.ape", output_red="reduce.red", **kwargs):
        """
        Runs the HiPERCAM 'genred' command to create a reduction control file (.red).
        """
        # Define paths
        args = ["genred", self.aperture, output_red, "none", "none", "none", "none", "0", "none", "ultraspec-tnt"]

        self.red_file = output_red
        print(f"{BLUE}#### Generating Reduction File (.red) ####{RESET}")
        try:
            # Using the scripts module imported from hipercam
            scripts.genred(args)
            print(f" Successfully created: {output_red}")
        except Exception as e:
            print(f"{RED}Error running genred: {e}{RESET}")
        if kwargs:
            with open(output_red, 'r') as f:
                content = f.read()
            for key, value in kwargs.items():
                pattern = rf"^{key}\s*=\s*.*"
                replacement = f"{key} = {value}"
                
                if re.search(pattern, content, re.MULTILINE):
                    content = re.sub(pattern, replacement, content, flags=re.MULTILINE)
                    print(f"   [Override] {key} -> {value}")
                else:
                    print(f"   [Warning] Key '{key}' not found.")
            with open(output_red, 'w') as f:
                f.write(content)

                
        # Preview the file
        if os.path.exists(output_red):
            print(f"\n{BLUE}--- Reduction File Preview (First 12 lines) ---{RESET}")
            with open(output_red, 'r') as f:
                for _ in range(12):
                    try: print(next(f), end='')
                    except StopIteration: break
    # def modify_red(self, red_file="reduce.red", **kwargs):
    #     if not os.path.exists(red_file):
    #         print(f"{RED}Error: File '{red_file}' not found.{RESET}")
    #         return

    #     print(f"{BLUE}#### Modifying {red_file} ####{RESET}")
        
    #     with open(red_file, 'r') as f:
    #         lines = f.readlines()

    #     new_lines = []
    #     modified_keys = set()

    #     for line in lines:
    #         updated = False
           
    #         for key, value in kwargs.items():
                
    #             pattern = rf"^{key}\s*=\s*"
    #             if re.match(pattern, line.strip()):
    #                 comment = ""
    #                 if "#" in line:
    #                     comment = " # " + line.split("#", 1)[1].strip()
                    
    #                 new_lines.append(f"{key} = {value}{comment}\n")
    #                 modified_keys.add(key)
    #                 updated = True
    #                 break
            
    #         if not updated:
    #             new_lines.append(line)

    #     with open(red_file, 'w') as f:
    #         f.writelines(new_lines)

    #     for k in kwargs.keys():
    #         if k in modified_keys:
    #             print(f"   [Modified] {k} -> {kwargs[k]}")
    #         else:
    #             print(f"   [Skip] Key '{k}' not found in file (No changes made).")

    #     print(f" {red_file} updated successfully.")




    def modify_red(self, red_file="reduce.red", target_section=None, **kwargs):
        if not os.path.exists(red_file):
            print(f"{RED}Error: File '{red_file}' not found.{RESET}")
            return

        print(f"{BLUE}#### Modifying {red_file} [{target_section if target_section else 'All'}] ####{RESET}")
        
        with open(red_file, 'r') as f:
            lines = f.readlines()

        new_lines = []
        modified_keys = set()
        current_section = None

        for line in lines:
            strip_line = line.strip()
            
            # ตรวจจับ Section ปัจจุบัน
            if strip_line.startswith('[') and strip_line.endswith(']'):
                current_section = strip_line[1:-1].strip().lower()

            updated = False
            # สร้าง Regex ที่ยืดหยุ่น: ข้ามช่องว่างหน้าบรรทัด และเช็ค Key ให้ตรงตัวเป๊ะ (\b)
            # รองรับทั้งคีย์ที่เป็นตัวหนังสือ (เช่น fit_method) และตัวเลข (เช่น 1)
            keys_pattern = '|'.join([re.escape(str(k)) for k in kwargs.keys()])
            pattern = rf"^\s*\b({keys_pattern})\b\s*="
            
            match = re.match(pattern, line)
            if match:
                key = match.group(1)

                target = target_section.lower() if target_section else None
                
                if target is None or current_section == target:
                    value = kwargs[key]
                    comment = ""
                    if "#" in line:
                        comment = " # " + line.split("#", 1)[1].strip()
                    
                    # เขียนใหม่ใน Format: key = value (เว้นวรรคให้ถูกต้อง)
                    new_lines.append(f"{key} = {value}{comment}\n")
                    modified_keys.add(key)
                    updated = True
            
            if not updated:
                new_lines.append(line)

        with open(red_file, 'w') as f:
            f.writelines(new_lines)

        for k in kwargs.keys():
            if k in modified_keys:
                print(f"   [Modified] {k}")
            else:
                print(f"   [Skip/Fail] {k} (Check section or key name)")

        

    def reduce(self, log_file="all.log"):
        if not self.list_file or not self.red_file:
            print(f"{RED}Error: Missing .lis or .red file!{RESET}"); return

        os.chdir(self.save_path)
        reduce_args = ["reduce", "source=hf", f"flist={self.list_file}", "trim=False", 
                       f"rfile={self.red_file}", f"log={log_file}", "tkeep=1", 
                       "lplot=False", "implot=False"]
        try:
            print(f" Starting FULLY AUTOMATED Reduction...")
            scripts.reduce(reduce_args)
            print(f" SUCCESS! Log saved to: {log_file}")
        except Exception as e:
            print(f" Reduction Failed: {e}")

    
    def plot_with_log(self, log_filename, ccd_num='1'):
        df_log, n_aps = read_hipercam_log(log_filename)
        if df_log is None:
            print(f"{RED}Error: Cannot read log file.{RESET}"); return

        lis_file = self.list_file
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
        plt.savefig(f'reduction_check_{ccd_num}.png', dpi=200)
        plt.show()

    def plot_with_zoom(self, log_filename, ccd_num='1', zoom_box=15):
        """
        Plots a zoomed-in cutout of EVERY aperture for EACH frame.
        Includes X and Y marginal profile histograms for each star.
        """
     
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        df_log, n_aps = read_hipercam_log(log_filename)
        if df_log is None:
            print(f"{RED}Error: Cannot read log file.{RESET}"); return

        lis_file = self.list_file
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
                        # สร้างแกนกราฟขนาด 25% แปะที่ขอบ
                        ax_histx = divider.append_axes("top", size="25%", pad=0.05, sharex=ax)
                        ax_histy = divider.append_axes("right", size="25%", pad=0.05, sharey=ax)
                        
                        # ดึงข้อมูล Pixel ในกรอบที่ซูมเพื่อมาคำนวณกราฟ
                        x_start = max(0, int(x - zoom_box))
                        x_end = min(data.shape[1], int(x + zoom_box))
                        y_start = max(0, int(y - zoom_box))
                        y_end = min(data.shape[0], int(y + zoom_box))
                        cutout = data[y_start:y_end, x_start:x_end]
                        
                        if cutout.size > 0:
                            # 1. หาค่า Median ของรูปที่ซูม เพื่อใช้เป็นตัวแทนของพื้นหลัง (Sky Background)
                            local_median = np.nanmedian(cutout)
                            
                            # 2. นำข้อมูลมาลบด้วย Median (Subtract with median)
                            cutout_subbed = cutout - local_median
                            
                            x_range = np.arange(x_start, x_end)
                            y_range = np.arange(y_start, y_end)
                            
                            # หาค่าเฉลี่ยแสงตามแนวแกนจากข้อมูลที่ลบพื้นหลังแล้ว
                            profile_x = np.nanmean(cutout_subbed, axis=0)
                            profile_y = np.nanmean(cutout_subbed, axis=1)
                            
                            # 3. พลอตกราฟสไตล์ Histogram (ใช้ drawstyle='steps-mid')
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
                
                # เซฟไฟล์ (Adding bbox_inches='tight' helps cleanly crop the saved image)
                save_file = os.path.join(self.dir_figs, f'zoom_F{i+1}_{frame_name}.png')
                plt.savefig(save_file, dpi=150, bbox_inches='tight')
                plt.show() 
                plt.close(fig)

            except Exception as e:
                print(f"{RED}Error plotting frame {i}: {e}{RESET}")