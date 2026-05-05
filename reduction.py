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
import hipercam.spooler as spooler
import hipercam.mpl as hmpl
import hipercam.ucam as ucam
# ANSI Color Codes for terminal text formatting
BLUE = '\033[94m'
RED = '\033[91m'
RESET = '\033[0m'

class Reduction:
    def __init__(self, save_path):
        self.save_path = save_path
        self.dir_figs = os.path.join(save_path,'fig')
        os.makedirs(self.dir_figs, exist_ok= True)
        # Base directory for saving files
        os.chdir(self.save_path)
        
    def setaper(self, list_file, ccd_label='1', win_label='1', SIGMA_THRESHOLD=1.5,
                output_plot="detection_labeled.png", output_ape_name="ape.ape",
                MARGIN_LEFT=15, MARGIN_RIGHT=15, MARGIN_BOTTOM=5, MARGIN_TOP=27,
                SKIP_BRIGHTEST=5, R_TARG=None, R_SKY1=15, R_SKY2=20):
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
        target_file = hcm_files[0]
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
            # sources.reverse()
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
                    "rsky1": (R_TARG+5) * binning, 
                    "rsky2": (R_TARG+10) * binning,
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
                ax.add_patch(patches.Circle((x, y), R_TARG, edgecolor='lime', facecolor='none', lw=1.5))
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
    def modify_red(self, red_file="reduce.red", **kwargs):
        if not os.path.exists(red_file):
            print(f"{RED}Error: File '{red_file}' not found.{RESET}")
            return

        print(f"{BLUE}#### Modifying {red_file} ####{RESET}")
        
        with open(red_file, 'r') as f:
            lines = f.readlines()

        new_lines = []
        modified_keys = set()

        for line in lines:
            updated = False
           
            for key, value in kwargs.items():
                
                pattern = rf"^{key}\s*=\s*"
                if re.match(pattern, line.strip()):
                    comment = ""
                    if "#" in line:
                        comment = " # " + line.split("#", 1)[1].strip()
                    
                    new_lines.append(f"{key} = {value}{comment}\n")
                    modified_keys.add(key)
                    updated = True
                    break
            
            if not updated:
                new_lines.append(line)

        with open(red_file, 'w') as f:
            f.writelines(new_lines)

        for k in kwargs.keys():
            if k in modified_keys:
                print(f"   [Modified] {k} -> {kwargs[k]}")
            else:
                print(f"   [Skip] Key '{k}' not found in file (No changes made).")

        print(f" {red_file} updated successfully.")

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
                            fwhm = frame_data[f'fwhm_{ap}'].values[0]/self.binning
                            flag = frame_data[f'flag_{ap}'].values[0]
                       
                            color = 'lime' if flag == 0 else 'red'
                            plot_radius = fwhm if (pd.notna(fwhm) and fwhm > 0) else 5.0
                            
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


    # def read_hipercam_log(filename):
    # """read .log then return DataFrame """
    # if not os.path.exists(filename):
    #     return None, 0

    # base_cols = ['CCD', 'nframe', 'MJD', 'MJDok', 'Exptim', 'mfwhm', 'mbeta']
    # ap_cols_template = ['x', 'xe', 'y', 'ye', 'fwhm', 'fwhme', 'beta', 'betae', 
    #                     'counts', 'countse', 'sky', 'skye', 'nsky', 'nrej', 'cmax', 'flag']

    # max_aperture = 0
    # with open(filename, 'r') as f:
    #     for line in f:
    #         if line.startswith('#') and 'flag_' in line:
    #             parts = line.split()
    #             for part in parts:
    #                 if part.startswith('flag_'):
    #                     try:
    #                         ap = int(part.replace('flag_', ''))
    #                         if ap > max_aperture: max_aperture = ap
    #                     except: pass
    #         if not line.startswith('#'): break 

    # all_cols = base_cols.copy()
    # for i in range(1, max_aperture + 1):
    #     all_cols.extend([f"{col}_{i}" for col in ap_cols_template])

    # df = pd.read_csv(filename, comment='#', sep=r'\s+', header=None, names=all_cols)
    # return df, max_aperture
    