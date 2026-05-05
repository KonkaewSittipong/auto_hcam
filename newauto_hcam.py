import os
import math
import shutil
import glob
import json
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from datetime import datetime
from typing import Optional
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.visualization import ZScaleInterval
from photutils.detection import DAOStarFinder

import hipercam as hcam
from hipercam import MCCD
import hipercam.scripts as scripts

"""
auto_hcam

"""
import os
import math
import shutil
import glob
import json
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from datetime import datetime
from typing import Optional
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.visualization import ZScaleInterval
from photutils.detection import DAOStarFinder

import hipercam as hcam
from hipercam import MCCD
import hipercam.scripts as scripts

BLUE = '\033[94m'
RED = '\033[91m'
RESET = '\033[0m'

class hipercam_setup:
    """
    Parameters:
       source : str
           Data source, five options:

              | 'hs' : HiPERCAM server
              | 'hl' : local HiPERCAM FITS file
              | 'us' : ULTRACAM server
              | 'ul' : local ULTRACAM .xml/.dat files
              | 'hf' : list of HiPERCAM hcm FITS-format files

    """
    def __init__(self, source, input_run_ul=None, ccd='1', save_dir=None, raw_dir=None, diagnostics=False):
        # Calibration tracking ###file###
        self.bias = {'file': [], 'master': [], 'lis': []}
        self.dark = {'file': [], 'master': [], 'lis': []}
        self.flat = {'file': [], 'master': [], 'lis': []}
        self.data = {'file': [], 'lis': []}
        
        # Path Setup           #Collect path of ###dir###
        self.bias_dir = None
        self.dark_dir = None
        self.flat_dir = None
        
        self.raw_dir = raw_dir
        self.base_dir = save_dir
        self.source = source
        
        if self.base_dir : ######### 
            os.makedirs(self.base_dir, exist_ok=True)
        
        
        # Initialize logic based on source
        if source == 'ul' :            
            if input_run_ul is not None:
                print('Work with source : ul')
                self.input_run_ul = input_run_ul
                self.source_ul()
                
            else:
                print('Source : ul must to input --input_run_ul-- ')
                return
        self.reduction()
        
        self.diagnostics = diagnostics
        if self.diagnostics:
            self.figs_dir = os.path.join(self.base_dir, 'figs')
            os.makedirs(self.figs_dir, exist_ok=True)
            self.diagnostic()
            
    def diagnostic(self):
        """Runs diagnostics for Bias, Dark, and Flat frames (Both Raw and Master)."""
        print(f"{BLUE}---  Running Diagnostics: Calibration Analysis ---{RESET}")

        categories = [
            ('bias', self.bias, "Bias"),
            ('dark', self.dark, "Dark"),
            ('flat', self.flat, "Flat")
        ]

        for key, data_dict, label in categories:
            # --- A. Plot Raw Frames & Trend Lines ---
            if data_dict['file']:
                self.plot_run(
                    data_dict['file'], 
                    save_name=f"diagnostic_raw_{key}.png", 
                )

            # --- B. Plot Master Frame ---
            if data_dict['master']:
                print(f" Plotting Master {label}...")
                self.plot_run(
                    data_dict['master'], 
                    save_name=f"diagnostic_master_{key}.png", 
                    show_median=True
                )

        print(f"{BLUE}--- Diagnostics Complete! Files saved in: {self.figs_dir} ---{RESET}\n")

    def plot_run(self, files, save_name, target_ccd='1', df=None, aps_count=0, show_median=False, ):
        """
        Plots a grid of HCM images and generates a median level trend plot.
        
        Parameters:
            files (list): List of full paths to .hcm files.
            save_name (str): Name for the output grid image.
            target_ccd (str): CCD label to plot (e.g., '1', '2', '3').
            df (DataFrame): Optional pandas DF containing aperture coordinates.
            aps_count (int): Number of apertures to plot from the DF.
            show_median (bool): If True, calculates and displays median ADU in titles.
        """
        if not files:
            return

        # --- Setup Data Collection for Trend Plot ---
        frame_numbers = []
        median_levels = []
        
        # --- Grid Calculations ---
        ncols = 5
        nrows = math.ceil(len(files) / ncols)
        if nrows == 0: nrows = 1 
        
        fig_grid, axs = plt.subplots(nrows, ncols, figsize=(6*ncols, 6*nrows), 
                                     gridspec_kw={'hspace': 0.25, 'wspace': 0.05}, facecolor='white')
        
        # Bulletproof flattening (works for 1 plot, 1D array, or 2D array)
        if hasattr(axs, 'flatten'):
            axs_flat = axs.flatten()
        else:
            axs_flat = [axs]
        
        zscale = ZScaleInterval()

        for i, fname in enumerate(files):
            ax = axs_flat[i]
            try:
                mccd = MCCD.read(fname)
                data = mccd[str(target_ccd)]['1'].data 
                
                # --- NEW ROBUST STATISTICS ---
                # 1. Clip the junk (hides dead pixels and hot pixels beyond 3 sigma)
                clipped_data = sigma_clip(data, sigma=3.0)
                
                # 2. Calculate stats only on the good pixels
                median_val = np.ma.median(clipped_data)
                min_val = np.ma.min(clipped_data)
                max_val = np.ma.max(clipped_data)
                
                frame_numbers.append(i + 1)
                median_levels.append(median_val)
                
                # Display Image
                vmin, vmax = zscale.get_limits(data)
                ax.imshow(data, origin='lower', cmap='gray_r', vmin=vmin, vmax=vmax)
                
                # Title Logic
                title_str = f"[#{i+1}] {os.path.basename(fname)}"
                if show_median:
                    title_str += f"\nMed: {median_val:.5f} | Min: {min_val:.5f} | Max: {max_val:.5f}"
                ax.set_title(title_str, fontsize=8)
                ax.axis('off')
                
                
                # --- Aperture Overlay (For Science Mode Later) ---
                if df is not None and i < len(df):
                    row = df.iloc[[i]]
                    if not row.empty:
                        for ap in range(1, aps_count + 1):
                            col_x, col_y = f'x_{ap}', f'y_{ap}'
                            if col_x in row.columns and col_y in row.columns:
                                x_, y_ = row.iloc[0][col_x], row.iloc[0][col_y]
                                if np.isfinite(x_) and np.isfinite(y_):
                                    targ = patches.Circle((x_-0.5, y_-0.5), 10, 
                                                          edgecolor='lime', facecolor='none', lw=1.0)
                                    ax.add_patch(targ)
                                    ax.text(x_, y_ + 15, f"{ap}", color='red', fontsize=6, ha='center')

            except Exception as e:
                ax.set_title(f"Read Error", color='red', fontsize=9)
                ax.axis('off')
                print(f"{RED}Error plotting {fname}: {e}{RESET}")

        # Hide empty subplot slots so the final image is clean
        for j in range(i + 1, len(axs_flat)):
            axs_flat[j].axis('off')

        # Save the Image Grid
        grid_path = os.path.join(self.figs_dir, save_name)
        plt.savefig(grid_path, bbox_inches='tight', dpi=100)
        plt.close(fig_grid)

        # --- Generate Trend Plot (Only if there are multiple frames) ---
        if len(frame_numbers) > 1:
            plt.figure(figsize=(12, 5))
            plt.plot(frame_numbers, median_levels, 'o-', markersize=4, color='tab:blue', label='Median ADU')
            
            plt.title(f"Median Level Trend: {os.path.basename(save_name).replace('.png', '')}")
            plt.xlabel("Frame Index")
            plt.ylabel("Counts (ADU)")
            plt.grid(True, alpha=0.3)
            plt.legend()
            
            trend_path = os.path.join(self.figs_dir, save_name.replace('.png', '_trend.png'))
            plt.savefig(trend_path, bbox_inches='tight')
            plt.close()
            print(f"   -> Saved {save_name} and trend plot.")
        else:
            print(f"   -> Saved {save_name}.")
        
    def grab(self, work_dir, run_path, prefix_out,  **kwargs):
        """Grab raw frames from a single run."""
        if not os.path.exists(work_dir):
            os.makedirs(work_dir, exist_ok=True)
            
        os.chdir(work_dir)
        run = os.path.basename(run_path)
        output_name = f"{prefix_out}_{run}"
        args = ["nodef", 
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
                # f"output={output_name}"
               ]
        try:
            print(f"{BLUE}#### Grabbing {run} from {run_path} ####{RESET}")
            scripts.grab(args)
            # raw_files = glob.glob(f"{run}_*.hcm")
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

    def makebias(self, lis_filename, output_name, sigma=3,):
        """Creates a Master Bias."""      
        args = [
            "nodef", 
            "source=hf", 
            f"flist={lis_filename}", 
            f"sigma={sigma}", 
            "plot=no", 
            f"output={output_name}"  
        ]
        try:
            print(f"{BLUE}--- Making Master Bias: {os.path.basename(output_name)} ---{RESET}")
            scripts.makebias(args)
            
            full_master_path = os.path.abspath(f'{output_name}.hcm')
            self.bias['master'].append(full_master_path)
            
        except Exception as e:
            print(f"{RED}Error processing bias: {e}{RESET}")


            
    def makeflat(self, lis_filename, output_name, sigma=3, chosen_bias = 'none', chosen_dark='none',):
        """Creates a Master Flat."""      
        args = [
            "nodef", 
            "source=hf", 
            "ngroup=5", 
            f"flist={lis_filename}", 
            f"sigma={sigma}",
            "lower=10000", 
            "upper=35000", 
            "plot=no", 
            f"bias={chosen_bias}", 
            # f"flat={chosen_flat}", 
            f"dark={chosen_dark}", 
            f"output={output_name}"  
        ]
        try:
            full_master_path = os.path.abspath(f'{output_name}.hcm')
            if os.path.exists(full_master_path):
                os.remove(full_master_path)
                print(f"Removed existing file: {os.path.basename(full_master_path)}")
                
            print(f"{BLUE}--- Making Master Flat: {os.path.basename(output_name)} ---{RESET}")
            scripts.makeflat(args)
            self.flat['master'].append(full_master_path)
            
        except Exception as e:
            print(f"{RED}Error processing bias: {e}{RESET}")

    def source_ul(self,):
        """Grab bias to .hcm"""

        self.bias_dir = os.path.join(self.base_dir, 'bias') 
        os.makedirs(self.base_dir, exist_ok=True) # Create the specific bias subdirectory
        
        # Filter keys starting with bias_
        bias_subset = {k: v for k, v in self.input_run_ul.items() if k.startswith('bias_')}
        for key, value in bias_subset.items():
            if value is not None:
                rawdir = value.get('rawdir')
                if not rawdir:
                    rawdir = self.raw_dir
                    
                for run_info in value.get('runs', []):
                    run, f1, f2 = run_info
                    run_path = os.path.join(rawdir, run)
                    self.grab(
                        work_dir=self.bias_dir,
                        run_path=run_path,
                        prefix_out=f"{key}",
                        f1=f1,
                        f2=f2
                    )
                relevant_files_path = [f for f in self.bias['file'] if f"{key}_" in os.path.basename(f)]
                if relevant_files_path:
                    lis_filename = os.path.join(self.bias_dir, f"{key}.lis")
                    with open(lis_filename, 'w') as f:
                        for filepath in relevant_files_path:
                            f.write(f"{filepath}\n") 
                    full_lis_path = os.path.abspath(lis_filename)
                    self.bias['lis'].append(full_lis_path)

                    master_out_path = os.path.join(self.bias_dir, f"master_{key}")
                    self.makebias(lis_filename=full_lis_path, output_name=master_out_path)
                print(f"Processing {key} finished...\n")

        """Grab flat to .hcm"""
        self.flat_dir = os.path.join(self.base_dir, 'flat') 
        os.makedirs(self.flat_dir, exist_ok=True) # Create the specific flat subdirectory
        flat_subset = {k: v for k, v in self.input_run_ul.items() if k.startswith('flat_')}
        for key, value in flat_subset.items():
            # find mater bias shoude in self.bias['master'].value() if bias with flat if not use bias with data
            if value is not None:
                #################
                chosen_bias = 'none'
                for b_path in self.bias['master']:
                    if f'master_{key}' in os.path.basename(b_path):
                        chosen_bias = b_path
                        break
                if chosen_bias == 'none':
                    for b_path in self.bias['master']:
                        if 'master_bias_all' in os.path.basename(b_path) or 'master_bias_data' in os.path.basename(b_path):
                            chosen_bias = b_path
                            break
                #################

                rawdir = value.get('rawdir')
                if not rawdir:
                    rawdir = self.raw_dir
                
                for run_info in value.get('runs', []):
                    run, f1, f2 = run_info
                    print("x=",run, f1, f2)
                    run_path = os.path.join(rawdir, run)
                    self.grab(
                        work_dir=self.flat_dir,
                        run_path=run_path,
                        prefix_out=f"{key}",
                        f1=f1,
                        f2=f2,
                        bias=chosen_bias
                    )
                relevant_files_path = [f for f in self.flat['file'] if f"{key}_" in os.path.basename(f)]
                if relevant_files_path:
                    lis_filename = os.path.join(self.flat_dir, f"{key}.lis")
                    with open(lis_filename, 'w') as f:
                        for filepath in relevant_files_path:
                            f.write(f"{filepath}\n") 
                    full_lis_path = os.path.abspath(lis_filename)
                    self.flat['lis'].append(full_lis_path)

                    master_out_path = os.path.join(self.flat_dir, f"master_{key}")
                    self.makeflat(lis_filename=full_lis_path, output_name=master_out_path)
                    
                # print(f"Processing {key} finished...\n")
 
    def reduction(self,):
        if self.source == 'ul':
            """Grab Data with reduce noise to .hcm"""
    
            self.data_dir = os.path.join(self.base_dir, 'data') 
            os.makedirs(self.data_dir, exist_ok=True) # Create the specific bias subdirectory
            
            data_subset = {k: v for k, v in self.input_run_ul.items() if k.startswith('data')}
            for key, value in data_subset.items():
                if value is not None:
                    rawdir = value.get('rawdir')
                    if not rawdir:
                        rawdir = self.raw_dir
                    for run_info in value.get('runs', []):
                        run, f1, f2 = run_info
                        run_path = os.path.join(rawdir, run)
                        
                        ##############
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
                        ###############
                        self.grab(
                            work_dir=self.data_dir,
                            run_path=run_path,
                            prefix_out=f"{key}",
                            f1=f1,
                            f2=f2,
                            bias=chosen_bias,
                            flat=chosen_flat)
                    relevant_files_path = [f for f in self.data['file'] if f"{key}_" in os.path.basename(f)]
                    if relevant_files_path:
                        lis_filename = os.path.join(self.data_dir, f"{key}.lis")
                        with open(lis_filename, 'w') as f:
                            for filepath in relevant_files_path:
                                f.write(f"{filepath}\n")

                        full_lis_path = os.path.abspath(lis_filename)
                        self.data['lis'].append(full_lis_path)


class photometry:
    def __init__(self, base_dir, lis, diagnostics=False, single_run = True):
        self.diagnostics = diagnostics
        self.base_dir = base_dir
        self.lis = lis
        self.single_run = single_run
        self.data = []
    def solvwcs():
        print('x')
        

    def setaper(self, #files_lis, 
                ccd_label='1', win_label='1', SIGMA_THRESHOLD=1.5,
                output_plot="detection_labeled.png", SKIP_BRIGHTEST=5,
                MARGIN_LEFT=15, MARGIN_RIGHT=15, MARGIN_BOTTOM=5, MARGIN_TOP=27,
                 R_TARG=None, R_SKY1=16, R_SKY2=24, frame=5, diagnostics=False):
        
        """
        Detects stars from the image in self.lis and creates 
        a HiPERCAM aperture (.ape) file.
        """

        
        hcm_files = []
        for lis_path in self.lis:
            if os.path.exists(lis_path):
                with open(lis_path, 'r') as f:
                    hcm_files.extend([line.strip() for line in f if line.strip()]) 
            
        if self.single_run:
            target_file = hcm_files[1]
            out_ape = os.path.splitext(target_file)[0] + '.ape'
            if not target_file or not os.path.exists(target_file):

                print(f"{RED}Error: List file '{target_file}' not found.{RESET}")
                return
    
            try:
                mccd = hcam.MCCD.read(target_file)
                window = mccd[ccd_label][win_label]
                data = data = window.data
                binning = getattr(window, 'xbin', 1)
                self.binning = binning 
                h, w = data.shape
            except Exception as e:
                print(f"{RED}Critical Error reading HCM file: {e}{RESET}")
                return
    
            
            _, median, std = sigma_clipped_stats(data, sigma=3)
            daofind = DAOStarFinder(fwhm=4.0, threshold=SIGMA_THRESHOLD * std)
            sources = daofind(data - median)
            if sources:
                sources = sources[::-1]
    
                # Filter stars based on defined margins
                mask = ((sources['xcentroid'] > MARGIN_LEFT) & 
                        (sources['xcentroid'] < w - MARGIN_RIGHT) & 
                        (sources['ycentroid'] > MARGIN_BOTTOM) & 
                        (sources['ycentroid'] < h - MARGIN_TOP))
                sources = sources[mask]
    
                sources = sources[SKIP_BRIGHTEST:]
                sources['fwhm'] = np.sqrt(sources['npix'] / np.pi) * 2.35
                # print(f"Detecting sources on {os.path.basename(target_file)} with threshold {SIGMA_THRESHOLD} sigma. Skipping the first {SKIP_BRIGHTEST} brightest stars. Total usable stars detected: {len(sources)}")
                self.num_aps = len(sources)  
             
                ccd_ape = []
                for i, source in enumerate(sources):
                    x, y, fwhm = source['xcentroid'], source['ycentroid'], source['fwhm']
                    ap_id = str(i + 1)
                    is_ref = (i == 0) # The first star in the remaining list becomes the reference
                    # if R_TARG  is None : R_TARG = fwhm * 1
                    r_t = (R_TARG if R_TARG is not None else fwhm) * binning
                    
                    ccd_ape.append([ap_id, {
                        "Comment": "hipercam.Aperture", 
                        "x": float(x) * binning, 
                        "y": float(y) * binning,
                        "rtarg": r_t, 
                        "rsky1": (r_t+4) * binning, 
                        "rsky2": (r_t+8) * binning,
                        "ref": is_ref, 
                        "compo": False, 
                        "mask": [], 
                        "extra": [], 
                        "link": "" if i == 0 else "1",
                    }])
                    
                ape_json = ["hipercam.MccdAper",
                    [ccd_label, ["hipercam.CcdAper"] + ccd_ape]]
                with open(out_ape, 'w') as f:
                    json.dump(ape_json, f, indent=2)
                # print(f"Aperture file: {out_ape}")

                self.data.append({
                    'target': target_file,
                    'out_ape': out_ape, 
                    'sources':len(sources)
                })
             
        else:
            for j, target_file in enumerate(hcm_files):
                
                out_ape = os.path.splitext(target_file)[0] + '.ape'
                if not target_file or not os.path.exists(target_file):
    
                    print(f"{RED}Error: List file '{target_file}' not found.{RESET}")
                    return
        
                try:
                    mccd = hcam.MCCD.read(target_file)
                    window = mccd[ccd_label][win_label]
                    data = data = window.data
                    binning = getattr(window, 'xbin', 1)
                    self.binning = binning 
                    h, w = data.shape
                except Exception as e:
                    print(f"{RED}Critical Error reading HCM file: {e}{RESET}")
                    return
        
                
                _, median, std = sigma_clipped_stats(data, sigma=3)
                daofind = DAOStarFinder(fwhm=4.0, threshold=SIGMA_THRESHOLD * std)
                sources = daofind(data - median)
                if sources:
                    sources = sources[::-1]
        
                    # Filter stars based on defined margins
                    mask = ((sources['xcentroid'] > MARGIN_LEFT) & 
                            (sources['xcentroid'] < w - MARGIN_RIGHT) & 
                            (sources['ycentroid'] > MARGIN_BOTTOM) & 
                            (sources['ycentroid'] < h - MARGIN_TOP))
                    sources = sources[mask]
        
                    sources = sources[SKIP_BRIGHTEST:]
                    sources['fwhm'] = np.sqrt(sources['npix'] / np.pi) * 2.35
                    # print(f"Detecting sources on {os.path.basename(target_file)} with threshold {SIGMA_THRESHOLD} sigma. Skipping the first {SKIP_BRIGHTEST} brightest stars. Total usable stars detected: {len(sources)}")
                    self.num_aps = len(sources)  
                 
                    ccd_ape = []
                    for i, source in enumerate(sources):
                        x, y, fwhm = source['xcentroid'], source['ycentroid'], source['fwhm']
                        ap_id = str(i + 1)
                        is_ref = (i == 0) # The first star in the remaining list becomes the reference
                        # if R_TARG  is None : R_TARG = fwhm * 1
                        r_t = (R_TARG if R_TARG is not None else fwhm) * binning
                        
                        ccd_ape.append([ap_id, {
                            "Comment": "hipercam.Aperture", 
                            "x": float(x) * binning, 
                            "y": float(y) * binning,
                            "rtarg": r_t, 
                            "rsky1": (r_t+4) * binning, 
                            "rsky2": (r_t+8) * binning,
                            "ref": is_ref, 
                            "compo": False, 
                            "mask": [], 
                            "extra": [], 
                            "link": "" if i == 0 else "1",
                        }])
                        
                    ape_json = ["hipercam.MccdAper",
                        [ccd_label, ["hipercam.CcdAper"] + ccd_ape]]
                    with open(out_ape, 'w') as f:
                        json.dump(ape_json, f, indent=2)
                    # print(f"Aperture file: {out_ape}")
    
                    self.data.append({
                        'target': target_file,
                        'out_ape': out_ape, 
                        'sources':len(sources)
                    })
    def genred(self, apertures=None, instru='ultraspec-tnt', **kwargs):
        """
        Runs the HiPERCAM 'genred' command to create a reduction control file (.red).
        """
        # 1. Figure out which entries to process and save it to the class (self)
        if self.single_run:
            self.entries_to_process = [self.data[0]] 
        else:
            self.entries_to_process = self.data
            
        # 2. Process them in a loop
        for entry in self.entries_to_process:
            print(entry)
            ape_file = entry.get('out_ape')
            print(ape_file)
            if not ape_file: 
                continue
            
            # Decide the output name
            if self.single_run:
                # Grab the folder path, and stick 'all.red' at the end
                out_red = os.path.join(os.path.dirname(ape_file), 'all.red')
            else:
                # Swap the .ape extension for .red
                out_red = os.path.splitext(ape_file)[0] + '.red'
                
            args = ["genred", ape_file, out_red, 
                    "none", "none", "none", "none", 
                    "0", "none", f"{instru}"]
            
            try:
                scripts.genred(args)
                
                # This appends the FULL PATH of the .red file to your entry dictionary
                entry['out_red'] = out_red 
                
                print(f"Successfully created: {os.path.basename(out_red)}")
            except Exception as e:
                print(f"{RED}Error running genred on {ape_file}: {e}{RESET}")
                
                

    def modify_red(self, target_section=None, **kwargs):
        for entry in self.entries_to_process:
            red_file = entry.get('out_red')
            if not red_file or not os.path.exists(red_file):
                continue
    
            print(f"{BLUE}#### Modifying {red_file} [{target_section if target_section else 'All'}] ####{RESET}")
            
            with open(red_file, 'r') as f:
                lines = f.readlines()
    
            new_lines = []
            modified_keys = set()
            current_section = None
    
            for line in lines:
                strip_line = line.strip()
                
                if strip_line.startswith('[') and strip_line.endswith(']'):
                    current_section = strip_line[1:-1].strip().lower()
    
                updated = False
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

    def reduce(self, sep_mode=False, diagnostics=False):
        if self.single_run:
            log_file = "all.log"
            
            # 1. Safely change to your base directory
            os.chdir(self.base_dir)
            
            # 2. Safely grab the .lis file path
            list_file = self.lis[0] 
            print(list_file)
            
            # 3. Safely grab the .red file from the bridge we built earlier
            red_file = self.entries_to_process[0].get('out_red')
            
            reduce_args = ["reduce", "source=hf", f"flist={list_file}", "trim=False", 
                           f"rfile={red_file}", f"log={log_file}", "tkeep=1", 
                           "lplot=False", "implot=False"]
            try:
                print(f"{BLUE} Starting FULLY AUTOMATED Reduction...{RESET}")
                scripts.reduce(reduce_args)
                print(f" SUCCESS! Log saved to: {log_file}")
            except Exception as e:
                print(f"{RED} Reduction Failed: {e}{RESET}")

        else:
            print("Not def yet")

            
                
"""
input_run = {
    'bias_data': {
        'runs': [['run002', 10, 0]], 
        'rawdir': '/lustre/MSSP/sittipong/all/data/2023_11_07',
        'ccd': '1'
    },
    'bias_dark': None, 
    'bias_flat': None, # {
        # 'runs': [['run002', 2, 0]], 
        # 'ccd': '1', 
    #     # 'rawdir': '/lustre/MSSP/sittipong/all/data/2023_11_07'
    # },
    'dark_data': None,
    'dark_flat': None,
    'flat_data': {
        'runs': [['run009', 2, 57]], 
        'ccd': '1', 
        'rawdir': '/lustre/MSSP/sittipong/all/data/2023_11_07'
    },
    'data': {
        'runs': [
            ['run010', 1, 30], ['run014', 1, 20], ['run018', 1, 20],
            ['run022', 1, 20], ['run026', 1, 20], ['run030', 1, 20],
            ['run037', 1, 20], ['run041', 1, 15]
        ],
        'ccd': '3', 
        'rawdir': '/lustre/MSSP/sittipong/all/data/2023_11_07'
    }
}

# ==========================================
# Jupyter Notebook
# ==========================================

save_directory = os.path.join('/lustre/MSSP/sittipong/temp/', 'hcam_reduction')

default_raw_directory = '/lustre/MSSP/sittipong/temp'

print(f"{BLUE}---  Initializing Pipeline ---{RESET}")
pipeline = hipercam_setup(
    source='ul', 
    input_run_ul=input_run, 
    ccd='1', 
    save_dir=save_directory, 
    raw_dir=default_raw_directory,
    diagnostics=True
)

# print(pipeline.bias['lis'])
# print(pipeline.dark['lis'])
# print(pipeline.flat['lis'])


photo_reduction = photometry(base_dir=save_directory, lis=pipeline.data['lis'])
print('impot photometry done')
# print(f"\n{BLUE}---  End of Process ---{RESET}")

R_SKY1=12 
R_SKY2=15        
photo_reduction.setaper(
    ccd_label='1', 
    SKIP_BRIGHTEST=10,
    # list_file=science_lis,  # <--- Pass the .lis file here
    SIGMA_THRESHOLD=7,
    frame=5,
    R_SKY1=R_SKY1,
    R_SKY2=R_SKY2
    
)
# print(photo_reduction.data)
photo_reduction.genred()

photo_reduction.modify_red(
    fit_method = 'moffat',
    fit_height_min_ref = 5.0,
    fit_height_min_nrf = 3.0,
    fit_half_width = 25.0,
    search_half_width = 15.0,
    )

photo_reduction.modify_red(
    target_section="extraction",
    **{"1": f"variable normal 1.80 3.0 10.0 2.5 {R_SKY1} {R_SKY1} 3.0 {R_SKY2} {R_SKY2}"}
)

photo_reduction.modify_red(fit_max_shift = 8, search_smooth_fft= 'yes')
photo_reduction.reduce()

"""
