# import os
# import sys
# import importlib
# import concurrent.futures

# # =========================================================
# # 1. Receive Config File Path from Command Line Argument
# # =========================================================
# if len(sys.argv) < 2:
#     print("[Error] Missing config file argument.")
#     print("[Usage] python main.py <path_to_config.py>")
#     sys.exit(1)

# config_path = os.path.abspath(sys.argv[1])

# # Check if the file exists
# if not os.path.exists(config_path):
#     print(f"[Error] Config file not found at '{config_path}'")
#     sys.exit(1)

# # Separate directory and filename
# config_dir = os.path.dirname(config_path)
# config_filename = os.path.basename(config_path)
# module_name = os.path.splitext(config_filename)[0]  # Remove .py extension

# # Add path to system to allow importing
# if config_dir not in sys.path:
#     sys.path.insert(0, config_dir)

# # Dynamically load the Config file
# try:
#     config = importlib.import_module(module_name)
#     print(f"[Success] Loaded config from: {config_path}")
# except Exception as e:
#     print(f"[Error] Failed to load config file: {e}")
#     sys.exit(1)

# # =========================================================
# # 2. Import main modules
# # =========================================================
# from orderHIPERCAM import Reduction, Photometry
# from ensemble import Ensemble


# def process_single_run(run_config):
#     """
#     Function to process a single run (will be executed on 1 Core).
#     """
#     run_name = run_config.get('run_name', 'default_run')
#     print(f"\n[Core Task Started] RUN: {run_name}")

#     current_save_dir = os.path.join(config.BASE_SAVE_DIR, f"hcam_reduction_{run_name}")
#     ensemble_save_path = os.path.join(current_save_dir, 'results')
#     log_file_path = os.path.join(current_save_dir, 'data.log') 

#     try:
#         # Pipeline Reduction
#         pipeline = Reduction(
#             source=config.REDUCTION_SOURCE,
#             input_run_ul=run_config,
#             save_dir=current_save_dir,
#             raw_dir=getattr(config, 'DEFAULT_RAW_DIR', None),
#             fix_pixel=config.FIX_PIXEL,
#             diagnostics=config.DIAGNOSTICS
#         )

#         # Photometry
#         if not pipeline.data.get('lis'):
#             return f"[Skipped] {run_name}: No data list generated."

#         photo_reduction = Photometry(
#             base_dir=current_save_dir,
#             lis=pipeline.data['lis'],
#             diagnostics=True
#         )

#         photo_reduction.setaper(
#             ccd_label=config.CCD_LABEL, 
#             SKIP_BRIGHTEST=config.SKIP_BRIGHTEST, 
#             SIGMA_THRESHOLD=config.SIGMA_THRESHOLD,
#             frame=config.FRAME,
#             R_SKY1=config.R_SKY1, 
#             R_SKY2=config.R_SKY2,
#             MARGIN_LEFT=config.MARGIN_LEFT, 
#             MARGIN_RIGHT=config.MARGIN_RIGHT, 
#             MARGIN_BOTTOM=config.MARGIN_BOTTOM, 
#             MARGIN_TOP=config.MARGIN_TOP,
#             diagnostics=True,
#             ref_index=config.REF_INDEX
#         )

#         photo_reduction.genred()

#         photo_reduction.modify_red(
#             fit_method='moffat',
#             fit_height_min_ref=5.0,
#             fit_height_min_nrf=3.0,
#             fit_half_width=25.0,
#             search_half_width=15.0,
#         )
        
#         extract_vals = " ".join(map(str, config.R_EXTRACT))
#         photo_reduction.modify_red(
#             target_section="extraction",
#             **{"1": f"variable normal {extract_vals}"}
#         )

#         photo_reduction.reduce(
#             R_EXTARCT=config.R_EXTRACT,
#             plot_with_zoom=False,
#             plot_with_log=True
#         )

#         # Ensemble Calibration
#         actual_log_files = pipeline.data.get('log_file', [])
#         log_to_use = actual_log_files[0] if actual_log_files else log_file_path

#         if os.path.exists(log_to_use):
#             ens = Ensemble(
#                 file=log_to_use,
#                 save_path=ensemble_save_path,
#                 solvwcs=True,
#             )

#             ens.run(
#                 target_rms=config.TARGET_RMS,
#                 numstars=config.NUM_STARS,
#                 ra_center=config.RA_CENTER,
#                 dec_center=config.DEC_CENTER,
#                 radius=config.RADIUS,
#                 scale_low=config.SCALE_LOW,
#                 scale_high=config.SCALE_HIGH,
#                 single_frame=config.SINGLE_FRAME,
#                 astrometry_cache=config.ASTROMETRY_CACHE,
#             )
#             return f"[Success] {run_name} completed."
#         else:
#             return f"[Warning] Cannot find log file for {run_name}. Skipping Ensemble."
            
#     except Exception as e:
#         return f"[Failed] {run_name} encountered an error: {e}"


# def main():
#     print('xx' * 50)
#     print(f"--- Initializing Multiprocessing Pipeline ---")
    
#     num_tasks = len(config.INPUT_RUNS)
#     max_cores = min(num_tasks, os.cpu_count() or 4)
#     print(f"Total Tasks: {num_tasks} | Using Cores: {max_cores}")
#     print("=" * 60)

#     with concurrent.futures.ProcessPoolExecutor(max_workers=max_cores) as executor:
#         results = executor.map(process_single_run, config.INPUT_RUNS)

#     print("\n" + "=" * 60)
#     print("BATCH PIPELINE SUMMARY:")
#     for res in results:
#         print(res)
#     print("=" * 60)

# if __name__ == "__main__":
#     main()



import os
import sys
import importlib
import concurrent.futures

# =========================================================
# 1. Receive Config File Path from Command Line Argument
# =========================================================
if len(sys.argv) < 2:
    print("[Error] Missing config file argument.")
    print("[Usage] python main.py <path_to_config.py>")
    sys.exit(1)

config_path = os.path.abspath(sys.argv[1])

# Check if the file exists
if not os.path.exists(config_path):
    print(f"[Error] Config file not found at '{config_path}'")
    sys.exit(1)

# Separate directory and filename
config_dir = os.path.dirname(config_path)
config_filename = os.path.basename(config_path)
module_name = os.path.splitext(config_filename)[0]  # Remove .py extension

# Add path to system to allow importing
if config_dir not in sys.path:
    sys.path.insert(0, config_dir)

# Dynamically load the Config file
try:
    config = importlib.import_module(module_name)
    print(f"[Success] Loaded config from: {config_path}")
except Exception as e:
    print(f"[Error] Failed to load config file: {e}")
    sys.exit(1)

# =========================================================
# 2. Import main modules
# =========================================================
from orderHIPERCAM import Reduction, Photometry
from ensemble import Ensemble


def process_single_run(run_config):
    """
    Function to process a single run (will be executed on 1 Core).
    """
    run_name = run_config.get('run_name', 'default_run')

    # ---------------------------------------------------------
    # Custom Logger to prefix ALL prints from this specific core
    # ---------------------------------------------------------
    class ProcessLogger:
        def __init__(self, name):
            self.name = name
            self.terminal = sys.__stdout__

        def write(self, message):
            if message.strip():  # Ignore empty newlines
                for line in message.splitlines():
                    if line.strip():
                        self.terminal.write(f"[{self.name}] {line}\n")
                        self.terminal.flush()

        def flush(self):
            self.terminal.flush()

    # Override standard output and error for this specific process
    sys.stdout = ProcessLogger(run_name)
    sys.stderr = ProcessLogger(f"{run_name} | ERROR")

    print(f"Core Task Started")

    current_save_dir = os.path.join(config.BASE_SAVE_DIR, f"hcam_reduction_{run_name}")
    ensemble_save_path = os.path.join(current_save_dir, 'results')
    log_file_path = os.path.join(current_save_dir, 'data.log') 

    try:
        # Pipeline Reduction
        pipeline = Reduction(
            source=config.REDUCTION_SOURCE,
            input_run_ul=run_config,
            save_dir=current_save_dir,
            raw_dir=getattr(config, 'DEFAULT_RAW_DIR', None),
            fix_pixel=config.FIX_PIXEL,
            diagnostics=config.DIAGNOSTICS
        )

        # Photometry
        if not pipeline.data.get('lis'):
            return f"[{run_name}] Skipped: No data list generated."

        photo_reduction = Photometry(
            base_dir=current_save_dir,
            lis=pipeline.data['lis'],
            diagnostics=True
        )

        photo_reduction.setaper(
            ccd_label=config.CCD_LABEL, 
            SKIP_BRIGHTEST=config.SKIP_BRIGHTEST, 
            SIGMA_THRESHOLD=config.SIGMA_THRESHOLD,
            frame=config.FRAME,
            R_SKY1=config.R_SKY1, 
            R_SKY2=config.R_SKY2,
            MARGIN_LEFT=config.MARGIN_LEFT, 
            MARGIN_RIGHT=config.MARGIN_RIGHT, 
            MARGIN_BOTTOM=config.MARGIN_BOTTOM, 
            MARGIN_TOP=config.MARGIN_TOP,
            diagnostics=True,
            ref_index=config.REF_INDEX
        )

        photo_reduction.genred()

        photo_reduction.modify_red(
            fit_method='moffat',
            fit_height_min_ref=5.0,
            fit_height_min_nrf=3.0,
            fit_half_width=25.0,
            search_half_width=15.0,
        )
        
        extract_vals = " ".join(map(str, config.R_EXTRACT))
        photo_reduction.modify_red(
            target_section="extraction",
            **{"1": f"variable normal {extract_vals}"}
        )

        photo_reduction.reduce(
            R_EXTARCT=config.R_EXTRACT,
            plot_with_zoom=False,
            plot_with_log=True
        )

        # Ensemble Calibration
        actual_log_files = pipeline.data.get('log_file', [])
        log_to_use = actual_log_files[0] if actual_log_files else log_file_path

        if os.path.exists(log_to_use):
            ens = Ensemble(
                file=log_to_use,
                save_path=ensemble_save_path,
                solvwcs=True,
            )

            ens.run(
                target_rms=config.TARGET_RMS,
                numstars=config.NUM_STARS,
                ra_center=config.RA_CENTER,
                dec_center=config.DEC_CENTER,
                radius=config.RADIUS,
                scale_low=config.SCALE_LOW,
                scale_high=config.SCALE_HIGH,
                single_frame=config.SINGLE_FRAME,
                astrometry_cache=config.ASTROMETRY_CACHE,
                ignor_stars=config.IGNOR_STARS
            )
            return f"[{run_name}] Success: Run completed."
        else:
            return f"[{run_name}] Warning: Cannot find log file. Skipping Ensemble."
            
    except Exception as e:
        return f"[{run_name}] Failed: Encountered an error: {e}"


def main():
    print('x' * 90)
    print(f"--- Initializing Multiprocessing Pipeline ---")
    
    num_tasks = len(config.INPUT_RUNS)
    max_cores = min(num_tasks, os.cpu_count() or 4)
    print(f"Total Tasks: {num_tasks} | Using Cores: {max_cores}")
    print("=" * 90)

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_cores) as executor:
        results = executor.map(process_single_run, config.INPUT_RUNS)

    print("\n" + "=" * 90)
    print("BATCH PIPELINE SUMMARY:")
    for res in results:
        print(res)
    print("=" * 90)

if __name__ == "__main__":
    main()