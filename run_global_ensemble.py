"""
Run multi-epoch ensemble calibration.

Usage:
    python run_multi_epoch.py
"""

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from global_ensemble import MultiEpochEnsemble, find_epoch_dirs

# =====================================================================
# 1. Collect run directories
#    Each path must contain calibrated_lightcurves.txt with RA_N/DEC_N
# =====================================================================

# Option A — list them manually
epoch_dirs = [
    '/lustre/MSSP/sittipong/reduce2/2018_11_07/hcam_reduction_2018_11_07_r/results/run001',
    '/lustre/MSSP/sittipong/reduce2/2018_11_08/hcam_reduction_2018_11_08_r/results/run001',
    '/lustre/MSSP/sittipong/reduce2/2023-11-06/hcam_reduction_2023_11_06_r/results/run001',
    '/lustre/MSSP/sittipong/reduce2/2023-11-07/hcam_reduction_2023_11_07_r/results/run001',
    '/lustre/MSSP/sittipong/reduce2/2023-11-08/hcam_reduction_2023_11_08_r/results/run001',
    '/lustre/MSSP/sittipong/reduce2/2023-11-09/hcam_reduction_2023_11_09_r/results/run001',
    '/lustre/MSSP/sittipong/reduce2/2024-05-12/hcam_reduction_2024_05_12_r/results/run001',

]

# Option B — auto-discover by glob pattern (uncomment to use)
# epoch_dirs = find_epoch_dirs(
#     '/lustre/MSSP/sittipong/reduce2/*/hcam_reduction_*_i/results/run*'
# )

# =====================================================================
# 2. Output directory
# =====================================================================
OUTPUT_DIR = '/lustre/MSSP/sittipong/multi_epoch/r_band'

# =====================================================================
# 3. Target & field parameters
# =====================================================================
TARGET_RA  = 322.4377   # deg  (science target)
TARGET_DEC = -4.4853    # deg

TOLERANCE_ARCSEC = 10  # sky-matching radius for cross-epoch star ID

# =====================================================================
# 4. Ensemble parameters
# =====================================================================
TARGET_RMS = 0.02
NUM_STARS  = 12         # minimum reference stars to keep

# =====================================================================
# 5. SDSS calibration  (set to None to skip)
# =====================================================================
SDSS_FILTER              = 'r'   # 'g', 'r', 'i', 'z'  or None
SDSS_SEARCH_RADIUS_ARCSEC = 10
SDSS_DATA_RELEASE         = 18

# =====================================================================
# Run
# =====================================================================
if __name__ == '__main__':

    print(f"Epoch directories found: {len(epoch_dirs)}")
    for d in epoch_dirs:
        print(f"  {d}")
    print()

    me = MultiEpochEnsemble(
        epoch_dirs       = epoch_dirs,
        output_dir       = OUTPUT_DIR,
        tolerance_arcsec = TOLERANCE_ARCSEC,
        target_star_ra   = TARGET_RA,
        target_star_dec  = TARGET_DEC,
        diagnostics      = False,
    )

    me.run_multi_epoch(
        target_rms               = TARGET_RMS,
        numstars                 = NUM_STARS,
        sdss_filter              = SDSS_FILTER,
        sdss_search_radius_arcsec = SDSS_SEARCH_RADIUS_ARCSEC,
        sdss_data_release        = SDSS_DATA_RELEASE,
        period      = 0.63522741310,   # days
        t0  = 55702.111161463
    )
