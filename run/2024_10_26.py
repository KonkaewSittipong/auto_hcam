import os

# =====================================================================
# 1. Directory & Paths Setup
# =====================================================================
RAW_DIR_26 = '/lustre/MSSP/sittipong/all/data/2024_10_26'
RAW_DIR_12 = '/lustre/MSSP/sittipong/all/data/2024_05_12' 
BASE_SAVE_DIR = '/lustre/MSSP/sittipong/reduce2/2024-10-26'  
ASTROMETRY_CACHE = '/lustre/MSSP/sittipong/astrometry_cache'

# =====================================================================
# 2. Input Runs Definition (2024_10_26)
# =====================================================================
INPUT_RUNS = [
    # ------------------ G-BAND ------------------
    {
        'run_name': '2024_10_26_g',
        'bias_data': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # 2x2 Bias
        'bias_dark': None, 
        'bias_flat': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run040', 7, 47]], 'rawdir': RAW_DIR_12, 'ccd': '1'},
        'data': {
            'runs': [
                # --- Full Frame (1056 x 1072) ---
                ['run010', 1, 20], # Note: observer noted bad alignment/elongated profile
                ['run011', 1, 14], # Focus adjustment
                ['run015', 1, 20],
                ['run019', 1, 20],
                # --- Cropped Window (600 x 1072) ---
                ['run020', 1, 20],
                ['run025', 1, 20],
                ['run029', 1, 20],
                ['run032', 1, 20],
                ['run039', 1, 20],
                ['run043', 1, 20],
                ['run047', 1, 20]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_26
        }
    },
    # ------------------ R-BAND ------------------
    {
        'run_name': '2024_10_26_r',
        'bias_data': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run009', 2, 20]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # Good 2x2 flat
        'data': {
            'runs': [
                # --- Full Frame (1056 x 1072) ---
                ['run012', 1, 20],
                ['run016', 1, 20],
                # --- Cropped Window (600 x 1072) ---
                ['run021', 1, 20],
                ['run026', 1, 20],
                ['run030', 1, 20],
                ['run033', 1, 20],
                ['run034', 1, 20], # Test offset
                ['run035', 1, 14], # Test offset
                ['run036', 1, 20],
                ['run040', 1, 20],
                ['run044', 1, 20],
                ['run048', 1, 3]   # Short run at end
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_26
        }
    },
    # ------------------ I-BAND ------------------
    {
        'run_name': '2024_10_26_i',
        'bias_data': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run008', 2, 20]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # Good 2x2 flat
        'data': {
            'runs': [
                # --- Full Frame (1056 x 1072) ---
                ['run013', 1, 20],
                ['run017', 1, 20],
                # --- Cropped Window (600 x 1072) ---
                ['run022', 1, 20],
                ['run023', 1, 20],
                ['run027', 1, 20],
                ['run031', 1, 20],
                ['run037', 1, 20],
                ['run041', 1, 20],
                ['run045', 1, 20]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_26
        }
    },
    # ------------------ Z-BAND ------------------
    {
        'run_name': '2024_10_26_z',
        'bias_data': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run007', 2, 20]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # Good 2x2 flat
        'data': {
            'runs': [
                # --- Full Frame (1056 x 1072) ---
                ['run014', 1, 20],
                ['run018', 1, 20],
                # --- Cropped Window (600 x 1072) ---
                ['run024', 1, 3],  # Note: only 3 frames
                ['run028', 1, 20],
                ['run038', 1, 20],
                ['run042', 1, 20],
                ['run046', 1, 20]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_26
        }
    }
]
# =====================================================================
# 3. Reduction & Photometry Parameters
# =====================================================================
REDUCTION_SOURCE = 'ul'
FIX_PIXEL = False
DIAGNOSTICS = True

# Photometry Aper Setup
CCD_LABEL = '1'
SKIP_BRIGHTEST = 10
SIGMA_THRESHOLD = 7
FRAME = 5
R_SKY1 = 12
R_SKY2 = 15
MARGIN_LEFT = 20
MARGIN_RIGHT = 20
MARGIN_BOTTOM = 20
MARGIN_TOP = 50
REF_INDEX = 2

# genred & extract Setup
R_EXTRACT = [1.8, 3.0, 10.0, 2.5, 12.0, 15.0, 3.0, 15.0, 18.0]

# =====================================================================
# 4. Ensemble Parameters
# =====================================================================
TARGET_RMS = 0.02
NUM_STARS = 12
RA_CENTER = 322.4377
DEC_CENTER = -4.4853
RADIUS = 0.5
SCALE_LOW = 0.40
SCALE_HIGH = 0.990
SINGLE_FRAME = True