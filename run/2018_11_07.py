import os

# =====================================================================
# 1. Directory & Paths Setup
# =====================================================================
RAW_DIR_07 = '/lustre/MSSP/sittipong/all/data/2018_11_07'
BASE_SAVE_DIR = '/lustre/MSSP/sittipong/reduce2/2018_11_07'  
ASTROMETRY_CACHE = '/lustre/MSSP/sittipong/astrometry_cache'

# =====================================================================
# 2. Input Runs Definition (2018_11_07)
# =====================================================================
INPUT_RUNS = [
    # ------------------ G-BAND (Science is 2x2) ------------------
    {
        'run_name': '2018_11_07_g',
        'bias_data': {'runs': [['run064', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'}, # 2x2 bias
        'bias_dark': None, 
        'bias_flat': {'runs': [['run064', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run009', 2, 20]], 'rawdir': RAW_DIR_07, 'ccd': '1'}, # WARNING: This is a 1x1 flat!
        'data': {
            'runs': [
                ['run017', 1, 50],
                ['run020', 1, 27],
                ['run023', 1, 29],
                ['run026', 1, 32],
                ['run029', 1, 39],
                ['run032', 1, 34],
                ['run035', 1, 32],
                ['run038', 1, 31]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_07
        }
    },
    # ------------------ R-BAND (Science is 2x2) ------------------
    {
        'run_name': '2018_11_07_r',
        'bias_data': {'runs': [['run064', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run064', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run063', 15, 98]], 'rawdir': RAW_DIR_07, 'ccd': '1'}, # Excellent 2x2 flat
        'data': {
            'runs': [
                ['run018', 1, 38],
                ['run021', 1, 31],
                ['run024', 1, 34],
                ['run027', 1, 33],
                ['run030', 1, 26],
                ['run033', 1, 29],
                ['run036', 1, 31],
                ['run039', 1, 35]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_07
        }
    },
    # ------------------ Z-BAND (Science is 2x2) ------------------
    {
        'run_name': '2018_11_07_z',
        'bias_data': {'runs': [['run064', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run064', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run060', 19, 37]], 'rawdir': RAW_DIR_07, 'ccd': '1'}, # WARNING: This is a 1x1 flat!
        'data': {
            'runs': [
                ['run019', 1, 33],
                ['run022', 1, 26],
                ['run025', 1, 29],
                ['run028', 1, 29],
                ['run031', 1, 31],
                ['run034', 1, 34],
                ['run037', 1, 29],
                ['run040', 1, 30]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_07
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
IGNOR_STARS = [4]