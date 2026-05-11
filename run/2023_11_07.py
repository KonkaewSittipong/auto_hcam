import os

# =====================================================================
# 1. Directory & Paths Setup
# =====================================================================
DEFAULT_RAW_DIR = '/lustre/MSSP/sittipong/all/data/2023_11_07'
BASE_SAVE_DIR = '/lustre/MSSP/sittipong/reduce2/2023-11-07'
ASTROMETRY_CACHE = '/lustre/MSSP/sittipong/astrometry_cache'

# =====================================================================
# 2. Input Runs Definition
# =====================================================================
INPUT_RUNS = [
    {
        'run_name': '2023_11_07_g',
        'bias_data': {'runs': [['run002', 2, 0]], 'rawdir': DEFAULT_RAW_DIR, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': None,
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run007', 2, 40]], 'ccd': '1', 'rawdir': DEFAULT_RAW_DIR},
        'data': {
            'runs': [
                ['run013', 1, 50], 
                ['run018', 15, 20], 
                ['run022', 1, 30], 
                ['run025', 1, 30],
                ['run028', 1, 20],
                ['run031', 1, 20], 
                ['run034', 13, 20],
                ['run037', 1, 20],
                ['run042', 1, 20]
            ], 
            'ccd': '1', 
            'rawdir': DEFAULT_RAW_DIR
        }
    },
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    {
        'run_name': '2023_11_07_i',  
        'bias_data': {'runs': [['run002', 10, 0]], 'ccd': '1', 'rawdir': DEFAULT_RAW_DIR},
        'bias_dark': None, 
        'bias_flat': None,
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run009', 2, 57]], 'ccd': '1', 'rawdir': DEFAULT_RAW_DIR},
        'data': {
            'runs': [
                ['run017', 1, 20], 
                ['run021', 1, 30], 
                ['run024', 1, 30], 
                ['run027', 1, 30],
                ['run030', 1, 20], 
                ['run033', 1, 20], 
                ['run036', 1, 20], 
                ['run041', 1, 0],
                ['run044', 1, 14]
            ], 
            'ccd': '1', 
            'rawdir': DEFAULT_RAW_DIR
        }
    },
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    {
        'run_name': '2023_11_07_r',
        'bias_data': {'runs': [['run002', 2, 0]], 'rawdir': DEFAULT_RAW_DIR, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': None,
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run007', 2, 32]], 'ccd': '1', 'rawdir': DEFAULT_RAW_DIR},
        'data': {
            'runs': [
                ['run014', 1, 11],   # run014: Short run
                ['run016', 1, 20],   # run016
                ['run020', 1, 30],   # run020
                ['run023', 1, 30],   # run023
                ['run026', 1, 30],   # run026
                ['run029', 1, 20],   # run029
                ['run032', 1, 20],   # run032
                ['run035', 1, 20],   # run035
                ['run043', 1, 20]    # run043
            ], 
            'ccd': '1', 
            'rawdir': DEFAULT_RAW_DIR
        }
    }
]

# =====================================================================
# 3. Reduction & Photometry Parameters
# =====================================================================
REDUCTION_SOURCE = 'ul'
FIX_PIXEL = False
DIAGNOSTICS = False

# Photometry Aper Setup
CCD_LABEL = '1'
SKIP_BRIGHTEST = 10
SIGMA_THRESHOLD = 7
FRAME = 3
R_SKY1 = 15
R_SKY2 = 20
MARGIN_LEFT = 20
MARGIN_RIGHT = 20
MARGIN_BOTTOM = 20
MARGIN_TOP = 50
REF_INDEX = 2

# genred & extract Setup
R_EXTRACT = [1.8, 3.0, 15.0, 2.5, 15.0, 18.0, 3.0, 18.0, 20.0]

# =====================================================================
# 4. Ensemble Parameters
# =====================================================================
TARGET_RMS = 0.02
NUM_STARS = 12
RA_CENTER = 322.437695506
DEC_CENTER = -4.485270444
RADIUS = 1.0
SCALE_LOW = 0.40
SCALE_HIGH = 0.990
SINGLE_FRAME = True