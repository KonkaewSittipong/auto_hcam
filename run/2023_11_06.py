import os

# =====================================================================
# 1. Directory & Paths Setup
# =====================================================================
RAW_DIR_06 = '/lustre/MSSP/sittipong/all/data/2023_11_06'
RAW_DIR_07 = '/lustre/MSSP/sittipong/all/data/2023_11_07'

BASE_SAVE_DIR = '/lustre/MSSP/sittipong/reduce2/2023-11-06'
ASTROMETRY_CACHE = '/lustre/MSSP/sittipong/astrometry_cache'

# =====================================================================
# 2. Input Runs Definition
# =====================================================================
INPUT_RUNS = [
    {
        'run_name': '2023_11_06_g',
        'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_06, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run002', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run007', 2, 40]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'data': {
            'runs': [
                ['run007', 1, 30], 
                ['run010', 1, 30], 
                ['run013', 1, 30],
                ['run016', 1, 30],
                ['run019', 1, 30], 
                ['run022', 1, 30], 
                ['run025', 1, 30],
                ['run028', 1, 30], 
                ['run031', 1, 30], 
                ['run032', 1, 15]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_06
        }
    },
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    {
        'run_name': '2023_11_06_i',  
        'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_06, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run002', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run009', 2, 57]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'data': {
            'runs': [
                ['run009', 1, 30], 
                ['run012', 1, 30], 
                ['run015', 1, 30],
                ['run018', 1, 30], 
                ['run021', 1, 30], 
                ['run024', 1, 30],
                ['run027', 1, 30], 
                ['run030', 1, 30], 
                ['run034', 1, 15]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_06
        }
    },
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    {
        'run_name': '2023_11_06_r',
        'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_06, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run002', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run008', 2, 32]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'data': {
            'runs': [
                ['run008', 1, 30],
                ['run011', 1, 30], 
                ['run014', 1, 30], 
                # ['run017', 1, 30], # Note: marked as "cloudy" in log
                ['run020', 1, 30], 
                ['run023', 1, 30], 
                ['run026', 1, 30], 
                ['run029', 1, 30], 
                ['run033', 1, 15]  # Note: exposure changed to 20s here
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_06
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