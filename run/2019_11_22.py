import os

# =====================================================================
# 1. Directory & Paths Setup
# =====================================================================
RAW_DIR_22 = '/lustre/MSSP/sittipong/all/data/2019_11_22'
RAW_DIR_26 = '/lustre/MSSP/sittipong/all/data/2019_11_26' # Added for the cross-night flats

BASE_SAVE_DIR = '/lustre/MSSP/sittipong/reduce2/2019-11-22'  
ASTROMETRY_CACHE = '/lustre/MSSP/sittipong/astrometry_cache'

# =====================================================================
# 2. Input Runs Definition (2019_11_26)
# =====================================================================
INPUT_RUNS = [
    # ------------------ G-BAND ------------------
    {
        'run_name': '2019_11_22_g',
        'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # Master 2x2 bias
        'bias_dark': None, 
        'bias_flat': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run007', 3, 27]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # WARNING: This is a 1x1 flat!
        'data': {
            'runs': [
                ['run010', 1, 80],
                ['run014', 1, 80],
                ['run017', 1, 80],
                ['run020', 1, 80],
                ['run023', 1, 80],
                ['run026', 1, 80],
                ['run029', 1, 70]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_22
        }
    },
    # ------------------ R-BAND ------------------
    {
        'run_name': '2019_11_22_r',
        'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run009', 2, 23]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # WARNING: This is a 1x1 flat!
        'data': {
            'runs': [
                ['run011', 1, 80],
                ['run015', 1, 80],
                ['run018', 1, 80],
                ['run021', 1, 80],
                ['run024', 1, 80],
                ['run027', 1, 80],
                ['run030', 1, 70]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_22
        }
    },
    # ------------------ Z-BAND ------------------
    {
        'run_name': '2019_11_22_z',
        'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_19, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_19, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run014', 2, 20]], 'rawdir': RAW_DIR_19, 'ccd': '1'}, # Good 2x2 flat
        'data': {
            'runs': [
                ['run012', 1, 80],
                ['run016', 1, 80],
                ['run019', 1, 80],
                ['run022', 1, 80],
                ['run025', 1, 80],
                ['run028', 1, 80],
                ['run031', 1, 70],
                ['run033', 1, 6],
                ['run034', 1, 8],
                ['run035', 1, 12],
                # ['run036', 1, 219] # "Fringe frame for z filter"
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_22
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