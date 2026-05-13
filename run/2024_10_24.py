import os

# =====================================================================
# 1. Directory & Paths Setup
# =====================================================================
RAW_DIR_24 = '/lustre/MSSP/sittipong/all/data/2024_10_24'
RAW_DIR_26 = '/lustre/MSSP/sittipong/all/data/2024_10_26' # Sourced for bias, r, i, z flats
RAW_DIR_12 = '/lustre/MSSP/sittipong/all/data/2024_05_12' # Sourced for g flat

BASE_SAVE_DIR = '/lustre/MSSP/sittipong/reduce2/2024-10-24'  
ASTROMETRY_CACHE = '/lustre/MSSP/sittipong/astrometry_cache'

# =====================================================================
# 2. Input Runs Definition (2024_10_24)
# =====================================================================
INPUT_RUNS = [
    # ------------------ G-BAND (Science is 2x2) ------------------
    {
        'run_name': '2024_10_24_g',
        'bias_data': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # 2x2 Bias from Oct 26
        'bias_dark': None, 
        'bias_flat': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run040', 7, 47]], 'rawdir': RAW_DIR_12, 'ccd': '1'}, # 2x2 Flat borrowed from May 12
        'data': {
            'runs': [
                ['run020', 1, 20],
                ['run024', 1, 20],
                ['run028', 1, 20],
                ['run032', 1, 20],
                ['run036', 1, 20],
                ['run040', 1, 20]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_24
        }
    },
    # ------------------ R-BAND (Science is 2x2) ------------------
    {
        'run_name': '2024_10_24_r',
        'bias_data': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # 2x2 Bias from Oct 26
        'bias_dark': None, 
        'bias_flat': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run009', 2, 20]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # 2x2 Flat from Oct 26
        'data': {
            'runs': [
                ['run021', 1, 20],
                ['run025', 1, 20],
                ['run029', 1, 20],
                ['run033', 1, 20],
                ['run037', 1, 20],
                ['run041', 1, 20]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_24
        }
    },
    # ------------------ I-BAND (Science is 2x2) ------------------
    {
        'run_name': '2024_10_24_i',
        'bias_data': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # 2x2 Bias from Oct 26
        'bias_dark': None, 
        'bias_flat': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run008', 2, 20]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # 2x2 Flat from Oct 26
        'data': {
            'runs': [
                ['run022', 1, 20],
                ['run026', 1, 20],
                ['run030', 1, 20],
                ['run034', 1, 20],
                ['run038', 1, 20],
                ['run042', 1, 20]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_24
        }
    },
    # ------------------ Z-BAND (Science is 2x2) ------------------
    {
        'run_name': '2024_10_24_z',
        'bias_data': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # 2x2 Bias from Oct 26
        'bias_dark': None, 
        'bias_flat': {'runs': [['run004', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run007', 2, 20]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # 2x2 Flat from Oct 26
        'data': {
            'runs': [
                ['run023', 1, 20],
                ['run027', 1, 20],
                ['run031', 1, 20],
                ['run035', 1, 20],
                ['run039', 1, 20],
                ['run043', 1, 20]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_24
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