import os

# =====================================================================
# 1. Directory & Paths Setup
# =====================================================================
RAW_DIR_08 = '/lustre/MSSP/sittipong/all/data/2018_11_08' # Today's Science Data
RAW_DIR_07 = '/lustre/MSSP/sittipong/all/data/2018_11_07' # Yesterday's Calibrations


BASE_SAVE_DIR = '/lustre/MSSP/sittipong/reduce2/2018_11_08'  
ASTROMETRY_CACHE = '/lustre/MSSP/sittipong/astrometry_cache'

# =====================================================================
# 2. Input Runs Definition (2018_11_08)
# =====================================================================
INPUT_RUNS = [
    # ------------------ G-BAND (Science is 2x2) ------------------
    {
        'run_name': '2018_11_08_g',
        'bias_data': {'runs': [['run064', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'}, # 2x2 bias
        'bias_dark': None, 
        'bias_flat': {'runs': [['run064', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run009', 2, 20]], 'rawdir': RAW_DIR_07, 'ccd': '1'}, # WARNING: This is a 1x1 flat!
        'data': {
            'runs': [
                ['run014', 1, 33],
                ['run018', 1, 38],
                ['run021', 1, 32],
                ['run024', 1, 36],
                ['run027', 1, 36],
                ['run030', 1, 31],
                ['run033', 1, 40],
                ['run036', 1, 36],
                ['run039', 1, 34],
                ['run042', 1, 35],
                ['run045', 1, 41]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_08
        }
    },
    # ------------------ R-BAND (Science is 2x2) ------------------
    {
        'run_name': '2018_11_08_r',
        'bias_data': {'runs': [['run064', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run064', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run063', 15, 98]], 'rawdir': RAW_DIR_07, 'ccd': '1'}, # Excellent 2x2 flat
        'data': {
            'runs': [
                ['run015', 1, 36], # Note: "satellite moving across"
                ['run019', 1, 56],
                ['run022', 1, 37],
                ['run025', 1, 37],
                ['run028', 1, 37],
                ['run031', 1, 36],
                ['run034', 1, 40],
                ['run037', 1, 32],
                ['run040', 1, 35],
                ['run043', 1, 41],
                ['run046', 1, 36]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_08
        }
    },
    # ------------------ Z-BAND (Science is 2x2) ------------------
    {
        'run_name': '2018_11_08_z',
        'bias_data': {'runs': [['run064', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run064', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run060', 19, 37]], 'rawdir': RAW_DIR_07, 'ccd': '1'}, # WARNING: This is a 1x1 flat!
        'data': {
            'runs': [
                ['run017', 1, 29], # Note: run016 excluded (aborted run)
                ['run020', 1, 38],
                ['run023', 1, 34],
                ['run026', 1, 32],
                ['run029', 1, 38],
                ['run032', 1, 33],
                ['run035', 1, 39],
                ['run038', 1, 34],
                ['run041', 1, 34],
                ['run044', 1, 31],
                ['run047', 1, 35]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_08
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
IGNOR_STARS = [5,8]