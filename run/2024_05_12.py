import os

# =====================================================================
# 1. Directory & Paths Setup
# =====================================================================
RAW_DIR_12 = '/lustre/MSSP/sittipong/all/data/2024_05_12'
RAW_DIR_12 = '/lustre/MSSP/sittipong/all/data/2024_05_12'    ####

BASE_SAVE_DIR = '/lustre/MSSP/sittipong/reduce2/2024-05-12'  
ASTROMETRY_CACHE = '/lustre/MSSP/sittipong/astrometry_cache'

# =====================================================================
# 2. Input Runs Definition
# =====================================================================
INPUT_RUNS = [
    {
        'run_name': '2024_05_12_g',
        'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_12, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_12, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run040', 7, 47]], 'rawdir': RAW_DIR_12, 'ccd': '1'},
        'data': {
            'runs': [
                ['run014', 1, 55],
                ['run018', 1, 20],
                ['run022', 1, 20],
                ['run026', 1, 20],
                ['run030', 1, 20],
                ['run034', 1, 20],
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_12
        }
    },
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    {
        'run_name': '2024_05_12_i',
        'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_12, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_12, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run043', 3, 34]], 'rawdir': RAW_DIR_12, 'ccd': '1'},
        'data': {
            'runs': [
                ['run016', 1, 20],
                ['run020', 1, 20],
                ['run024', 1, 20],
                ['run028', 1, 20],
                ['run032', 1, 20],
                ['run036', 1, 20],
                
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_12
        }
    },
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    {
        'run_name': '2024_05_12_r',
        'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_12, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_12, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run041', 2, 40]], 'rawdir': RAW_DIR_12, 'ccd': '1'},
        'data': {
            'runs': [
                ['run015', 1, 20], #
                ['run019', 1, 20],
                ['run023', 1, 20],
                ['run027', 1, 20],
                ['run031', 1, 20],
                ['run035', 1, 20],
                
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_12
        }
    },
    # {
    #     'run_name': '2024_05_12_z',
    #     'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_12, 'ccd': '1'},
    #     'bias_dark': None, 
    #     'bias_flat': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_12, 'ccd': '1'},
    #     'dark_data': None,
    #     'dark_flat': None,
    #     'flat_data': {'runs': [['run045', 2, 68]], 'rawdir': RAW_DIR_12, 'ccd': '1'}, # Selected optimal flat run
    #     'data': {
    #         'runs': [
    #             ['run017', 1, 20],
    #             ['run021', 1, 20],
    #             ['run025', 1, 20],
    #             ['run029', 1, 20],
    #             ['run033', 1, 20]
    #         ], 
    #         'ccd': '1', 
    #         'rawdir': RAW_DIR_12
    #     }
    # }
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