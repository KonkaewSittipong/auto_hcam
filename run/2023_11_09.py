import os

# =====================================================================
# 1. Directory & Paths Setup
# =====================================================================
RAW_DIR_09 = '/lustre/MSSP/sittipong/all/data/2023_11_09'
RAW_DIR_07 = '/lustre/MSSP/sittipong/all/data/2023_11_07' # For missing flats

BASE_SAVE_DIR = '/lustre/MSSP/sittipong/reduce2/2023-11-09'
ASTROMETRY_CACHE = '/lustre/MSSP/sittipong/astrometry_cache'

# =====================================================================
# 2. Input Runs Definition
# =====================================================================
INPUT_RUNS = [
    {
        'run_name': '2023_11_09_g',
        'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_09, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run002', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run007', 2, 40]], 'rawdir': RAW_DIR_07, 'ccd': '1'}, # From Nov 7
        'data': {
            'runs': [
                ['run008', 1, 20], 
                ['run013', 1, 20],
                ['run017', 1, 20], 
                ['run022', 1, 20], 
                ['run026', 1, 20],
                ['run030', 1, 20], 
                ['run035', 1, 20], 
                ['run040', 1, 20]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_09
        }
    },
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    {
        'run_name': '2023_11_09_r',  
        'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_09, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run002', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run008', 2, 32]], 'rawdir': RAW_DIR_07, 'ccd': '1'}, # From Nov 7
        'data': {
            'runs': [
                ['run009', 1, 20], 
                ['run014', 1, 20], 
                ['run018', 1, 20],
                ['run019', 1, 20], 
                ['run023', 1, 20], 
                ['run027', 1, 20],
                ['run031', 1, 20], 
                ['run036', 1, 20], 
                ['run041', 1, 20]
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_09
        }
    },
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    {
        'run_name': '2023_11_09_i',
        'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_09, 'ccd': '1'},
        'bias_dark': None, 
        'bias_flat': {'runs': [['run002', 2, 0]], 'rawdir': RAW_DIR_07, 'ccd': '1'},
        'dark_data': None,
        'dark_flat': None,
        'flat_data': {'runs': [['run009', 2, 57]], 'rawdir': RAW_DIR_07, 'ccd': '1'}, # From Nov 7
        'data': {
            'runs': [
                ['run010', 1, 20],
                ['run015', 1, 20], # Note: Comment said g but filter listed as i
                ['run020', 1, 20], 
                ['run024', 1, 20], 
                ['run028', 1, 20], 
                ['run032', 1, 20], 
                ['run033', 1, 20], 
                ['run037', 1, 20], 
                ['run042', 1, 20] 
            ], 
            'ccd': '1', 
            'rawdir': RAW_DIR_09
        }
    },
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    # {
    #     'run_name': '2023_11_09_z',
    #     'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_09, 'ccd': '1'},
    #     'bias_dark': None, 
    #     # Using Nov 9 bias for the Nov 9 flat
    #     'bias_flat': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_09, 'ccd': '1'}, 
    #     'dark_data': None,
    #     'dark_flat': None,
    #     'flat_data': {'runs': [['run005', 2, 100]], 'rawdir': RAW_DIR_09, 'ccd': '1'}, # From Nov 9
    #     'data': {
    #         'runs': [
    #             ['run011', 1, 20],
    #             ['run012', 1, 20], 
    #             ['run016', 1, 20], 
    #             ['run021', 1, 20], 
    #             ['run025', 1, 20], 
    #             ['run029', 1, 20], 
    #             ['run034', 1, 20], 
    #             ['run039', 1, 20], 
    #             ['run043', 1, 20] 
    #         ], 
    #         'ccd': '1', 
    #         'rawdir': RAW_DIR_09
    #     }
    # },
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
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