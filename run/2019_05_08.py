# import os

# # =====================================================================
# # 1. Directory & Paths Setup
# # =====================================================================
# RAW_DIR_08 = '/lustre/MSSP/sittipong/all/data/2019_05_08'
# RAW_DIR_26 = '/lustre/MSSP/sittipong/all/data/2019_11_26' # Added for the cross-night flats
# BASE_SAVE_DIR = '/lustre/MSSP/sittipong/reduce2/2019-05-08'  
# ASTROMETRY_CACHE = '/lustre/MSSP/sittipong/astrometry_cache'

# # =====================================================================
# # 2. Input Runs Definition (2019_11_26)
# # =====================================================================
# INPUT_RUNS = [
#     # ------------------ G-BAND ------------------
# {
#         'run_name': '2019_05_08_g',
#         'bias_data': {'runs': [['run044', 1, 0]], 'rawdir': RAW_DIR_08, 'ccd': '1'}, # Master 1x1 bias
#         'bias_dark': None, 
#         'bias_flat': {'runs': [['run044', 1, 0]], 'rawdir': RAW_DIR_08, 'ccd': '1'},
#         'dark_data': None,
#         'dark_flat': None,
#         'flat_data': {'runs': [['run007', 3, 27]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # WARNING: This is a 1x1 flat!
#         'data': {
#             'runs': [
#                 # run014, run015, run016 skipped (caution/test)
#                 ['run019', 1, 3],
#                 ['run022', 1, 3]
#             ], 
#             'ccd': '1', 
#             'rawdir': RAW_DIR_08
#         }
#     },
#     # ------------------ R-BAND ------------------
#     {
#         'run_name': '2019_05_08_r',
#         'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
#         'bias_dark': None, 
#         'bias_flat': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_26, 'ccd': '1'},
#         'dark_data': None,
#         'dark_flat': None,
#         'flat_data': {'runs': [['run009', 2, 23]], 'rawdir': RAW_DIR_26, 'ccd': '1'}, # WARNING: This is a 1x1 flat!
#         'data': {
#             'runs': [
#                 ['run017', 1, 3],
#                 ['run020', 1, 3],
#                 ['run023', 1, 3]
#             ], 
#             'ccd': '1', 
#             'rawdir': RAW_DIR_08
#         }
#     },
#     # ------------------ Z-BAND ------------------
#     {
#         'run_name': '2019_05_08_z',
#         'bias_data': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_19, 'ccd': '1'},
#         'bias_dark': None, 
#         'bias_flat': {'runs': [['run003', 2, 0]], 'rawdir': RAW_DIR_19, 'ccd': '1'},
#         'dark_data': None,
#         'dark_flat': None,
#         'flat_data': {'runs': [['run014', 2, 20]], 'rawdir': RAW_DIR_19, 'ccd': '1'}, # Good 2x2 flat
#         'data': {
#             'runs': [
#                 ['run018', 1, 3],
#                 ['run021', 1, 3],
#                 ['run024', 1, 3],
#             ], 
#             'ccd': '1', 
#             'rawdir': RAW_DIR_08
#         }
#     }
# ]
# # =====================================================================
# # 3. Reduction & Photometry Parameters
# # =====================================================================
# REDUCTION_SOURCE = 'ul'
# FIX_PIXEL = False
# DIAGNOSTICS = True

# # Photometry Aper Setup
# CCD_LABEL = '1'
# SKIP_BRIGHTEST = 10
# SIGMA_THRESHOLD = 7
# FRAME = 5
# R_SKY1 = 12
# R_SKY2 = 15
# MARGIN_LEFT = 20
# MARGIN_RIGHT = 20
# MARGIN_BOTTOM = 20
# MARGIN_TOP = 50
# REF_INDEX = 2

# # genred & extract Setup
# R_EXTRACT = [1.8, 3.0, 10.0, 2.5, 12.0, 15.0, 3.0, 15.0, 18.0]

# # =====================================================================
# # 4. Ensemble Parameters
# # =====================================================================
# TARGET_RMS = 0.02
# NUM_STARS = 12
# RA_CENTER = 322.4377
# DEC_CENTER = -4.4853
# RADIUS = 0.5
# SCALE_LOW = 0.40
# SCALE_HIGH = 0.990
# SINGLE_FRAME = True