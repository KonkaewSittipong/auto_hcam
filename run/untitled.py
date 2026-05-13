# import os

# # =====================================================================
# # 1. Directory & Paths Setup
# # =====================================================================
# RAW_DIR_24 = '/lustre/MSSP/sittipong/all/data/2024_10_24'
# RAW_DIR_26 = '/lustre/MSSP/sittipong/all/data/2024_10_26' # Sourced for bias, r, i, z flats
# RAW_DIR_12 = '/lustre/MSSP/sittipong/all/data/2024_05_12' # Sourced for g flat

# BASE_SAVE_DIR = '/lustre/MSSP/sittipong/reduce2/2024-10-24'  
# ASTROMETRY_CACHE = '/lustre/MSSP/sittipong/astrometry_cache'

# # =====================================================================
# # 2. Input Runs Definition (2024_10_24)
# # =====================================================================

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