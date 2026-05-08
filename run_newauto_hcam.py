# from orderHIPERCAM import Hipercam_setup, Photometry
# import os

# # # set inputs
# input_run = {
#     'bias_data': {
#         'runs': [['run002', 10, 0]], 
#         'rawdir': '/lustre/MSSP/sittipong/all/data/2023_11_07',
#         'ccd': '1'
#     },
#     'bias_dark': None, 
#     'bias_flat': None, # {
#         # 'runs': [['run002', 2, 0]], 
#         # 'ccd': '1', 
#     #     # 'rawdir': '/lustre/MSSP/sittipong/all/data/2023_11_07'
#     # },
#     'dark_data': None,
#     'dark_flat': None,
#     'flat_data': {
#         'runs': [['run009', 2, 57]], 
#         'ccd': '1', 
#         'rawdir': '/lustre/MSSP/sittipong/all/data/2023_11_07'
#     },
#     'data': {
#         'runs': [
#             # ['run010', 1, 30], 
#             ['run014', 1, 11], 
#             # ['run018', 1, 20],
#             ['run022', 1, 20], ['run026', 1, 20], ['run030', 1, 20],
#             ['run037', 1, 20], ['run041', 1, 15]
#         ],
#         'ccd': '1', 
#         'rawdir': '/lustre/MSSP/sittipong/all/data/2023_11_07'
#     }
# }
# save_directory = os.path.join('/lustre/MSSP/sittipong/buildmodule/temps', 'hcam_reduction')
# default_raw_directory = '/lustre/MSSP/sittipong/all/data/2023_11_07'
# ASTROMETRY_CACHE = '/lustre/MSSP/sittipong/astrometry_cache'
# R_EXTARCT = [2,  3.0, 10.0,
#              2.5, 12.0, 18.0,
#              3.0, 18.0, 20.0 ]  # unbin radius
# # # Reduction path
# print(f"---  Initializing Pipeline ---")
# pipeline = Hipercam_setup(
#     source='ul', 
#     input_run_ul=input_run, 
#     ccd='1', 
#     save_dir=save_directory, 
#     raw_dir=default_raw_directory,
#     diagnostics=True
# )

# # # Photometry path
# photo_reduction = Photometry(base_dir=save_directory, lis=pipeline.data['lis'], diagnostics=True)

# photo_reduction.setaper(
#     ccd_label='1', SKIP_BRIGHTEST=10, SIGMA_THRESHOLD=7,
#     frame=5, R_SKY1=15, R_SKY2=20, frame=3,
#     MARGIN_LEFT=17, MARGIN_RIGHT=17, MARGIN_BOTTOM=10, MARGIN_TOP=30,
#     diagnostics=True
# )

# photo_reduction.genred()

# photo_reduction.modify_red(
#     fit_method = 'moffat',
#     fit_height_min_ref = 5.0,
#     fit_height_min_nrf = 3.0,
#     fit_half_width = 25.0,
#     search_half_width = 15.0,
#     )

# photo_reduction.modify_red(
#     target_section="extraction",
#     **{"1": f"variable normal {R_EXTARCT[0]} {R_EXTARCT[1]} {R_EXTARCT[2]} {R_EXTARCT[3]} {R_EXTARCT[4]} {R_EXTARCT[5]} {R_EXTARCT[6]} {R_EXTARCT[7]} {R_EXTARCT[8]}"}
# )

# photo_reduction.modify_red(fit_max_shift=8, search_smooth_fft='yes')

# photo_reduction.reduce(R_EXTARCT = R_EXTARCT,
#                        # plot_all = True,
#                        plot_zoom = True,
#                        diagnostics=True)

# photo_reduction.solvwcs(          # reads all.log -> wcs_radec.csv
#     ccd_label='1',
#     # scale_low=0.3,
#     # scale_high=0.9,
#     ra_center=322.4377,
#     dec_center=-4.4853,
#     radius=0.5,
#     astrometry_cache= ASTROMETRY_CACHE,
# )




print('xx'*100)
from orderHIPERCAM import Reduction, Photometry

import os

# set inputs
input_run = {
    'bias_data': {
        'runs': [['run002', 10, 0]], 
        'rawdir': '/lustre/MSSP/sittipong/all/data/2023_11_07',
        'ccd': '1'
    },
    'bias_dark': None, 
    'bias_flat': None,
    'dark_data': None,
    'dark_flat': None,
    'flat_data': {
        'runs': [['run009', 2, 57]], 
        'ccd': '1', 
        'rawdir': '/lustre/MSSP/sittipong/all/data/2023_11_07'
    },
    'data': {
        'runs': [
            ['run014', 1, 11], 
            ['run022', 1, 20], ['run026', 1, 20], ['run030', 1, 20],
            ['run037', 1, 20], ['run041', 1, 15]
        ],
        'ccd': '1', 
        'rawdir': '/lustre/MSSP/sittipong/all/data/2023_11_07'
    }
}

save_directory = os.path.join('/lustre/MSSP/sittipong/buildmodule/temps', 'hcam_reduction')
default_raw_directory = '/lustre/MSSP/sittipong/all/data/2023_11_07'
ASTROMETRY_CACHE = '/lustre/MSSP/sittipong/astrometry_cache'
R_EXTARCT = [2.2,  3.0, 10.0,
             2.5, 12.0, 18.0,
             3.0, 18.0, 20.0]

print(f"---  Initializing Pipeline ---")

pipeline = Reduction(
    source='ul',
    input_run_ul=input_run,
    save_dir=save_directory,
    raw_dir=default_raw_directory,
    fix_pixel=True,          # detect + fix bad pixels
    diagnostics=True
)

photo_reduction = Photometry(
    base_dir=save_directory,
    lis=pipeline.data['lis'],
    bad_pixel_masks=pipeline.bad_pixel_masks,  # pass mask for overlay
    diagnostics=True
)

photo_reduction.setaper(
    ccd_label='1', SKIP_BRIGHTEST=10, SIGMA_THRESHOLD=7,
    frame=3,                           # fix 3: removed duplicate frame=5
    R_SKY1=15, R_SKY2=20,
    MARGIN_LEFT=17, MARGIN_RIGHT=17, MARGIN_BOTTOM=10, MARGIN_TOP=30,
    diagnostics=True
)

photo_reduction.genred()

photo_reduction.modify_red(
    fit_method='moffat',
    fit_height_min_ref=5.0,
    fit_height_min_nrf=3.0,
    fit_half_width=25.0,
    search_half_width=15.0,
)

photo_reduction.modify_red(
    target_section="extraction",
    **{"1": f"variable normal {R_EXTARCT[0]} {R_EXTARCT[1]} {R_EXTARCT[2]} {R_EXTARCT[3]} {R_EXTARCT[4]} {R_EXTARCT[5]} {R_EXTARCT[6]} {R_EXTARCT[7]} {R_EXTARCT[8]}"}
)

photo_reduction.modify_red(fit_max_shift=8, search_smooth_fft='yes')

photo_reduction.reduce(
    R_EXTARCT=R_EXTARCT,
    plot_with_zoom=True,               # fix 4: plot_zoom → plot_with_zoom
    diagnostics=True
)