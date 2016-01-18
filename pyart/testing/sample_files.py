"""
pyart.testing.sample_files
==========================

Sample radar files in a number of formats.  Many of these files
are incomplete, they should only be used for testing, not production.

.. autosummary::
    :toctree: generated/

    MDV_PPI_FILE
    MDV_RHI_FILE
    CFRADIAL_PPI_FILE
    CFRADIAL_RHI_FILE
    CHL_RHI_FILE
    SIGMET_PPI_FILE
    SIGMET_RHI_FILE
    NEXRAD_ARCHIVE_MSG31_FILE
    NEXRAD_ARCHIVE_MSG31_COMPRESSED_FILE
    NEXRAD_ARCHIVE_MSG1_FILE
    NEXRAD_LEVEL3_MSG19
    NEXRAD_LEVEL3_MSG163
    NEXRAD_CDM_FILE
    UF_FILE
    INTERP_SOUNDE_FILE

"""

import os

DATA_PATH = os.path.join(os.path.dirname(__file__), 'data')

MDV_PPI_FILE = os.path.join(DATA_PATH, 'example_mdv_ppi.mdv')
MDV_RHI_FILE = os.path.join(DATA_PATH, 'example_mdv_rhi.mdv')
MDV_GRID_FILE = os.path.join(DATA_PATH, 'example_mdv_grid.mdv')
CFRADIAL_PPI_FILE = os.path.join(DATA_PATH, 'example_cfradial_ppi.nc')
CFRADIAL_RHI_FILE = os.path.join(DATA_PATH, 'example_cfradial_rhi.nc')
CHL_RHI_FILE = os.path.join(DATA_PATH, 'example_chl_rhi.chl')
SIGMET_PPI_FILE = os.path.join(DATA_PATH, 'example_sigmet_ppi.sigmet')
SIGMET_RHI_FILE = os.path.join(DATA_PATH, 'example_sigmet_rhi.sigmet')
NEXRAD_ARCHIVE_MSG31_FILE = os.path.join(
    DATA_PATH, 'example_nexrad_archive_msg31.bz2')
NEXRAD_ARCHIVE_MSG31_COMPRESSED_FILE = os.path.join(
    DATA_PATH, 'example_nexrad_archive_msg31_compressed.ar2v')
NEXRAD_ARCHIVE_MSG1_FILE = os.path.join(
    DATA_PATH, 'example_nexrad_archive_msg1.bz2')
# NEXRAD Level 3 file downloaded from NCDC with filenames:
# KBMX_SDUS54_N0RBMX_201501020205
# KBMX_SDUS84_N0KBMX_201501020205
NEXRAD_LEVEL3_MSG19 = os.path.join(DATA_PATH, 'example_nexrad_level3_msg19')
NEXRAD_LEVEL3_MSG163 = os.path.join(DATA_PATH, 'example_nexrad_level3_msg163')
NEXRAD_CDM_FILE = os.path.join(DATA_PATH, 'example_nexrad_cdm.bz2')
UF_FILE = os.path.join(DATA_PATH, 'example_uf_ppi.uf')
INTERP_SOUNDE_FILE = os.path.join(DATA_PATH, 'example_interpolatedsonde.cdf')
_EXAMPLE_RAYS_FILE = os.path.join(DATA_PATH, 'example_rays.npz')
