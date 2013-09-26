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
    SIGMET_PPI_FILE
    SIGMET_RHI_FILE
    NEXRAD_ARCHIVE_FILE
    NEXRAD_ARCHIVE_COMPRESSED_FILE
    NEXRAD_CDM_FILE
    INTERP_SOUNDE_FILE

"""

import os

DATA_PATH = os.path.join(os.path.dirname(__file__), 'data')

MDV_PPI_FILE = os.path.join(DATA_PATH, 'example_mdv_ppi.mdv')
MDV_RHI_FILE = os.path.join(DATA_PATH, 'example_mdv_rhi.mdv')
CFRADIAL_PPI_FILE = os.path.join(DATA_PATH, 'example_cfradial_ppi.nc')
CFRADIAL_RHI_FILE = os.path.join(DATA_PATH, 'example_cfradial_rhi.nc')
SIGMET_PPI_FILE = os.path.join(DATA_PATH, 'example_sigmet_ppi.sigmet')
SIGMET_RHI_FILE = os.path.join(DATA_PATH, 'example_sigmet_rhi.sigmet')
NEXRAD_ARCHIVE_FILE = os.path.join(DATA_PATH, 'example_nexrad_archive.bz2')
NEXRAD_ARCHIVE_COMPRESSED_FILE = os.path.join(
    DATA_PATH, 'example_nexrad_archive_compressed.ar2v')
NEXRAD_CDM_FILE = os.path.join(DATA_PATH, 'example_nexrad_cdm.bz2')
INTERP_SOUNDE_FILE = os.path.join(DATA_PATH, 'example_interpolatedsonde.cdf')
_EXAMPLE_RAYS_FILE = os.path.join(DATA_PATH, 'example_rays.npz')
