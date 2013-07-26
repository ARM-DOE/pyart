"""
pyart.testing.sample_files
==========================

"""

import os

DATA_PATH = os.path.join(os.path.dirname(__file__), 'data')

MDV_PPI_FILE = os.path.join(DATA_PATH, 'example_mdv_ppi.mdv')
MDV_RHI_FILE = os.path.join(DATA_PATH, 'example_mdv_rhi.mdv')
NETCDF_PPI_FILE = os.path.join(DATA_PATH, 'example_netcdf_ppi.nc')
NETCDF_RHI_FILE = os.path.join(DATA_PATH, 'example_netcdf_rhi.nc')
SIGMET_PPI_FILE = os.path.join(DATA_PATH, 'example_sigmet_ppi.sigmet')
SIGMET_RHI_FILE = os.path.join(DATA_PATH, 'example_sigmet_rhi.sigmet')
NEXRAD_ARCHIVE_FILE = os.path.join(DATA_PATH, 'example_nexrad_archive.bz2')
NEXRAD_ARCHIVE_COMPRESSED_FILE = os.path.join(
    DATA_PATH, 'example_nexrad_archive_compressed.ar2v')
NEXRAD_CDM_FILE = os.path.join(DATA_PATH, 'example_nexrad_cdm.bz2')
INTERP_SOUNDE_FILE = os.path.join(DATA_PATH, 'example_interpolatedsonde.cdf')
_EXAMPLE_RAYS_FILE = os.path.join(DATA_PATH, 'example_rays.npz')
