"""
Utilities helpful when writing and running unit tests such as sample files
and sample objects.

"""

from .sample_files import MDV_PPI_FILE, MDV_RHI_FILE, MDV_GRID_FILE
from .sample_files import CFRADIAL_PPI_FILE, CFRADIAL_RHI_FILE, CFRADIAL_CR_RASTER_FILE
from .sample_files import CHL_RHI_FILE, UF_FILE
from .sample_files import SIGMET_PPI_FILE, SIGMET_RHI_FILE
from .sample_files import INTERP_SOUNDE_FILE, SONDE_FILE
from .sample_files import NEXRAD_ARCHIVE_MSG31_FILE, NEXRAD_ARCHIVE_MSG1_FILE
from .sample_files import NEXRAD_CDM_FILE
from .sample_files import NEXRAD_ARCHIVE_MSG31_COMPRESSED_FILE
from .sample_files import NEXRAD_LEVEL3_MSG19, NEXRAD_LEVEL3_MSG163
from .sample_files import NEXRAD_LEVEL3_MSG176
from .sample_objects import make_empty_ppi_radar, make_target_radar
from .sample_objects import make_single_ray_radar, make_velocity_aliased_radar
from .sample_objects import make_empty_grid
from .sample_objects import make_target_grid, make_storm_grid
from .sample_objects import make_empty_rhi_radar
from .sample_objects import make_velocity_aliased_rhi_radar
from .sample_objects import make_empty_spectra_radar
from. sample_objects import make_target_spectra_radar
from .tmpdirs import InTemporaryDirectory
from .sample_objects import make_normal_storm

__all__ = [s for s in dir() if not s.startswith('_')]
