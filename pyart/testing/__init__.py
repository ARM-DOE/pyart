"""
========================================
Testing Utilities (:mod:`pyart.testing`)
========================================

.. currentmodule:: pyart.testing

Utilities helpful when writing and running unit tests.

Testing functions
=================

.. autosummary::
    :toctree: generated/

    make_empty_ppi_radar
    make_target_radar
    make_single_ray_radar
    make_velocity_aliased_radar
    make_empty_grid
    make_target_grid

Testing classes
===============

.. autosummary::
    :toctree: generated/

    InTemporaryDirectory

"""

from .sample_files import MDV_PPI_FILE, MDV_RHI_FILE, MDV_GRID_FILE
from .sample_files import CFRADIAL_PPI_FILE, CFRADIAL_RHI_FILE
from .sample_files import CHL_RHI_FILE, UF_FILE
from .sample_files import SIGMET_PPI_FILE, SIGMET_RHI_FILE
from .sample_files import INTERP_SOUNDE_FILE
from .sample_files import NEXRAD_ARCHIVE_MSG31_FILE, NEXRAD_ARCHIVE_MSG1_FILE
from .sample_files import NEXRAD_CDM_FILE
from .sample_files import NEXRAD_ARCHIVE_MSG31_COMPRESSED_FILE
from .sample_files import NEXRAD_LEVEL3_MSG19, NEXRAD_LEVEL3_MSG163
from .sample_objects import make_empty_ppi_radar, make_target_radar
from .sample_objects import make_single_ray_radar, make_velocity_aliased_radar
from .sample_objects import make_empty_grid
from .sample_objects import make_target_grid, make_storm_grid
from .sample_objects import make_empty_rhi_radar
from .sample_objects import make_velocity_aliased_rhi_radar
from .tmpdirs import InTemporaryDirectory

__all__ = [s for s in dir() if not s.startswith('_')]
