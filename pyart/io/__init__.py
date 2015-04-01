"""
==================================
Input and output (:mod:`pyart.io`)
==================================

.. currentmodule:: pyart.io

Py-ART has modules, classes and functions which are able to read data
from and write data to a number of file formats.

.. autosummary::
    :toctree: generated/

    read
    read_rsl
    read_mdv
    read_sigmet
    read_cfradial
    read_chl
    read_nexrad_archive
    read_nexrad_cdm
    read_nexrad_level3
    write_cfradial
    read_grid
    write_grid

"""

try:
    from .rsl import read_rsl
    _RSL_AVAILABLE = True
except ImportError:
    _RSL_AVAILABLE = False
from .mdv import read_mdv
from .sigmet import read_sigmet
from .chl import read_chl
from .cfradial import read_cfradial, write_cfradial
from .nexrad_archive import read_nexrad_archive
from .nexrad_cdm import read_nexrad_cdm
from .nexradl3_read import read_nexrad_level3
from .grid_io import read_grid, write_grid
from .auto_read import read

__all__ = [s for s in dir() if not s.startswith('_')]
