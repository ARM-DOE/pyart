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
    read_nexrad_archive
    read_nexrad_cdm
    read_gamic
    write_cfradial
    read_grid
    write_grid
    is_vpt
    to_vpt

    Radar
    Grid

"""

try:
    from .rsl import read_rsl
    _RSL_AVAILABLE = True
except ImportError:
    _RSL_AVAILABLE = False
from .mdv import read_mdv
from .sigmet import read_sigmet
from .cfradial import read_cfradial, write_cfradial
from .nexrad_archive import read_nexrad_archive
from .nexrad_cdm import read_nexrad_cdm
try:
    from .gamic_hdf5 import read_gamic
    _HDF5_AVAILABLE = True
except:
    _HDF5_AVAILABLE = False
from .radar import Radar, is_vpt, to_vpt
from .grid import read_grid, write_grid, Grid
from .auto_read import read

__all__ = [s for s in dir() if not s.startswith('_')]
