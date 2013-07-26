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
    read_netcdf
    read_nexrad_archive
    read_nexrad_cdm
    write_netcdf
    read_grid
    write_grid

    Radar
    Grid

"""

from .auto_read import read
from .rsl import read_rsl
from .mdv import read_mdv
from .sigmet import read_sigmet
from .netcdf import read_netcdf, write_netcdf
from .nexrad_archive import read_nexrad_archive
from .nexrad_cdm import read_nexrad_cdm
from .radar import Radar
from .grid import read_grid, write_grid, Grid

__all__ = [s for s in dir() if not s.startswith('_')]
