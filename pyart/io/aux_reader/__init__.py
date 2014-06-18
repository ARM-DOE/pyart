"""
==================================
Input and output (:mod:`pyart.io.aux_reader`)
==================================

.. currentmodule:: pyart.io.auxreader

Py-ART has modules, classes and functions which are able to read data
from and write data to a number of file formats.


    read_gamic

"""

try:
    from .gamic_hdf5 import read_gamic
    _HDF5_AVAILABLE = True
except:
    _HDF5_AVAILABLE = False

__all__ = [s for s in dir() if not s.startswith('_')]
