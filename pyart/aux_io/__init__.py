"""
================================================
Auxiliary input and output (:mod:`pyart.aux_io`)
================================================

.. currentmodule:: pyart.aux_io

Additional classes and functions for reading and writing data from a number
of file formats.  These auxiliary input/output routines are not as well
polished as those in :mod:`pyart.io`.  They may require addition dependencies
beyond those required for a standard Py-ART install, use non-standard
function parameter and naming, are not supported by the
:py:func:`pyart.io.read` function and are not fully tested if tested at all.
Please use these at your own risk.  Bugs in these function should be
reported but fixing them will not be a priority.

.. autosummary::
    :toctree: generated/

    read_d3r_gcpex_nc
    read_gamic
    read_pattern
    read_radx

"""

from .pattern import read_pattern
from .radx import read_radx
from .d3r_gcpex_nc import read_d3r_gcpex_nc
try:
    from .gamic_hdf5 import read_gamic
    from .odim_h5 import read_odim_h5
    _HDF5_AVAILABLE = True
except:
    _HDF5_AVAILABLE = False

__all__ = [s for s in dir() if not s.startswith('_')]
