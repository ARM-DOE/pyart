"""
=============================
Graphing (:mod:`pyart.graph`)
=============================

.. currentmodule:: pyart.graph

.. autosummary::
    :toctree: generated/

    RadarDisplay
    GridMapDisplay
    MdvDisplay
    RslDisplay
    NetcdfDisplay


"""

from .radar_display import RadarDisplay
from .plot_mdv import MdvDisplay
from .plot_netcdf import NetcdfDisplay

try:
    from .gridmapdisplay import GridMapDisplay
except ImportError:
    import warnings
    warnings.warn('No grid plotting support, requires matplotlib and pyproj')

try:
    from .plot_rsl import RslDisplay
except ImportError:
    pass

__all__ = [s for s in dir() if not s.startswith('_')]
