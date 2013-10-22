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
    CFRadialDisplay


"""

from .radar_display import RadarDisplay
from .plot_mdv import MdvDisplay
from .plot_cfradial import CFRadialDisplay
from . import cm

try:
    from .gridmapdisplay import GridMapDisplay
except ImportError:
    import warnings
    warnings.warn('No grid plotting support, requires basemap and pyproj')

try:
    from .plot_rsl import RslDisplay
except ImportError:
    pass

__all__ = [s for s in dir() if not s.startswith('_')]
