"""
=============================
Graphing (:mod:`pyart.graph`)
=============================

.. currentmodule:: pyart.graph

Creating plots of Radar and Grid fields.

Plotting radar data
===================

.. autosummary::
    :toctree: generated/

    RadarDisplay
    RadarMapDisplay
    RadarDisplay_Airborne

Plotting grid data
==================

.. autosummary::
    :toctree: generated/

    GridMapDisplay

"""

from .radardisplay import RadarDisplay
from . import cm
from .radardisplay_airborne import RadarDisplay_Airborne

try:
    from .gridmapdisplay import GridMapDisplay
except ImportError:
    import warnings
    warnings.warn('No grid plotting support, requires basemap.')

try:
    from .radarmapdisplay import RadarMapDisplay
except ImportError:
    import warnings
    warnings.warn('No grid plotting support, requires basemap.')

__all__ = [s for s in dir() if not s.startswith('_')]
