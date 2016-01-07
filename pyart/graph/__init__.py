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
from .radardisplay_airborne import AirborneRadarDisplay
from .radardisplay_airborne import  AirborneRadarDisplay as RadarDisplay_Airborne
from .gridmapdisplay import GridMapDisplay
from .radarmapdisplay import RadarMapDisplay

__all__ = [s for s in dir() if not s.startswith('_')]
