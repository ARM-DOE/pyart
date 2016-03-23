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
    AirborneRadarDisplay

Plotting grid data
==================

.. autosummary::
    :toctree: generated/

    GridMapDisplay

"""

from .radardisplay import RadarDisplay
from . import cm
from .radardisplay_airborne import AirborneRadarDisplay
from .gridmapdisplay import GridMapDisplay
from .radarmapdisplay import RadarMapDisplay

import warnings as _warnings


class RadarDisplay_Airborne(AirborneRadarDisplay):
    """ Deprecated name for the AirborneRadarDisplay class. """

    def __init__(self, *args, **kwargs):
        _warnings.warn(
            ("'RadarDisplay_Airborne' is deprecated and will be removed in"
             "future versions of Py-ART, please use 'AirborneRadarDisplay'"),
            DeprecationWarning)
        AirborneRadarDisplay.__init__(self, *args, **kwargs)

__all__ = [s for s in dir() if not s.startswith('_')]
