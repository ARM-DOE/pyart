"""
Creating plots of Radar and Grid fields.

There are also Radar related colormaps and colorblind friendly radar
colormaps for plotting.

Available colormaps, reversed versions (_r) are also provided, these
colormaps are available within matplotlib with names 'pyart_COLORMAP':

    * BlueBrown10
    * BlueBrown11
    * BrBu10
    * BrBu12
    * Bu10
    * Bu7
    * BuDOr12
    * BuDOr18
    * BuDRd12
    * BuDRd18
    * BuGr14
    * BuGy8
    * BuOr10
    * BuOr12
    * BuOr8
    * BuOrR14
    * Carbone11
    * Carbone17
    * Carbone42
    * Cat12
    * EWilson17
    * GrMg16
    * Gray5
    * Gray9
    * NWSRef
    * NWSVel
    * NWS_SPW
    * PD17
    * RRate11
    * RdYlBu11b
    * RefDiff
    * SCook18
    * StepSeq25
    * SymGray12
    * Theodore16
    * Wild25

Colorblind friendly

    * LangRainbow12
    * HomeyerRainbow
    * balance
    * ChaseSpectral
    * SpectralExtended

"""

from . import cm, cm_colorblind  # noqa
from .convstrat_scheme_plot import plot_convstrat_scheme  # noqa
from .gridmapdisplay import GridMapDisplay  # noqa
from .gridmapdisplay_basemap import GridMapDisplayBasemap  # noqa
from .radardisplay import RadarDisplay  # noqa
from .radardisplay_airborne import AirborneRadarDisplay  # noqa
from .radarmapdisplay import RadarMapDisplay  # noqa
from .radarmapdisplay_basemap import RadarMapDisplayBasemap  # noqa

__all__ = [s for s in dir() if not s.startswith("_")]
