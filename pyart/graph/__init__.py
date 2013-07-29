"""
=============================
Graphing (:mod:`pyart.graph`)
=============================

.. currentmodule:: pyart.graph

.. autosummary::
    :toctree: generated/

    RadarDisplay
    MdvDisplay
    RslDisplay
    NetcdfDisplay


"""

__all__ = ['RadarDisplay', 'MdvDisplay', 'RslDisplay', 'NetcdfDisplay']

from radar_display import RadarDisplay
from plot_mdv import MdvDisplay
try:
    from plot_rsl import RslDisplay
except ImportError:
    import warnings
    warnings.warn('RSL not installed, RslDisplay not available')
from plot_netcdf import NetcdfDisplay
