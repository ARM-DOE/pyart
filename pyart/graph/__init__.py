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
    mapplotgrid_3p


"""

__all__ = ['RadarDisplay', 'MdvDisplay', 'RslDisplay', 'NetcdfDisplay', 
           'mapplotgrid_3p']

from radar_display import RadarDisplay
from plot_mdv import MdvDisplay
try:
    from plot_rsl import RslDisplay
except ImportError:
    pass
from plot_netcdf import NetcdfDisplay
try:
    from grid_plotting import mapplotgrid_3p
except ImportError:
    print('Skipping grid plotting, requires matplotlib and pyproj')
