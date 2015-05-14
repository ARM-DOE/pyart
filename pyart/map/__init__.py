"""
==========================
Mapping (:mod:`pyart.map`)
==========================

.. current modules:: pyart.map

Py-ART has a robust function for mapping radar data from the collected radar
coordinates to Cartesian coordinates.

.. autosummary::
    :toctree: generated/

    grid_from_radars
    map_to_grid
    map_gates_to_grid
    example_roi_func_constant
    example_roi_func_dist
    example_roi_func_dist_beam

"""

from .grid_mapper import map_to_grid, grid_from_radars
from .grid_mapper import example_roi_func_constant
from .grid_mapper import example_roi_func_dist
from .grid_mapper import example_roi_func_dist_beam
from .gates_to_grid import map_gates_to_grid

__all__ = [s for s in dir() if not s.startswith('_')]
