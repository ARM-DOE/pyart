"""
=================================
Grid Data Using Xradar and Py-ART
=================================

An example which uses xradar and Py-ART to grid a PPI file.

"""

# Author: Max Grover (mgrover@anl.gov)
# License: BSD 3 clause


import xradar as xd

import pyart
from pyart.testing import get_test_data

# Locate the test data and read in using xradar
filename = get_test_data("swx_20120520_0641.nc")
tree = xd.io.open_cfradial1_datatree(filename)

# Give the tree Py-ART radar methods
radar = pyart.xradar.Xradar(tree)

# Grid using 11 vertical levels, and 101 horizontal grid cells at a resolution on 1 km
grid = pyart.map.grid_from_radars(
    (radar,),
    grid_shape=(11, 101, 101),
    grid_limits=(
        (0.0, 10_000),
        (-50_000.0, 50_000.0),
        (-50_000, 50_000.0),
    ),
)

display = pyart.graph.GridMapDisplay(grid)
display.plot_grid(
    "reflectivity_horizontal", level=0, vmin=-20, vmax=60, cmap="pyart_ChaseSpectral"
)
