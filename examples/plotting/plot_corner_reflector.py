"""
=======================
Plot a Corner Reflector
=======================

This is an example of how to plot a corner reflector.

"""

print(__doc__)

# Author: Max Grover (mgrover@anl.gov)
# License: BSD 3 clause

from open_radar_data import DATASETS

import pyart

filename = DATASETS.fetch("sgpkasacrcrrasterC1.a1.20130419.012153.nc")

# Read the data
radar = pyart.io.read(filename)

# Set up the display object
display = pyart.graph.RadarDisplay(radar)

# Set the corner reflector to be 478 m away, with elevations between -0.5 to 2.5 deg
display.plot_cr_raster(target_range=478.0, el_limits=[-0.5, 2.5], cmap="HomeyerRainbow")
