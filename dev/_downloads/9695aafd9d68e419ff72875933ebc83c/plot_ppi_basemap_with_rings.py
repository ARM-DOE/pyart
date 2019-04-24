"""
==============================
Create a PPI plot on a basemap
==============================

An example which creates a PPI plot of a file with a basemap background
and range rings

"""
print(__doc__)

# Author: Scott Collis (scollis@anl.gov)
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt
import pyart

# read in the file, create a RadarMapDisplayBasemap object
filename = 'nsaxsaprppiC1.a1.20140201.184802.nc'
radar = pyart.io.read(filename)
display = pyart.graph.RadarMapDisplayBasemap(radar)

# plot the second tilt
display.plot_ppi_map('reflectivity_horizontal', 1, vmin=-20, vmax=20,
                     min_lon=-157.1, max_lon=-156, min_lat=71.2, max_lat=71.6,
                     lon_lines=np.arange(-158, -154, .2), projection='lcc',
                     lat_lines=np.arange(69, 72, .1), resolution='h',
                     lat_0=radar.latitude['data'][0],
                     lon_0=radar.longitude['data'][0])

# plot range rings at 10, 20, 30 and 40km
display.plot_range_ring(10., line_style='k-')
display.plot_range_ring(20., line_style='k--')
display.plot_range_ring(30., line_style='k-')
display.plot_range_ring(40., line_style='k--')

# plots cross hairs
display.plot_line_xy(np.array([-40000.0, 40000.0]), np.array([0.0, 0.0]),
                     line_style='k-')
display.plot_line_xy(np.array([0.0, 0.0]), np.array([-20000.0, 200000.0]),
                     line_style='k-')

# Indicate the radar location with a point
display.plot_point(radar.longitude['data'][0], radar.latitude['data'][0])

plt.show()
