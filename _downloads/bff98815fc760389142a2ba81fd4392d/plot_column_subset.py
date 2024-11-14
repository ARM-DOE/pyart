"""
====================================
Extract a radar column above a point
====================================

Given a radar and a point, extract the column of radar data values above
a point
"""

# Author: Maxwell Grover (mgrover@anl.gov)
# License: BSD 3 clause

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

import pyart
from pyart.testing import get_test_data

# Read in some test data
filename = get_test_data("swx_20120520_0641.nc")
radar = pyart.io.read(filename)

######################################
# **Plot the first sweep and our desired point**
#
# Let's visualize our radar data from a single sweep, and plot
# the location of our desired point on a map.
# This will provide some context as to where we are extracting our
# column of values.

site_lon = -97.73  # longitude in degrees
site_lat = 36.41  # latitdue in degrees

# Setup the RadarMapDisplay and add our projection
display = pyart.graph.RadarMapDisplay(radar)
ax = plt.subplot(111, projection=ccrs.PlateCarree())

# Visualize the reflectivity field, using the lowest sweep with
# latitude and longitude lines
display.plot_ppi_map(
    "reflectivity_horizontal",
    0,
    ax=ax,
    vmin=-32,
    vmax=64.0,
    lon_lines=np.arange(-98, -97, 0.2),
    lat_lines=np.arange(36, 37, 0.2),
)

# Plot our site location on top of the radar image
ax.scatter(site_lon, site_lat, color="black")

######################################
# Now that we have our point defined, and our radar object, we can use the following
# utility function in Py-ART to subset a column
ds = pyart.util.columnsect.get_field_location(radar, site_lat, site_lon)

######################################
# This function returns an xarray dataset, with all of our data fields!
print(ds)

######################################
# **Visualize the Reflectivity Values in the Column**
#
# Let's visualize the reflectivity values in the column
# above our point, which is stored in our new dataset
ds.corrected_reflectivity_horizontal.plot(y="height")
