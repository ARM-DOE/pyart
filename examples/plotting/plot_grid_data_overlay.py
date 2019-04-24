"""
====================================
Create a grid plot with data overlay
====================================

An example which creates a plot of a gridded NEXRAD radar on a map
with latitude and NCEP North American regional reanalysis (NARR) pressure
is plotted on top of the grid.

"""
print(__doc__)

# Author Jonathan J. Helmus, Cory Weber
# License: BSD 3 clause

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from netCDF4 import num2date, date2num, Dataset
import pyart


# read in the NEXRAD data, create the display
fname = '20110520100000_nexrad_grid.nc'
grid = pyart.io.read_grid(fname)
display = pyart.graph.GridMapDisplayBasemap(grid)

# create the figure
font = {'size': 10}
matplotlib.rc('font', **font)
fig = plt.figure(figsize=[10, 8])

# Add Basic Title
title = 'Basic Plot with Overlay Example Title'
#     Xleft%, ybot%
fig.text(0.5, 0.9, title, horizontalalignment='center', fontsize=24)

# panel sizes      xleft%, ybot%, xright% ,ytop%
map_panel_axes = [0.05, 0.15, 0.9, 0.7]
colorbar_panel_axes = [0.15, 0.09, 0.7, .010]

# parameters
level = 5
vmin = -8
vmax = 64
lat = 36.5
lon = -98.5

# panel 1, basemap, radar reflectivity and NARR overlay
ax1 = fig.add_axes(map_panel_axes)
display.plot_grid('REF', level=level, vmin=vmin, vmax=vmax, title_flag=False,
                  colorbar_flag=False)

# load overlay data
url = 'narr-a_221_20110520_0000_000.nc'
data = Dataset(url)

# extract data at correct time
grid_date = num2date(grid.time['data'], grid.time['units'])[0]
data_time = data.variables['time']
t_idx = abs(data_time[:] - date2num(grid_date, data_time.units)).argmin()
prmsl = 0.01 * data.variables['prmsl'][t_idx]

# plot the reanalysis on the basemap
lons, lats = np.meshgrid(data.variables['lon'], data.variables['lat'][:])

x, y = display.basemap(lons, lats)
clevs = np.arange(900, 1100., 1.)

display.basemap.contour(x, y, prmsl, clevs, colors='k', linewidths=1.)

# colorbar
cbax = fig.add_axes(colorbar_panel_axes)
display.plot_colorbar(cax=cbax)

# enable below to add cross hairs
# display.plot_crosshairs(lon=lon, lat=lat)

plt.show()
