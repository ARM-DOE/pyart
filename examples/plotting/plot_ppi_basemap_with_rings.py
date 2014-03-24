"""
==============================
Create a PPI plot on a basemap
==============================

An example which creates a PPI plot of a file with a basemap background
and range rings

"""
print __doc__

# Author: Scott Collis (scollis@anl.gov)
# License: BSD 3 clause

import matplotlib.pyplot as plt
import pyart

filename = 'nsaxsaprppiC1.a1.20140201.184802.nc'

#Radar gate spacing is wrong in the file is should be 60m

radar.range['data'] = 60.0*radar.range['data']/50.0

font = {'size' : 16}
rc('font', **font)

#create the figure

f = plt.figure(figsize=[15,8])

#create the display

myd = pyart.graph.RadarMapDisplay(radar)

#plot a single tilt, in this case the second tilt

myd.plot_ppi_map('reflectivity_horizontal', 1, vmin=-20, vmax=20,
                 min_lon=-157.1, max_lon=-156, min_lat=71.2, max_lat=71.6,
                 lon_lines = np.arange(-158, -154, .2),
                 lat_lines = np.arange(69, 72, .1), resolution = 'h',
                 auto_range=False)

#plot range rings at 10,20,30 and 40km

myd.plot_range_ring(10. * 1000.0, line_style = 'k-')
myd.plot_range_ring(20. * 1000.0, line_style = 'k--')
myd.plot_range_ring(30. * 1000.0, line_style = 'k-')
myd.plot_range_ring(40. * 1000.0, line_style = 'k--')

#the next two lines plots in radar centered coordinates
#and plots some cross hairs

myd.plot_line_xy([-40000.0, 40000.0], [0.0,0.0], line_style = 'k-')
myd.plot_line_xy([0.0,0.0],[-20000.0, 200000.0], line_style = 'k-')

#put a point where the radar is

myd.plot_point(radar.longitude['data'], radar.latitude['data'])

#square up the plot

ax1 = plt.gca()
x0,x1 = ax1.get_xlim()
y0,y1 = ax1.get_ylim()
ax1.set_aspect((x1-x0)/(y1-y0) - .05)

plt.show()
