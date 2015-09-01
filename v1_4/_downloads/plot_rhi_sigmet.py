"""
====================================
Create a RHI plot from a Sigmet file
====================================

An example which creates a RHI plot of a Sigmet file using a RadarDisplay
object.

"""
print(__doc__)

# Author: Jonathan J. Helmus (jhelmus@anl.gov)
# License: BSD 3 clause

import matplotlib.pyplot as plt
import pyart

filename = 'XSW110520113537.RAW7HHL'

# create the plot using RadarDisplay
radar = pyart.io.read_rsl(filename)
display = pyart.graph.RadarDisplay(radar)

fig = plt.figure(figsize=[10, 4])
ax = fig.add_subplot(111)

radar_name = radar.metadata['instrument_name']
time_text = ' ' + display.time_begin.isoformat() + 'Z '
azimuth = radar.fixed_angle['data'][0]
title = 'RHI ' + radar_name + time_text + 'Azimuth %.2f' % (azimuth)

display.plot('reflectivity', 0, vmin=-32, vmax=64,
             title=title, colorbar_flag=False, ax=ax)
display.set_limits(ylim=[0, 17])

cax = fig.add_axes([.9, .1, 0.02, .8])
colorbar_label = 'Eq refl fact (dBz)'
display.plot_colorbar(fig=fig, cax=cax, label=colorbar_label)

plt.show()
