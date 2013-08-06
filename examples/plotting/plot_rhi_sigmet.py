"""
====================================
Create a RHI plot from a Sigmet file
====================================

An example which creates a RHI plot of a Sigmet file using both a RslDisplay
and a RadarDisplay object.

"""
print(__doc__)

# Author: Jonathan J. Helmus (jhelmus@anl.gov)
# License: BSD 3 clause

import matplotlib.pyplot as plt
import pyart

filename = 'XSW110520113537.RAW7HHL'

# create the plot using RslDisplay (this method is not recommended)
rslradar = pyart.io._rsl_interface.RslFile(filename)
display = pyart.graph.RslDisplay(rslradar)

fig = plt.figure(figsize=[10, 4])
ax = fig.add_subplot(111)

radar_name = rslradar.get_radar_header()['name']
time_text = ' ' + display.time_begin.isoformat() + 'Z '
azimuth = rslradar.get_volume(4).get_sweep_azimuths()[0]
title = 'RHI ' + radar_name + time_text + 'Azimuth %.2f' % (azimuth)

display.plot_rhi('ZT', 0, vmin=-32, vmax=64, title=title,
                 colorbar_flag=False, ax=ax)
display.set_limits(ylim=[0, 17])

cax = fig.add_axes([.9, .1, 0.02, .8])
colorbar_label = 'Eq refl fact (dBz)'
display.plot_colorbar(fig=fig, cax=cax, label=colorbar_label)

plt.show()

# create the plot using RadarDisplay (recommended method)
radar = pyart.io.read_rsl(filename)
display = pyart.graph.RadarDisplay(radar)

fig = plt.figure(figsize=[10, 4])
ax = fig.add_subplot(111)

radar_name = radar.metadata['instrument_name']
time_text = ' ' + display.time_begin.isoformat() + 'Z '
azimuth = radar.fixed_angle['data'][0]
title = 'RHI ' + radar_name + time_text + 'Azimuth %.2f' % (azimuth)

display.plot_rhi('reflectivity_horizontal', 0, vmin=-32, vmax=64,
                 title=title, colorbar_flag=False, ax=ax)
display.set_limits(ylim=[0, 17])

cax = fig.add_axes([.9, .1, 0.02, .8])
colorbar_label = 'Eq refl fact (dBz)'
display.plot_colorbar(fig=fig, cax=cax, label=colorbar_label)

plt.show()
