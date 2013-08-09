"""
====================================
Create a plot of NEXRAD reflectivity
====================================

An example which creates a plot containing the first collected super resolution
and standard resolution reflectivity scans from a NEXRAD file.

"""
print __doc__

# Author: Jonathan J. Helmus (jhelmus@anl.gov)
# License: BSD 3 clause

import matplotlib.pyplot as plt
import pyart

# open the file, create the displays and figure
filename = 'Level2_KATX_20130717_1950.ar2v'
radars = pyart.io.read_nexrad_archive(filename)
hires_display = pyart.graph.RadarDisplay(radars[0])
stdres_display = pyart.graph.RadarDisplay(radars[2])
fig = plt.figure(figsize=(12, 5))

# plot super resolution reflectivity
ax = fig.add_subplot(121)
hires_display.plot_ppi('reflectivity', 0, title='Super Resolution',
                       colorbar_label='', ax=ax)
hires_display.plot_range_ring(radars[0].range['data'][-1]/1000., ax=ax)
hires_display.set_limits(xlim=(-500, 500), ylim=(-500, 500), ax=ax)

# plot standard resolution reflectivity
ax = fig.add_subplot(122)
stdres_display.plot_ppi('reflectivity', 0, title='Standard Resolution',
                        colorbar_label='', ax=ax)
stdres_display.plot_range_ring(radars[2].range['data'][-1]/1000., ax=ax)
stdres_display.set_limits(xlim=(-500, 500), ylim=(-500, 500), ax=ax)

plt.show()
