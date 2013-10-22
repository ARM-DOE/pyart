#! /usr/bin/env python
import pyart
import matplotlib.pyplot as plt

# plot quickly
pradar = pyart.io.read_cfradial('example_cfradial_ppi.nc')
display = pyart.graph.RadarDisplay(pradar)
fig = plt.figure()
ax = fig.add_subplot(111)
display.plot_ppi('reflectivity_horizontal', 0, vmin=-16, vmax=64)
fig.savefig('example_cfradial_ppi.png')
