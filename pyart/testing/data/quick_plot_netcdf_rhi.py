#! /usr/bin/env python
import pyart
import matplotlib.pyplot as plt

# plot quickly
pradar = pyart.io.read_netcdf('example_netcdf_rhi.nc')
display = pyart.graph.RadarDisplay(pradar)
fig = plt.figure()
ax = fig.add_subplot(111)
display.plot_rhi('reflectivity_horizontal', 0, vmin=-16, vmax=64)
fig.savefig('example_netcdf_rhi.png')
