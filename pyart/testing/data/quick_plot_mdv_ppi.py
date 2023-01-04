#! /usr/bin/env python
import matplotlib.pyplot as plt

import pyart

# plot quickly
pradar = pyart.io.read_mdv("example_mdv_ppi.mdv")
display = pyart.graph.RadarDisplay(pradar)
fig = plt.figure()
ax = fig.add_subplot(111)
display.plot_ppi("reflectivity_horizontal", 0, vmin=-16, vmax=64)
fig.savefig("example_mdv_ppi.png")
