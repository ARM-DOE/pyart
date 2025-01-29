"""
====================
Compare PPI vs CAPPI
====================

This example demonstrates how to create and compare PPI (Plan Position Indicator)
and CAPPI (Constant Altitude Plan Position Indicator) plots using radar data.

In this example, we load sample radar data, create a CAPPI at 2,000 meters
for the 'reflectivity' field, and then plot both the PPI and CAPPI for comparison.

"""

print(__doc__)

# Author: Hamid Ali Syed (syed44@purdue.edu)
# License: BSD 3 clause

import matplotlib.pyplot as plt
from open_radar_data import DATASETS

import pyart

# Load the sample radar data
file = DATASETS.fetch("RAW_NA_000_125_20080411190016")
radar = pyart.io.read(file)

# Apply gate filtering to exclude unwanted data
gatefilter = pyart.filters.GateFilter(radar)
gatefilter.exclude_transition()

# Create CAPPI at 2,000 meters for the 'reflectivity' field
cappi = pyart.retrieve.create_cappi(
    radar, fields=["reflectivity"], height=2000, gatefilter=gatefilter
)

# Create RadarMapDisplay objects for both PPI and CAPPI
radar_display = pyart.graph.RadarMapDisplay(radar)
cappi_display = pyart.graph.RadarMapDisplay(cappi)

# Plotting the PPI and CAPPI for comparison
fig, ax = plt.subplots(1, 2, figsize=(13, 5))

# Plot PPI for 'reflectivity' field
radar_display.plot_ppi("reflectivity", vmin=-10, vmax=60, ax=ax[0])
ax[0].set_title("PPI Reflectivity")

# Plot CAPPI for 'reflectivity' field
cappi_display.plot_ppi("reflectivity", vmin=-10, vmax=60, ax=ax[1])
ax[1].set_title("CAPPI Reflectivity at 2000 meters")

# Show the plots
plt.show()
