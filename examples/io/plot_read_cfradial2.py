"""
==========================================================
Read and Plot Cfradial2/FM301 data Using Xradar and Py-ART
==========================================================

An example which uses xradar and Py-ART to read and plot Cfradial2/FM301 data.

"""

# Author: Max Grover (mgrover@anl.gov)
# License: BSD 3 clause


import xarray as xr
from open_radar_data import DATASETS

import pyart

# Locate the test data and read in using xradar
filename = DATASETS.fetch("cfrad2.20080604_002217_000_SPOL_v36_SUR.nc")
tree = xr.open_datatree(filename)

# Give the tree Py-ART radar methods
radar = tree.pyart.to_radar()

# Plot the Reflectivity Field (corrected_reflectivity_horizontal)
display = pyart.graph.RadarMapDisplay(radar)
display.plot_ppi("DBZ", cmap="ChaseSpectral", vmin=-20, vmax=70)
