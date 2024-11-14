"""
=========================================
Calculate and Plot VAD profile
=========================================

Calculates a VAD and plots a vertical profile of wind
"""

# Author: Daniel Wolfensberger (daniel.wolfensberger@meteoswiss.ch)
# License: BSD 3 clause

import matplotlib.pyplot as plt
import numpy as np
from open_radar_data import DATASETS

import pyart

# Read in a sample file
filename = DATASETS.fetch("MLA2119412050U.nc")
radar = pyart.io.read_cfradial(filename)

# Loop on all sweeps and compute VAD
zlevels = np.arange(100, 5000, 100)  # height above radar
u_allsweeps = []
v_allsweeps = []

for idx in range(radar.nsweeps):
    radar_1sweep = radar.extract_sweeps([idx])
    vad = pyart.retrieve.vad_browning(
        radar_1sweep, "corrected_velocity", z_want=zlevels
    )
    u_allsweeps.append(vad.u_wind)
    v_allsweeps.append(vad.v_wind)

# Average U and V over all sweeps and compute magnitude and angle
u_avg = np.nanmean(np.array(u_allsweeps), axis=0)
v_avg = np.nanmean(np.array(v_allsweeps), axis=0)
orientation = np.rad2deg(np.arctan2(-u_avg, -v_avg)) % 360
speed = np.sqrt(u_avg**2 + v_avg**2)

# Display vertical profile of wind
fig, ax = plt.subplots(1, 2, sharey=True)
ax[0].plot(speed * 2, zlevels + radar.altitude["data"])
ax[1].plot(orientation, zlevels + radar.altitude["data"])
ax[0].set_xlabel("Wind speed [m/s]")
ax[1].set_xlabel("Wind direction [deg]")
ax[0].set_ylabel("Altitude [m]")
fig.suptitle("Wind profile obtained from VAD")
