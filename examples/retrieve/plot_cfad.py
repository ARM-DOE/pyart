"""
=======================================
Creating a CFAD diagram
=======================================
This example shows how to create a contoured frequency by altitude (CFAD) diagram
"""

print(__doc__)

# Author: Laura Tomkins (lmtomkin@ncsu.edu)
# License: BSD 3 clause


import matplotlib.pyplot as plt
import numpy as np

import pyart

######################################
# Example 1

# get test data
filename = pyart.testing.get_test_data("sgpxsaprrhicmacI5.c0.20110524.015604_NC4.nc")
radar = pyart.io.read_cfradial(filename)

# compute CFAD
# get reflectivity data and mask extremes
subset_slice = radar.get_slice(0)
ref_data = radar.fields['reflectivity_horizontal']['data'][subset_slice]
ref_data_masked = np.ma.masked_outside(ref_data, -15, 35)
# get altitude data
_, _, gate_z = radar.get_gate_x_y_z(0)
gate_z_masked = np.ma.masked_where(ref_data_masked.mask, gate_z)

freq_norm, height_edges, field_edges = createCFAD(ref_data_masked, gate_z_masked, field_bins=np.linspace(-15,35,100),
                                                  altitude_bins=np.linspace(0,15000,100), min_frac_thres=0.1)

# plot CFAD
freq_norm_masked = np.ma.masked_less_equal(freq_norm, 0)

fig = plt.figure()
ax = plt.axes()
cfad_pm = ax.pcolormesh(field_edges, height_edges, freq_norm_masked, cmap='plasma', vmin=0, vmax=0.10)
plt.colorbar(cfad_pm)
ax.set_xlabel('Reflectivity [dBZ]')
ax.set_ylabel('Height [m]')
ax.grid(ls='--', color='gray', lw=0.5, alpha=0.7)
plt.show()

# plot RHI data
display = pyart.graph.RadarDisplay(radar)

fig = plt.figure(figsize=[12,3])
ax = plt.axes()
plt.tight_layout()
display.plot('reflectivity_horizontal', 0, vmin=-15, vmax=35, mask_outside=True, cmap='magma_r', ax=ax)
display.set_limits(ylim=[0,15], ax=ax)
ax.set_aspect('equal')
plt.show()

######################################
# Example 2

# get test data
filename = pyart.testing.get_test_data("034142.mdv")
radar = pyart.io.read_mdv(filename)

# compute CFAD
# get reflectivity data and mask extremes
subset_slice = radar.get_slice(0)
ref_data = radar.fields['reflectivity']['data'][subset_slice]
ref_data_masked = np.ma.masked_outside(ref_data, -5, 60)
# get altitude data
_, _, gate_z = radar.get_gate_x_y_z(0)
gate_z_masked = np.ma.masked_where(ref_data_masked.mask, gate_z)

freq_norm, height_edges, field_edges = createCFAD(ref_data_masked, gate_z_masked, field_bins=np.linspace(-5,60,100),
                                                  altitude_bins=np.linspace(0,20000,100), min_frac_thres=0.1)

# plot CFAD
freq_norm_masked = np.ma.masked_less_equal(freq_norm, 0)

fig = plt.figure()
ax = plt.axes()
cfad_pm = ax.pcolormesh(field_edges, height_edges, freq_norm_masked, cmap='plasma', vmin=0, vmax=0.10)
plt.colorbar(cfad_pm)
ax.set_xlabel('Reflectivity [dBZ]')
ax.set_ylabel('Height [m]')
ax.grid(ls='--', color='gray', lw=0.5, alpha=0.7)
plt.show()

# plot RHI data
display = pyart.graph.RadarDisplay(radar)

fig = plt.figure(figsize=[12,3])
ax = plt.axes()
plt.tight_layout()
display.plot('reflectivity', 0, vmin=-5, vmax=60, cmap='magma_r', ax=ax)
display.set_limits(ylim=[0,20], xlim=[0,80], ax=ax)
ax.set_aspect('equal')
plt.show()