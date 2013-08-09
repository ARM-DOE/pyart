"""
===================================
Linear programming phase processing
===================================

An example of using linear processing to process the differential phase
fields of a ARM C-SAPR radar.

"""
print __doc__

# Author: Jonathan J. Helmus (jhelmus@anl.gov)
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt
import pyart


# perform LP phase processing (this takes a while)
radar = pyart.io.read_mdv('095636.mdv')

# the next line force only the first sweep to be processed, this
# significantly speeds up the calculation but should be commented out
# in production so that the entire volume is processed
radar.sweep_start_ray_index['data'] = np.array([0])

phidp, kdp = pyart.correct.phase_proc_lp(radar, 0.0, debug=True)
radar.fields['proc_dp_phase_shift'] = phidp
radar.fields['recalculated_diff_phase'] = kdp

# the following line can be used to save/load in preprocessed data
#pyart.io.write_netcdf('preprocessed.nc', radar)
#radar = pyart.io.read_netcdf('preprocessed.nc')

# create a plot of the various differential phase fields
display = pyart.graph.RadarDisplay(radar)
fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(221)
display.plot_ppi('dp_phase_shift', 0, ax=ax1,
                 title='Raw Differential Phase', colorbar_label='',
                 axislabels_flag=False)

ax2 = fig.add_subplot(222)
display.plot_ppi('proc_dp_phase_shift', 0, ax=ax2,
                 title='Processed Differential Phase', colorbar_label='',
                 axislabels_flag=False)

ax3 = fig.add_subplot(223)
display.plot_ppi('diff_phase', 0, ax=ax3,
                 title='Raw Specific Differential Phase', colorbar_label='',
                 axislabels_flag=False)

ax4 = fig.add_subplot(224)
display.plot_ppi('recalculated_diff_phase', 0, ax=ax4,
                 title='Processed Specific Differential Phase',
                 colorbar_label='',
                 axislabels_flag=False)
plt.show()

# plot a fields from a single ray
fig = plt.figure(figsize=[10, 5])
ax = fig.add_subplot(111)

ray_num = 191
range_km = radar.range['data'] / 1000.
phidp_ray = radar.fields['proc_dp_phase_shift']['data'][ray_num]
unfolded_phidp_ray = radar.fields['unf_dp_phase_shift']['data'][ray_num]
kdp_ray = radar.fields['recalculated_diff_phase']['data'][ray_num]
refl_ray = radar.fields['reflectivity_horizontal']['data'][ray_num]

# filtered phidp and unfolded phidp
p1, = ax.plot(range_km, phidp_ray, 'b-')
p2, = ax.plot(range_km, unfolded_phidp_ray, 'g-')

# set labels
ax.set_ylim(0, 250)
ax.set_ylabel('Differential phase shift (degrees)')
ax.set_xlabel('Range (km)')

# plot KDP and reflectivity on second axis
ax2 = ax.twinx()
p3, = ax2.plot(range_km, kdp_ray, 'r-')
p4, = ax2.plot(range_km, refl_ray/10.)

# decorate
ax2.yaxis.grid(color='gray', linestyle='dashed')
ax.legend([p1, p2, p3, p4],
          ["Filtered phiDP", "Unfolded phiDP", 'KDP', 'Z/10.0'],
          loc='upper left')
plt.show()
