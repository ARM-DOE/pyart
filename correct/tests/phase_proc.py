#! /usr/bin/env python
# creates the phase_processed_data.nc file needed for the unit tests
# this should be done using a version of py-ART which has passed the
# test_correct.py:test_phase_rsl test.

import numpy as np
import netCDF4
import pyart
from pyart.correct import phase_proc

# read in the data
FILENAME = "XSE110510113001.RAW7BPA"
radarobj = pyart.io.py4dd.RSL_anyformat_to_radar(FILENAME)
radar = pyart.io.radar.Radar(radarobj)

# process phase
gates = radar.range['data'][1] - radar.range['data'][0]
rge = 10.0 * 1000.0
ng = rge / gates
mylp = phase_proc.phase_proc(radar, 8.6, 
                             sys_phase=332.0,
                             overide_sys_phase=True, 
                             debug=True, nowrap=ng)
reproc_phase, sob_kdp = mylp(debug=True)

radar.fields.update({'recalculated_diff_phase': sob_kdp,
                          'proc_dp_phase_shift': reproc_phase})

ofilename = "phase_processed_data.nc"
netcdf_obj = netCDF4.Dataset(ofilename, 'w', format='NETCDF4')
pyart.io.nc_utils.write_radar4(netcdf_obj, radar)
netcdf_obj.close()
