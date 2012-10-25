#!/usr/bin/python
import sys
import os
pyart_dir=os.environ.get('PYART_DIR',os.environ['HOME']+'/python')
sys.path.append(pyart_dir)
from pyart.io import py_mdv, radar
from time import time
from pyart.graph import rayplot, radar_display
from pyart.correct import phase_proc
import numpy as np
import pylab
from pyart.io import nc_utils
import netCDF4

if __name__ == "__main__":
	filename=sys.argv[1]
	outdir=sys.argv[2]
	my_mdv_object=py_mdv.read_mdv(filename, debug=True)
	myradar=radar.Radar(my_mdv_object)
	mylp=phase_proc.phase_proc(myradar, -2.0, debug=True)
	my_new_radar=mylp(debug=True)
	netcdf_obj=netCDF4.Dataset(outdir+'/testcase.nc', 'w',format='NETCDF4')
	nc_utils.write_radar4(netcdf_obj, my_new_radar)
	netcdf_obj.close()