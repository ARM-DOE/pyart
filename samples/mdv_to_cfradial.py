#!/usr/bin/python
#Sample for processing phase using LP methods
import matplotlib
matplotlib.use('agg')

import sys
import os
pyart_dir=os.environ.get('PYART_DIR',os.environ['HOME']+'/python')
sys.path.append(pyart_dir)
from pyart.io import py_mdv, radar
from time import time
from pyart.graph import rayplot, radar_display
import numpy as np
import pylab
from pyart.io import nc_utils
import netCDF4

def dt_to_dict(dt, **kwargs):
        pref=kwargs.get('pref', '')
        return dict( [(pref+key, getattr(dt, key)) for key in    ['year', 'month', 'day', 'hour', 'minute', 'second']])


if __name__ == "__main__":
	filename=sys.argv[1]
	outdir=sys.argv[2]
	site_deseg=sys.argv[3]
	my_mdv_object=py_mdv.read_mdv(filename, debug=True)
	myradar=radar.Radar(my_mdv_object)
	mydatetime=netCDF4.num2date(myradar.time['data'][0], myradar.time['units'], calendar=myradar.time['calendar']) #append a datetime object
	mydict=dt_to_dict(mydatetime)
	mydict.update({'scanmode':{'ppi':'sur','rhi':'rhi'}[myradar.sweep_mode[0]], 'site_deseg':site_deseg})
	ofilename=outdir+'%(scanmode)scmac%(site_deseg)s.a1.%(year)04d%(month)02d%(day)02d.%(hour)02d%(minute)02d%(second)02d.nc' % mydict
	netcdf_obj=netCDF4.Dataset(ofilename, 'w',format='NETCDF4')
	nc_utils.write_radar4(netcdf_obj, myradar)
	netcdf_obj.close()