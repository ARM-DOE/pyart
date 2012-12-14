#!/usr/bin/python
#Sample for processing phase using LP methods
import matplotlib
matplotlib.use('agg')

import sys
import os
pyart_dir=os.environ.get('PYART_DIR',os.environ['HOME']+'/python')
sys.path.append(pyart_dir)
from pyart.io import py4dd, radar
from time import time
from pyart.graph import rayplot, radar_display
from pyart.correct import phase_proc
import numpy as np
import pylab
from pyart.io import nc_utils
import netCDF4
from pyart.correct import attenuation
from pyart.correct import dealias


def dt_to_dict(dt, **kwargs):
        pref=kwargs.get('pref', '')
        return dict( [(pref+key, getattr(dt, key)) for key in    ['year', 'month', 'day', 'hour', 'minute', 'second']])


if __name__ == "__main__":
	filename=sys.argv[1]
	outdir=sys.argv[2]
	offset=float(sys.argv[3])
	myradar=radar.Radar(py4dd.RSL_anyformat_to_radar(filename))
	#dealias
	is_dir=sys.argv[4]
	target=netCDF4.num2date(myradar.time['data'][0], myradar.time['units'])
	fname='sgpinterpolatedsondeC1.c1.%(year)04d%(month)02d%(day)02d.000000.cdf' %radar.dt_to_dict(target)
	print fname
	interp_sonde=netCDF4.Dataset(is_dir+fname)
	myheight,myspeed,mydirection=dealias.find_time_in_interp_sonde(interp_sonde, target)
	deal_obj=dealias.dealiaser(myradar, myheight*1000.0,myspeed,mydirection, target)
	my_new_field=deal_obj()
	myradar.fields.update({'corrected_mean_doppler_velocity':my_new_field})
	interp_sonde.close()
	#Process Phase
	gates=myradar.range['data'][1]-myradar.range['data'][0]
	rge=5.0*1000.0
	ng=rge/gates
	mydatetime=netCDF4.num2date(myradar.time['data'][0], myradar.time['units'], calendar=myradar.time['calendar']) #append a datetime object
	mydict=dt_to_dict(mydatetime)
	mydict.update({'scanmode':{'ppi':'sur','rhi':'rhi'}[myradar.sweep_mode[0]], 'fac':sys.argv[5]})
	ofilename=outdir+'%(scanmode)scmac%(fac)s.c0.%(year)04d%(month)02d%(day)02d.%(hour)02d%(minute)02d%(second)02d.nc' % mydict
	mylp=phase_proc.phase_proc(myradar, offset, debug=True, nowrap=ng)
	reproc_phase, sob_kdp=mylp(debug=True)
	myradar.fields.update({'recalculated_diff_phase':sob_kdp, 'proc_dp_phase_shift': reproc_phase})
	#attenuation correction
	spec_at, cor_z=attenuation.calculate_attenuation(myradar,offset, debug=True, a_coef=0.17)
	myradar.fields.update({'specific_attenuation':spec_at})
	myradar.fields.update({'corrected_reflectivity_horizontal':cor_z})
	netcdf_obj=netCDF4.Dataset(ofilename, 'w',format='NETCDF4')
	nc_utils.write_radar4(netcdf_obj, myradar)
	netcdf_obj.close()