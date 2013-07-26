#! /usr/bin/env python

import netCDF4

dset = netCDF4.Dataset('sgpinterpolatedsondeC1.c1.20110510.000000.cdf')
dvars = dset.variables

out = netCDF4.Dataset('example_interpolatedsonde.cdf', mode='w')

out.createDimension('time', None)
out.createDimension('height', 316)

time = out.createVariable('time', dvars['time'].dtype, ('time'))
height = out.createVariable('height', dvars['height'].dtype, ('height'))
wspd = out.createVariable('wspd', dvars['wspd'].dtype, ('time', 'height'))
wdir = out.createVariable('wdir', dvars['wdir'].dtype, ('time', 'height'))

time[:] = dvars['time'][685:695]
time.units = dvars['time'].units
height[:] = dvars['height'][:]
wspd[:] = dvars['wspd'][685:695, :]
wdir[:] = dvars['wdir'][685:695, :]

out.close()
