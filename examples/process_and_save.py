#! /usr/bin/env python
# Sample script for processing phase using LP methods
# usage : process_and_save.py filename.mdv outdir

import sys
import os
from time import time

import netCDF4
import numpy as np
import pylab
import matplotlib

from pyart.io import py_mdv, radar
from pyart.correct import phase_proc
from pyart.io import nc_utils
from pyart.correct import attenuation


def dt_to_dict(dt, **kwargs):
        pref = kwargs.get('pref', '')
        return dict([(pref+key, getattr(dt, key)) for key in
                    ['year', 'month', 'day', 'hour', 'minute', 'second']])


if __name__ == "__main__":

    # read in the command line arguments
    filename = sys.argv[1]
    outdir = sys.argv[2]

    # open the radar file
    my_mdv_object = py_mdv.read_mdv(filename, debug=True)
    myradar = radar.Radar(my_mdv_object)

    gates = myradar.range['data'][1] - myradar.range['data'][0]
    rge = 10.0 * 1000.0
    ng = rge / gates
    #append a datetime object
    mydatetime = netCDF4.num2date(myradar.time['data'][0],
                                  myradar.time['units'],
                                  calendar=myradar.time['calendar'])
    mydict = dt_to_dict(mydatetime)
    mydict.update({'scanmode': {'ppi': 'sur', 'rhi': 'rhi'}
                   [myradar.sweep_mode[0]]})

    ofilename = outdir + '%(scanmode)scmacI7.c0.%(year)04d%(month)02d%(day)02d.%(hour)02d%(minute)02d%(second)02d.nc' % mydict

    reproc_phase, sob_kdp = phase_proc.phase_proc(myradar, -2.0, debug=True, 
                                                  nowrap=ng)
    myradar.fields.update({'recalculated_diff_phase': sob_kdp,
                           'proc_dp_phase_shift': reproc_phase})
    spec_at, cor_z = attenuation.calculate_attenuation(
        myradar, -2.0, debug=True)
    myradar.fields.update({'specific_attenuation': spec_at})
    myradar.fields.update({'corrected_reflectivity_horizontal': cor_z})
    netcdf_obj = netCDF4.Dataset(ofilename, 'w', format='NETCDF4')
    nc_utils.write_radar4(netcdf_obj, myradar)
    netcdf_obj.close()
