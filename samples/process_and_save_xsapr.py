#! /usr/bin/env python
# Sample script for processing phase using LP methods
# usage : process_and_save_xsapr.py filename outdir offset is_dir fac

import sys
import copy

import netCDF4

from pyart.io import py4dd, radar, nc_utils
from pyart.correct import phase_proc, attenuation, dealias


def dt_to_dict(dt, **kwargs):
        pref = kwargs.get('pref', '')
        return dict([(pref+key, getattr(dt, key)) for key in
                    ['year', 'month', 'day', 'hour', 'minute', 'second']])


if __name__ == "__main__":

    # read the command line arguments
    filename = sys.argv[1]
    outdir = sys.argv[2]
    offset = float(sys.argv[3])
    is_dir = sys.argv[4]
    fac = sys.argv[5]

    # read in the data
    myradar = radar.Radar(py4dd.RSL_anyformat_to_radar(filename))

    # dealias
    target = netCDF4.num2date(myradar.time['data'][0], myradar.time['units'])
    fname = 'sgpinterpolatedsondeC1.c1.%(year)04d%(month)02d%(day)02d.000000.cdf' % radar.dt_to_dict(target)
    print fname
    interp_sonde = netCDF4.Dataset(is_dir + fname)
    myheight, myspeed, mydirection = dealias.find_time_in_interp_sonde(
        interp_sonde, target)

    my_new_field = dealias.dealias(myradar, myheight * 1000.0, myspeed,
                                   mydirection, target)
    myradar.fields.update({'corrected_mean_doppler_velocity': my_new_field})
    interp_sonde.close()

    # process phase
    gates = myradar.range['data'][1] - myradar.range['data'][0]
    rge = 10.0 * 1000.0
    ng = rge / gates
    mydatetime = netCDF4.num2date(
        myradar.time['data'][0], myradar.time['units'],
        calendar=myradar.time['calendar'])  # append a datetime object
    mydict = dt_to_dict(mydatetime)
    mydict.update({'scanmode': {'ppi': 'sur', 'rhi': 'rhi'}
                  [myradar.sweep_mode[0]], 'fac': fac})
    ofilename = outdir + '%(scanmode)scmac%(fac)s.c0.%(year)04d%(month)02d%(day)02d.%(hour)02d%(minute)02d%(second)02d.nc' % mydict
    mylp = phase_proc.phase_proc(myradar, offset, debug=True, nowrap=ng)
    reproc_phase, sob_kdp = mylp(debug=True)
    myradar.fields.update({'recalculated_diff_phase': sob_kdp,
                           'proc_dp_phase_shift': reproc_phase})

    # attenuation correction
    spec_at, cor_z = attenuation.calculate_attenuation(
        myradar, offset, debug=True, a_coef=0.17)
    myradar.fields.update({'specific_attenuation': spec_at})
    myradar.fields.update({'corrected_reflectivity_horizontal': cor_z})
    R = 51.3 * (myradar.fields['specific_attenuation']['data']) ** 0.81
    rainrate = copy.deepcopy(myradar.fields['diff_phase'])
    rainrate['data'] = R
    rainrate['valid_min'] = 0.0
    rainrate['valid_max'] = 400.0
    rainrate['standard_name'] = 'rainfall_rate'
    rainrate['long_name'] = 'rainfall_rate'
    rainrate['least_significant_digit'] = 1
    rainrate['units'] = 'mm/hr'
    myradar.fields.update({'rain_rate_A': rainrate})
    netcdf_obj = netCDF4.Dataset(ofilename, 'w', format='NETCDF4')
    nc_utils.write_radar4(netcdf_obj, myradar)
    netcdf_obj.close()
