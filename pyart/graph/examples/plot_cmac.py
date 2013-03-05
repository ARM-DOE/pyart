#!/usr/bin/env python
"""
Plotting utilities for working with the CMAC VAP

"""

import netCDF4
from pylab import *
from numpy import *


def dt_to_dict(dt, **kwargs):
    pref = kwargs.get('pref', '')
    return dict([(pref + key, getattr(dt, key)) for key in
                ['year', 'month', 'day', 'hour', 'minute', 'second']])


def radar_coords_to_cart(rng, az, ele, debug=False):
    """
    Asumes standard atmosphere, ie R = 4 * Re / 3
    """
    Re = 6371.0 * 1000.0
    p_r = 4.0 * Re / 3.0
    rm = rng * 1000.0
    z = (rm ** 2 + p_r ** 2 + 2.0 * rm * p_r *
         sin(ele * pi / 180.0)) ** 0.5 - p_r
    #arc length
    s = p_r * arcsin(rm * cos(ele * pi / 180.) / (p_r + z))
    if debug:
        print "Z=", z, "s=", s
    y = s * cos(az * pi / 180.0)
    x = s * sin(az * pi / 180.0)
    return x, y, z


def plot_multi_hsrhi(ncobj, var, **kwargs):
    ncvars = ncobj.variables
    sweeps = ncvars['sweep_number'][:]
    myfig = figure(figsize=[12, 17])
    myfig.subplots_adjust(hspace=0.4)
    for snum in sweeps:
        i1 = ncvars['sweep_start_ray_index'][snum]
        i2 = ncvars['sweep_end_ray_index'][snum]
        print snum
        print ncvars['range'].shape, ncvars['elevation'].shape
        rg, eleg = meshgrid(ncvars['range'][:] / 1000.0,
                            ncvars['elevation'][i1:i2])
        az = ones(rg.shape, dtype=float32) * ncvars['fixed_angle'][snum]
        x, y, z = radar_coords_to_cart(rg, az, eleg)
        subplot(len(sweeps), 1, snum + 1)
        pcolormesh(sqrt(x ** 2 + y ** 2) * sign(y) / 1000.0, z / 1000.0,
                   mynetcdf.variables[var][i1:i2, :], vmin=-20, vmax=20)
        ylim([0, 15])
        xlabel('Distance from radar (km)')
        ylabel('Height agl (km)')
        tt = 'HSRHI Az=%(az)f.3 %(extra)s' % {
            'az': ncvars['fixed_angle'][snum], 'extra': kwargs['extra']}
        title(tt)
        cb = colorbar()
        cb.set_label('Hz. Eq. Refl. Fac. (dBZ)')
    return myfig


if __name__ == "__main__":
    filename = sys.argv[1]
    outdir = sys.argv[2]
    mynetcdf = netCDF4.Dataset(filename, 'r')
    mydate = netCDF4.num2date(mynetcdf.variables['time'][:],
                              mynetcdf.variables['time'].units)
    fname = 'multi_panel_hsrhi_%(year)04d%(month)02d%(day)02d.%(hour)02d%(minute)02d%(second)02d.png' % dt_to_dict(mydate[0])
    myf = plot_multi_hsrhi(mynetcdf, 'reflectivity_horizontal',
                           extra=sys.argv[3])
    ttx = 'Time: %(year)04d%(month)02d%(day)02d %(hour)02d%(minute)02d%(second)02d Z' % dt_to_dict(mydate[0])
    myf.text(0.35, 0.92, ttx)
    savefig(outdir + '/' + sys.argv[3] + fname)
    close(myf)
    mynetcdf.close()
