#! /usr/bin/env python
# create plot of
# usage : plot_cf.py filename.nc

import sys

import netCDF4
import matplotlib as mpl
from pylab import *
from numpy import sqrt


def dt_to_dict(dt, **kwargs):
    pref = kwargs.get('pref', '')
    return dict([(pref + key, getattr(dt, key)) for key in
                ['year', 'month', 'day', 'hour', 'minute', 'second']])


def radar_coords_to_cart(rng, az, ele, debug=False):
    """
    Asumes standard atmosphere, ie R=4Re/3
    """
    Re = 6371.0 * 1000.0
    #h=(r^2 + (4Re/3)^2 + 2r(4Re/3)sin(ele))^1/2 -4Re/3
    #s=4Re/3arcsin(rcos(ele)/(4Re/3+h))
    p_r = 4.0 * Re / 3.0
    rm = rng * 1000.0
    z = (rm ** 2 + p_r ** 2 + 2.0 * rm * p_r *
         np.sin(ele * np.pi / 180.0)) ** 0.5 - p_r
    #arc length
    s = p_r * np.arcsin(rm * np.cos(ele * np.pi / 180.) / (p_r + z))
    if debug:
        print "Z=", z, "s=", s
    y = s * np.cos(az * pi / 180.0)
    x = s * np.sin(az * pi / 180.0)
    return x, y, z


mynetcdf = netCDF4.Dataset(sys.argv[1], 'r')
mydict = dt_to_dict(netCDF4.num2date(mynetcdf.variables['time'][0],
                    units=mynetcdf.variables['time'].units,
                    calendar=mynetcdf.variables['time'].calendar))
s = 0
i1 = mynetcdf.variables['sweep_start_ray_index'][s]
i2 = mynetcdf.variables['sweep_end_ray_index'][s]
print i2
if i2 > 200:
    rg, ele = meshgrid(mynetcdf.variables['range'],
                       mynetcdf.variables['elevation'][i1:i2])
    az = ones(rg.shape, dtype=float32) * mynetcdf.variables['fixed_angle'][s]
    x, y, z = radar_coords_to_cart(rg, az, ele)
    f = figure()
    pc1 = pcolormesh(sqrt(x ** 2 + y ** 2) / 1000.0, z / 1000.0,
                     mynetcdf.variables['reflectivity_horizontal'][i1:i2, :],
                     vmax=64, vmin=-8)
    title("CSAPR Manus RHI %(year)04d-%(month)02d-%(day)02d %(hour)02d:%(minute)02d" % mydict)
    gca().set_ylim([0, 17])
    ax3 = f.add_axes([.9, .1, 0.02, .8])
    cb = colorbar(cax=ax3, mappable=pc1)
    ylabel('Eq. refl. fac. (dBz)')
    fname = 'test_cmac_%(year)04d%(month)02d%(day)02d%(hour)02d%(minute)02d.png' % mydict
    savefig(fname)
    print "yobbo"
    close(f)
mynetcdf.close()
