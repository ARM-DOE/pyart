#! /usr/bin/env python
# Code to plot SUR file and search for an RHI to accompany it and plot it...

import sys
import getopt

import matplotlib as mpl
from pylab import *
import numpy as np
import netCDF4


def get_optargs(argv):
    try:
        opts, args = getopt.getopt(
            argv, "de:n:h:t:b:l:r:z:v:",
            ['debug', 'east', 'north', 'hieght', 'ylim1', 'ylim2', 'xlim1',
             'xlim2', 'zlim', 'var'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    return opts, args

if __name__ == "__main__":
    opts, args = get_optargs(sys.argv[1:])
    opt_dict = dict(opts)
    if ('--debug' in opt_dict.keys()) | ('-d' in opt_dict.keys()):
        debug = True
    else:
        debug = False
    east = int(opt_dict.get('-e', opt_dict.get('--east', '102')))
    north = int(opt_dict.get('-n', opt_dict.get('--north', '102')))
    height = int(opt_dict.get('-h', opt_dict.get('--height', '4')))
    y1 = float(opt_dict.get('-b', opt_dict.get('--ylim1', '-100')))
    y2 = float(opt_dict.get('-t', opt_dict.get('--ylim2', '100')))
    x1 = float(opt_dict.get('-l', opt_dict.get('--xlim1', '-100')))
    x2 = float(opt_dict.get('-r', opt_dict.get('--xlim2', '100')))
    zlim = float(opt_dict.get('-z', opt_dict.get('--zlim1', '17')))
    var = opt_dict.get('-v', opt_dict.get(
        '--var', 'attenuation_corrected_reflectivity_horizontal'))
    if debug:
        print "slicer.py "
    try:
        filename = args[0]
    except IndexError:
        usage()
        sys.exit(2)
    try:
        dirname = args[1]
    except IndexError:
        usage()
        sys.exit(2)
    print var
    ofile = netCDF4.Dataset(filename, 'r', format='NETCDF3_CLASSIC')
    times = ofile.variables['time']
    print "moo"
    dates = netCDF4.num2date(times[:], units=times.units)
    print "show"
    rad_head = dates[0]
    print "hi"
    print times[:].data
    ofname = dirname + filename.split('/')[-1] + '.png'
    print "foo"
    x = ofile.variables['x_disp'][:]
    y = ofile.variables['y_disp'][:]
    z = ofile.variables['z_disp'][:]
    f = figure(figsize=[14, 5])
    subplot(2, 2, 4)
    pcolormesh(x / 1000.0, z / 1000.0, ofile.variables[var][0, :, :, east],
               vmin=-4, vmax=72)
    xlim([x1, x2])
    ylim([0.0, zlim])
    cb = colorbar()
    cb.set_label('Eq Refl fac (dBz)')
    ylabel('Height ASL (km)')
    xlabel('Distance N of CF (km)')
    #title('Const Lon Slice')
    subplot(2, 2, 2)
    pcolormesh(x / 1000.0, z / 1000.0, ofile.variables[var][0, :, north, :],
               vmin=-4, vmax=72)
    xlim([y1, y2])
    ylim([0.0, zlim])
    cb = colorbar()
    cb.set_label('Eq Refl fac (dBz)')
    ylabel('Height ASL (km)')
    xlabel('Distance E of CF (km)')
    #title('Const Lat slice')
    subplot(1, 2, 1)
    pcolormesh(x / 1000.0, y / 1000.0, ofile.variables[var][0, height, :, :],
               vmin=-4, vmax=72)
    plot([x1, x2], [y[north] / 1000.0, y[north] / 1000.0], 'r--')
    plot([x[east] / 1000.0, x[east] / 1000.0], [y1, y2], 'r--')
    xlim([x1, x2])
    ylim([y1, y2])
    cb = colorbar()
    cb.set_label('Eq Refl fac (dBz)')
    ylabel('Distance N of CF (km)')
    xlabel('Distance E of CF (km)')
    title(dates[0])
    savefig(ofname)
    close(f)
    ofile.close()
