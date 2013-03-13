#!/usr/bin/env python
# Code to plot a Domer RHI file

import sys, os
import getopt

import pyart.io._rsl as _rsl
from pylab import *
import numpy as N

from pyart.graph.common import radar_coords_to_cart

def create_RHI_array(sweep):
    ppi = N.zeros([sweep.h.nrays, sweep.rays[0].h.nbins],
                  dtype=float) + 1.31072000e+05
    for raynum in range(sweep.h.nrays):
        data = sweep.rays[raynum].data
        ppi[raynum, 0:len(data)] = sweep.rays[raynum].data
    return ppi


def get_optargs(argv):
    try:
        opts, args = getopt.getopt(argv, "d", ['debug'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    return opts, args


def usage():
    print """
    plot_rhi.py -d filename outpath
    should work on UF, Sigmet/IRIS/88D (not that they do RHIs!)
    """


def plot_rhi_multi(xsapr, imagefilename, varss, sweep, **kwargs):
    debug = kwargs.get('debug', False)
    rges = kwargs.get('rges', {
        'CZ': [-22., 64.],
        'ZT': [-16.0, 64.0],
        'DZ': [-16.0, 64.0],
        'VR': [-16, 16],
        'RH': [0.9, 1.0],
        'KD': [-1, 5],
        'ZD': [-1, 6],
        'DR': [-1, 6],
        'PH': [100, 150]})
    var_lab = kwargs.get('var_lab', {
        'CZ': 'Eq refl fact (dBz)',
        'ZT': 'Eq refl fact (dBz)',
        'DZ': 'Eq refl fact (dBz)',
        'VR': 'Radial Velocity (+away) m/s',
        'RH': 'C. Coef. (frac)',
        'ZD': 'Dif Refl (dB)',
        'KD': 'Spec. Dif. Phase (deg/km)',
        'DR': 'Dif Refl (dB)',
        'PH': 'Phidp degs'})
    fsize = kwargs.get('fsize', [10, 10])
    print rges
    fields = [_rsl.fieldTypes().list.index(var) for var in varss]
    plain_data = [create_RHI_array(xsapr.contents.volumes[field].sweeps[sweep])
                  for field in fields]
    if 'mask' in kwargs.keys():
        plain_mask = create_RHI_array(xsapr.contents.volumes[
            _rsl.fieldTypes().list.index(kwargs['mask'][0])].sweeps[sweep])
    azmths = array(
        [xsapr.contents.volumes[fields[0]].sweeps[sweep].rays[i].h.azimuth
         for i in
         range(xsapr.contents.volumes[fields[0]].sweeps[sweep].h.nrays)])
    ranges = (xsapr.contents.volumes[fields[0]].sweeps[sweep].rays[0].dists /
              1000.0)
    elevs = array(
        [xsapr.contents.volumes[fields[0]].sweeps[sweep].rays[i].h.elev
         for i in
         range(xsapr.contents.volumes[fields[0]].sweeps[sweep].h.nrays)])
    rg, ele = meshgrid(ranges, elevs)
    azg = (ones(rg.shape, dtype=float) *
           xsapr.contents.volumes[fields[0]].sweeps[sweep].h.azimuth)
    x, y, z = radar_coords_to_cart(rg, azg, ele)
    titl = kwargs.get('titl',
                      "RHI %(radar_name)s %(year)04d-%(month)02d-%(day)02d %(hour)02d:%(min)02d:%(sec)02d Azimuth: %(az).2f")
    rad_head = xsapr.contents.h
    titl_comp = titl % {
        'radar_name': rad_head.radar_name,
        'year': rad_head.year,
        'month': rad_head.month,
        'day': rad_head.day,
        'hour': rad_head.hour,
        'min': rad_head.minute,
        'sec': rad_head.sec,
        'az': azmths[0]}
    f = figure(figsize=fsize)
    for ff in range(len(varss)):
        if 'bias' in kwargs.keys():
            if varss[ff] in kwargs['bias']:
                bias = kwargs['bias'][varss[ff]]
            else:
                bias = 0
        else:
            bias = 0
        print varss[ff], ' BIAS:', bias
        subplot(len(varss), 1, ff + 1)
        if 'mask' in kwargs.keys():
            pp = pcolormesh(
                sign(y) * sqrt(y ** 2 + x ** 2) / 1000.0,
                z / 1000.0,
                ma.masked_where(logical_or(plain_mask < kwargs['mask'][1],
                                plain_mask > kwargs['mask'][2]),
                                plain_data[ff]) +
                bias, vmin=rges[varss[ff]][0], vmax=rges[varss[ff]][1])
        else:
            pp = pcolormesh(
                sign(y) * sqrt(y ** 2 + x ** 2) / 1000.0,
                z / 1000.0,
                N.ma.masked_outside(plain_data[ff],
                                    rges[varss[ff]][0],
                                    rges[varss[ff]][1]) +
                bias, vmin=rges[varss[ff]][0], vmax=rges[varss[ff]][1])
        ylim(kwargs.get('ylim', [0, 17]))
        ylabel('Height above radar (km)')
        if ff == 0:
            title(titl_comp)
        cb = colorbar()  # (cax=cbax, mappable=pp)
        cb.set_label(var_lab[varss[ff]])
        if ff == len(varss) - 1:
            xlabel('Distance from radar (km)')
        else:
            gca().xaxis.set_ticklabels([])
    subplots_adjust(right=0.75, hspace=0.051)
    savefig(imagefilename)
    close(f)


if __name__ == "__main__":
    pwd = os.path.abspath(os.path.curdir) + '/'
    opts, args = get_optargs(sys.argv[1:])
    opt_dict = dict(opts)
    if ('--debug' in opt_dict.keys()) | ('-d' in opt_dict.keys()):
        debug = True
    else:
        debug = False
    if debug:
        print "plot_rhi.py "
    try:
        filename = args[0]
    except IndexError:
        usage()
        sys.exit(2)
    if debug:
        print "file name=", filename
    if debug:
        _rsl.RSL_radar_verbose_on()
    xsapr = _rsl.RSL_anyformat_to_radar(filename)
    plot_rhi_multi(xsapr, 'test_rhi0.png',
                   ['DZ'], 0, titl="RHI " + filename.split('/')[-1][0:3] +
                   "-%(radar_name)sp %(year)04d-%(month)02d-%(day)02d %(hour)02d:%(min)02d:%(sec)02d Azimuth: %(az).2f",
                   rges={
                       'VR': [-17, 17],
                       'CZ': [-22., 64.],
                       'ZT': [-16.0, 64.0],
                       'DZ': [-22.0, 34.0],
                       'PH': [90, 160],
                       'RH': [0.9, 1.0],
                       'KD': [-1, 5],
                       'DR': [-6, 6]},
                   fsize=[10, 5])
    print "hi"
    _rsl.RSL_free_radar(xsapr)
