"""
pyart.graph.plot_rsl
====================

Routines for plotting radar data from files readable by RSL.

.. autosummary::
    :toctree:: generated/


"""

#!/usr/bin/python
# Code to plot a PPI from p4dd

import sys
import getopt

import matplotlib
from pylab import *
import numpy as N

import pyart.io._rsl as _rsl
from .common import ax_radius, corner_to_point, radar_coords_to_cart, dms_to_d


def create_flat_array(sweep):
    ppi = N.zeros([sweep.h.nrays, sweep.rays[0].h.nbins],
                  dtype=float) + 1.31072000e+05
    for raynum in range(sweep.h.nrays):
        data = sweep.rays[raynum].data
        ppi[raynum, 0:len(data)] = sweep.rays[raynum].data
    return ppi


def plot_sur(xsapr, imagefilename, var, sweep, **kwargs):
    debug = kwargs.get('debug', False)
    rges = kwargs.get('rges', {
        'CZ': [-16., 64.],
        'DZ': [-16., 64.],
        'ZT': [-16.0, 64.0],
        'VE': [-40.0, 40.0],
        'VR': [-16, 16],
        'PH': [300, 350],
        'DR': [-4, 3]}[var])
    var_lab = kwargs.get('var_lab', {
        'CZ': 'Eq refl fact (dBz)',
        'DZ': 'Eq refl fact (dBz)',
        'ZT': 'Eq refl fact (dBz)',
        'PH': 'Phidp deg',
        'VE': 'Unfolded VR (m/s)',
        'DR': 'ZDR dB',
        'VR': 'Radial Velocity (+away) m/s'}[var])
    locs = kwargs.get('locs', [])
    labels = kwargs.get('labels', [])
    print rges
    field = _rsl.fieldTypes().list.index(var)
    plain_data = create_flat_array(xsapr.contents.volumes[field].sweeps[sweep])
    azmths = array(
        [xsapr.contents.volumes[field].sweeps[sweep].rays[i].h.azimuth
         for i in range(xsapr.contents.volumes[field].sweeps[sweep].h.nrays)])

    ranges = xsapr.contents.volumes[field].sweeps[sweep].rays[0].dists / 1000.0
    elevs = array(
        [xsapr.contents.volumes[field].sweeps[sweep].rays[i].h.elev
         for i in range(xsapr.contents.volumes[field].sweeps[sweep].h.nrays)])
    rg, azg = meshgrid(ranges, azmths)
    ele = (ones(rg.shape, dtype=float) *
           xsapr.contents.volumes[field].sweeps[sweep].h.elev)
    x, y, z = radar_coords_to_cart(rg, azg, ele)
    titl = kwargs.get('titl',
                      "PPI %(radar_name)s %(year)04d-%(month)02d-%(day)02d %(hour)02d:%(min)02d:%(sec)02d Elevation: %(el).2f")
    rad_head = xsapr.contents.h
    titl_comp = titl % {
        'radar_name': rad_head.radar_name,
        'year': rad_head.year,
        'month': rad_head.month,
        'day': rad_head.day,
        'hour': rad_head.hour,
        'min': rad_head.minute,
        'sec': rad_head.sec,
        'el': elevs[0]}
    f = figure(figsize=[7, 5.5])
    pp = pcolormesh(x / 1000.0, y / 1000.0,
                    N.ma.masked_outside(plain_data, rges[0], rges[1]),
                    vmin=rges[0], vmax=rges[1])
    ylabel('Distance North from radar (km)')
    xlabel('Distance South from radar (km)')
    gca().set_aspect('equal')
    title(titl_comp)
    radar_loc = [dms_to_d([rad_head.latd, rad_head.latm, rad_head.lats]),
                 dms_to_d([rad_head.lond, rad_head.lonm, rad_head.lons])]
    for i in range(len(locs)):
        carts = corner_to_point(radar_loc, loc[i])
        plot([carts[0] / 1000.0, carts[0] / 1000.0], [carts[1]/1000.0,
             carts[1]/1000.0], 'r+')
        text(carts[0] - 5.0, carts[1], label[i])
    cbax = f.add_axes([.9, .1, 0.02, .8])
    c2 = colorbar(cax=cbax, mappable=pp)
    ylabel(var_lab)
    savefig(imagefilename)
    close(f)


def plot_sur_masked(xsapr, imagefilename, var, masking_var, masking_value,
                    sweep, **kwargs):
    debug = kwargs.get('debug', False)
    rges = kwargs.get('rges', {
        'CZ': [-16., 64.],
        'ZT': [-16.0, 64.0],
        'VE': [-40.0, 40.0]}[var])
    var_lab = kwargs.get('var_lab', {
        'CZ': 'Eq refl fact (dBz)',
        'ZT': 'Eq refl fact (dBz)',
        'VE': 'Unfolded VR (m/s)'}[var])
    locs = kwargs.get('locs', [])
    labels = kwargs.get('labels', [])
    print rges
    field = _rsl.fieldTypes().list.index(var)
    mask_field = _rsl.fieldTypes().list.index(masking_var)
    plain_data = create_flat_array(xsapr.contents.volumes[field].sweeps[sweep])
    plain_mask = create_flat_array(
        xsapr.contents.volumes[mask_field].sweeps[sweep])
    azmths = array(
        [xsapr.contents.volumes[field].sweeps[sweep].rays[i].h.azimuth
         for i in range(xsapr.contents.volumes[field].sweeps[sweep].h.nrays)])
    ranges = xsapr.contents.volumes[field].sweeps[sweep].rays[0].dists / 1000.0
    elevs = array(
        [xsapr.contents.volumes[field].sweeps[sweep].rays[i].h.elev
         for i in range(xsapr.contents.volumes[field].sweeps[sweep].h.nrays)])
    rg, azg = meshgrid(ranges, azmths)
    ele = (ones(rg.shape, dtype=float) *
           xsapr.contents.volumes[field].sweeps[sweep].h.elev)
    x, y, z = radar_coords_to_cart(rg, azg, ele)
    titl = kwargs.get('titl',
                      "PPI %(radar_name)s %(year)04d-%(month)02d-%(day)02d %(hour)02d:%(min)02d:%(sec)02d Elevation: %(el).2f")
    rad_head = xsapr.contents.h
    titl_comp = titl % {
        'radar_name': rad_head.radar_name,
        'year': rad_head.year,
        'month': rad_head.month,
        'day': rad_head.day,
        'hour': rad_head.hour,
        'min': rad_head.minute,
        'sec': rad_head.sec,
        'el': elevs[0]}
    f = figure(figsize=[7, 5.5])
    pp = pcolormesh(x / 1000.0, y / 1000.0,
                    N.ma.masked_where(plain_mask < masking_value, plain_data),
                    vmin=rges[0], vmax=rges[1])
    ylabel('Distance North from radar (km)')
    xlabel('Distance East from radar (km)')
    gca().set_aspect('equal')
    title(titl_comp)
    radar_loc = [dms_to_d([rad_head.latd, rad_head.latm, rad_head.lats]),
                 dms_to_d([rad_head.lond, rad_head.lonm, rad_head.lons])]
    for i in range(len(locs)):
        carts = corner_to_point(radar_loc, locs[i])
        plot([carts[0] / 1000.0, carts[0] / 1000.0],
             [carts[1] / 1000.0, carts[1] / 1000.0], 'r+')
        text(carts[0] / 1000.0 - 5.0, carts[1] / 1000.0, labels[i])
    cbax = f.add_axes([.83, .1, 0.02, .8])
    c2 = colorbar(cax=cbax, mappable=pp)
    ylabel(var_lab)
    savefig(imagefilename)
    close(f)
