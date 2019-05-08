"""
General meteorological calculations useful to other modules.
"""

import os
import heapq

import numpy as np
import netCDF4

from pylab import datestr2num, date2num


def nth_smallest(n, iter):
    return heapq.nsmallest(n, iter)[-1]


def get_best_sounding(target, sdir, minl, maxl):
    sondes = os.listdir(sdir)
    sondes.sort()
    offsets = [np.abs(datestr2num(s[18:33].replace('.', ' ')) -
                      date2num(target)) for s in sondes]
    cont = True
    n = 1
    while cont:
        test_sonde = sondes[offsets.index(nth_smallest(n, offsets))]
        ncf_obj = netCDF4.Dataset(sdir + test_sonde, 'r')
        ncf_min = ncf_obj.variables['alt'][:].min()
        ncf_max = ncf_obj.variables['alt'][:].max()
        if ncf_min < minl and ncf_max > maxl:
            cont = False
            chosen_sonde = test_sonde
        ncf_obj.close()
        n = n + 1
    return chosen_sonde
