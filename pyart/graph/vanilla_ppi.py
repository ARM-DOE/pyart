#! /usr/bin/env python

import sys
import getops

import numpy as np
from pylab import *

import pyart.io.py_mdv as py_mdv
import pyart.io.py4dd as py4dd
from pyart.graph import plot_mdv as plot_mdv


def get_optargs(argv):
    try:
        opts, args = getopt.getopt(argv, "d", ['debug'])
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

    t = float(opt_dict.get('-z', opt_dict.get('--zlim1', '17')))
    var = opt_dict.get('-v', opt_dict.get('--var', 'DBZ_F'))
    if debug:
        print "Vanilla ppi "
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
    myfile = py_mdv.read_mdv(filename, debug=True)
    f = figure(figsize=[5, 5])
    plot_mdv.single_panel_ppi(myfile, 0, 'DBZ_F', mask=['NCP_F', 0.5],
                              ylim=[-120, 120], xlim=[-120, 120])
    gca().axes.get_xaxis().set_visible(False)
    gca().axes.get_yaxis().set_visible(False)
    f.get_axes()[1].set_visible(False)
    subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    savefig('/home/sc8/scratch/test_kml.png', transparent=True)
    close(f)
    myfile.close()
