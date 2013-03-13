#! /usr/bin/env python
# Create a PPI plot from a MDV file.
# Example usage: ./vanilla_ppi.py 110635.mdv test_mdv_ppi.png

import sys
import getopt

import numpy as np
import matplotlib.pyplot as plt

import pyart.io.mdv as mdv
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
    
    outfile = 'test_kml.png'

    # read in the data
    mdvfile = mdv.MdvFile(filename, debug=True)
    mdvdisplay = plot_mdv.MdvDisplay(mdvfile)   

    # create the figure
    fig = plt.figure(figsize=[5, 5])
    ax = fig.add_subplot(111)
    mdvdisplay.plot_ppi('DBZ_F', 0, mask=['NCP_F', 0.5])
    mdvdisplay.set_limits(ylim=[-120, 120], xlim=[-120, 120])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    fig.savefig(outfile, transparent=True)
    mdvfile.close()
