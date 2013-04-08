#! /usr/bin/env python
# sample command line utility for plotting two panels from a radar file
# usage : two_panel_plot.py filename

import sys

import pylab as plt

from pyart.graph.radar_display import RadarDisplay 
from pyart.io import  radar, rsl, mdv

if __name__ == "__main__":

    # parse command line arguments
    filename = sys.argv[1]

    print "plotting ", filename, " to test.png"
    if ".mdv" in filename:
        my_object = mdv.read_mdv(filename, debug=True)
    else:
        radar_object = rsl.read_rsl(filename)

    my_display = RadarDisplay(radar_object)
    fig = plt.figure()
    plt.subplot(1, 2, 1)
    my_display.plot_ppi('reflectivity_horizontal', 1)

    plt.subplot(1, 2, 2)

    my_display.plot_ppi('mean_doppler_velocity', 1)
    plt.savefig('test.png')
    plt.close(fig)
