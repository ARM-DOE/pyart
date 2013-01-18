#!/usr/bin/env python
"""
Example command line utility for plotting two panels from a radar file

usage : dynamics_plot.py filename outdir tilt

"""

import copy
import sys
import os

import matplotlib
from pylab import *
import netCDF4

from pyart.graph import radar_display
from pyart.io import radar, py4dd, py_mdv

if __name__ == "__main__":

    # read in command line arguments
    filename = sys.argv[1]
    outdir = sys.argv[2]
    tilt = int(sys.argv[3])

    print "plotting ", filename, " to ", outdir

    # read in the radar file
    if ".mdv" in filename:
        my_object = py_mdv.read_mdv(filename, debug=True)
        myradar = radar.Radar(my_object)

    elif ".nc" in filename:
        my_object = netCDF4.Dataset(filename)
        myradar = radar.Radar(my_object)
        my_object.close()

    else:
        py4dd.RSL_radar_verbose_on()
        my_object = py4dd.RSL_anyformat_to_radar(filename)
        myradar = radar.Radar(my_object)

    # create the plots
    myradar.metadata.update({'instrument_name': "Radar"})
    my_display = radar_display.radar_display(myradar)
    f = figure(figsize=[12, 5])
    subplot(1, 2, 1)
    my_display.plot_ppi('mean_doppler_velocity', tilt)
    my_display.add_cb()
    subplot(1, 2, 2)
    my_display.plot_ppi('corrected_mean_doppler_velocity', tilt)
    my_display.add_cb()
    figname = my_display.generate_filename('six_panel', tilt)
    savefig(outdir + '/' + figname.replace(' ', '_'))
    close(f)
