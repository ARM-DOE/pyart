#! /usr/bin/env python
"""
Example script for plotting RHI of multiple fields from a Sigmet file from the
SGP XSAPR.
"""

import argparse

import matplotlib.pyplot as plt

import pyart
from pyart.io import _rsl

if __name__ == "__main__":

    # parse the command line arguments
    parser = argparse.ArgumentParser(
        description="Plot the first RHI from a Sigmet file.")
    parser.add_argument("filename", type=str, help="Sigmet file to plot.")
    parser.add_argument("figurename", type=str,
                        help="Filename of figure to create.")
    parser.add_argument('-d', '--debug', help="Show debugging information",
                        action='store_true')
    args = parser.parse_args()

    if args.debug:
        _rsl.RSL_radar_verbose_on()

    # read the data and create the display object
    radar = pyart.io.read_rsl(args.filename)
    display = pyart.graph.RadarDisplay(radar)

    # fields to plot and ranges
    fields_to_plot = ['reflectivity_horizontal', 'mean_doppler_velocity']
    ranges = [(-32, 64), (-17.0, 17.0)]

    # plot the data
    nplots = len(fields_to_plot)
    plt.figure(figsize=[5*nplots, 4])

    # plot each field
    for plot_num in xrange(nplots):
        field = fields_to_plot[plot_num]
        vmin, vmax = ranges[plot_num]

        plt.subplot(1, nplots, plot_num + 1)
        display.plot_rhi(field, 0, vmin=vmin, vmax=vmax, title_flag=False)
        display.set_limits(ylim=[0, 17])

    # set the figure title
    radar_name = display.radar_name
    time_text = ' ' + display.time_begin.isoformat() + 'Z '
    azimuth = radar.sweep_info['fixed_angle']['data'][0]
    title = 'RHI ' + radar_name + time_text + 'Azimuth %.2f' % (azimuth)
    plt.suptitle(title)

    plt.savefig(args.figurename)
