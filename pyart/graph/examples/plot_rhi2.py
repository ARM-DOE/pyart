#! /usr/bin/env python
"""
Example script for plotting a RHI from a Sigmet file from the SGP X-SAPR.
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

    # create the figure
    fig = plt.figure(figsize=[10, 4])
    ax = fig.add_subplot(111)

    # plot the data
    radar_name = radar.metadata['instrument_name']
    time_text = ' ' + display.time_begin.isoformat() + 'Z '
    azimuth = radar.fixed_angle['data'][0]
    title = 'RHI ' + radar_name + time_text + 'Azimuth %.2f' % (azimuth)
    colorbar_label = 'Eq refl fact (dBz)'

    display.plot_rhi('reflectivity_horizontal', 0, vmin=-32, vmax=64,
                     title=title, colorbar_flag=False, ax=ax)
    display.set_limits(ylim=[0, 17])

    cax = fig.add_axes([.9, .1, 0.02, .8])
    display.plot_colorbar(fig=fig, cax=cax, label=colorbar_label)

    # save the figure
    fig.savefig(args.figurename)
