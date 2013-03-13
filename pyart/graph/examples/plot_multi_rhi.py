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
    rslradar = _rsl.RSL_anyformat_to_radar(args.filename)
    display = pyart.graph.RslDisplay(rslradar)

    # fields to plot and ranges
    fields_to_plot = ['ZT', 'VR']
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
    radar_name = rslradar.contents.h.radar_name
    time_text = ' ' + display.time_begin.isoformat() + 'Z '
    field_num = _rsl.fieldTypes().list.index('ZT')
    azimuth = rslradar.contents.volumes[field_num].sweeps[0].rays[0].azimuth
    title = 'RHI ' + radar_name + time_text + 'Azimuth %.2f' % (azimuth)
    plt.suptitle(title)

    plt.savefig(args.figurename)
