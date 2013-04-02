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
    rslradar = _rsl.RSL_anyformat_to_radar(args.filename)
    display = pyart.graph.RslDisplay(rslradar)

    # create the figure
    fig = plt.figure(figsize=[10, 4])
    ax = fig.add_subplot(111)

    # plot the data
    radar_name = rslradar.contents.h.radar_name
    time_text = ' ' + display.time_begin.isoformat() + 'Z '
    field_num = _rsl.fieldTypes().list.index('ZT')
    azimuth = rslradar.contents.volumes[field_num].sweeps[0].rays[0].azimuth
    title = 'RHI ' + radar_name + time_text + 'Azimuth %.2f' % (azimuth)
    colorbar_label = 'Eq refl fact (dBz)'

    display.plot_rhi('ZT', 0, vmin=-32, vmax=64, title=title,
                     colorbar_flag=False, ax=ax)
    display.set_limits(ylim=[0, 17])

    cax = fig.add_axes([.9, .1, 0.02, .8])
    display.plot_colorbar(fig=fig, cax=cax, label=colorbar_label)

    # save the figure and free the RSL radar
    fig.savefig(args.figurename)
    _rsl.RSL_free_radar(rslradar)
