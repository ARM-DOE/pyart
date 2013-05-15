#!/usr/bin/env python
"""
Example script for plotting RHIs from the CMAC VAP
"""

import argparse

import netCDF4
import matplotlib.pyplot as plt
import pyart

if __name__ == "__main__":

    # parse the command line arguments
    parser = argparse.ArgumentParser(
        description="Plot all RHIs in a CMAC VAP NetCDF file.  A figure "
                    "named prefix_%timestr.png is created.")
    parser.add_argument("filename", type=str,
                        help="CMAC VAP NetCDF file to plot.")
    parser.add_argument("prefix", type=str,
                        help="Prefix in filename, example 'figure_'.")
    args = parser.parse_args()

    # read the data and create the display object
    radar = pyart.io.read_netcdf(args.filename)
    radar.metadata['instrument_name'] = 'XSAPR'  # XXX Hack
    display = pyart.graph.RadarDisplay(radar)

    # create the figure
    fig = plt.figure(figsize=[12, 17])
    fig.subplots_adjust(hspace=0.4)
    xlabel = 'Distance from radar (km)'
    ylabel = 'Height agl (km)'
    colorbar_label = 'Hz. Eq. Refl. Fac. (dBZ)'

    # plot each RHI
    for snum in radar.sweep_number['data']:

        fixed_angle = radar.fixed_angle['data'][snum]
        title = 'HSRHI Az=%.3f' % (fixed_angle)
        ax = fig.add_subplot(radar.nsweeps, 1, snum+1)
        display.plot_rhi('reflectivity_horizontal', snum, vmin=-20, vmax=20,
                         mask_outside=False, title=title,
                         axislabels=(xlabel, ylabel),
                         colorbar_label=colorbar_label, ax=ax)
        display.set_limits(ylim=[0, 15], ax=ax)

    # add a figure title
    figure_title = 'Time: ' + display.time_begin.isoformat() + 'Z'
    fig.text(0.35, 0.92, figure_title)

    # save the figure
    time_text = display.time_begin.strftime('%Y%m%d.%H%M%S')
    fname = args.prefix + time_text + '.png'
    fig.savefig(fname)
