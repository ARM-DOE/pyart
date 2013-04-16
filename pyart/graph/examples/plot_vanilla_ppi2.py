#! /usr/bin/env python
"""
Example script for creating a PPI plot from a MDV file
"""

import argparse

import matplotlib.pyplot as plt

import pyart

if __name__ == "__main__":

    # parse the command line arguments
    parser = argparse.ArgumentParser(
        description="Create a PPI plot from a MDV file.")
    parser.add_argument("filename", type=str, help="MDV file to plot.")
    parser.add_argument("figurename", type=str,
                        help="Filename of figure to create.")
    parser.add_argument('-d', '--debug', help="Show debugging information",
                        action='store_true')
    args = parser.parse_args()

    # read in the data
    radar = pyart.io.read_mdv(args.filename)
    display = pyart.graph.RadarDisplay(radar)

    # create the figure
    fig = plt.figure(figsize=[5, 5])
    ax = fig.add_subplot(111, frameon=False)
    display.plot_ppi('reflectivity_horizontal', 0, vmin=-16.0, vmax=64.0,
                     mask_tuple=['norm_coherent_power', 0.5],
                     colorbar_flag=False, title_flag=False,
                     axislabels_flag=False)
    display.set_limits(ylim=[-120, 120], xlim=[-120, 120])
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    fig.savefig(args.figurename, transparent=True)
