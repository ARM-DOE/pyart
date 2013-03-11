#! /usr/bin/env python
# nose will check that figures can be created
# execute this script to create figure_plot_mdv_*.png files.
# To compare results
# display figure_plot_rsl_ppi.png & display figure_plot_rsl_ppi_radar.png &
# display figure_plot_rsl_rhi.png & display figure_plot_rsl_rhi_radar.png &

import matplotlib.pyplot as plt
import pyart

def test_plot_rsl_rhi(outfile=None):
    rslradar = pyart.io._rsl.RSL_anyformat_to_radar('XSW110520113537.RAW7HHL')
    display = pyart.graph.RslDisplay(rslradar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('DZ', 0)
    if outfile:
        fig.savefig(outfile)


def test_plot_rsl_ppi(outfile=None):
    rslradar = pyart.io._rsl.RSL_anyformat_to_radar('XSW110520105408.RAW7HHF')
    display = pyart.graph.RslDisplay(rslradar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('DZ', 0)
    if outfile:
        fig.savefig(outfile)


# version of above which use the RadarDisplay to return similar images.


def test_plot_rsl_radar_rhi(outfile=None):
    radar = pyart.io.read_rsl('XSW110520113537.RAW7HHL')
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('reflectivity_horizontal', 0, vmin=-16, vmax=64)
    if outfile:
        fig.savefig(outfile)


def test_plot_rsl_radar_ppi(outfile=None):
    radar = pyart.io.read_rsl('XSW110520105408.RAW7HHF')
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('reflectivity_horizontal', 0, vmin=-16, vmax=64)
    if outfile:
        fig.savefig(outfile)


if __name__ == "__main__":
    test_plot_rsl_rhi('figure_plot_rsl_rhi.png')
    test_plot_rsl_ppi('figure_plot_rsl_ppi.png')
    test_plot_rsl_radar_rhi('figure_plot_rsl_rhi_radar.png')
    test_plot_rsl_radar_ppi('figure_plot_rsl_ppi_radar.png')
