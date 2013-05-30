#! /usr/bin/env python
# nose will check that figures can be created
# execute this script to create figure_plot_mdv_*.png files.
# To compare results
# display figure_plot_rsl_ppi.png & display figure_plot_rsl_ppi_radar.png &
# display figure_plot_rsl_rhi.png & display figure_plot_rsl_rhi_radar.png &

import os.path

import matplotlib.pyplot as plt
import pyart


DIR = os.path.dirname(__file__)
RSL_RHI = os.path.join(DIR, 'XSW110520113537.RAW7HHL')
RSL_PPI = os.path.join(DIR, 'XSW110520105408.RAW7HHF')


def test_plot_rsl_rhi(outfile=None):
    rslfile = pyart.io._rsl_interface.RslFile(RSL_RHI)
    display = pyart.graph.RslDisplay(rslfile)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('DZ', 0)
    if outfile:
        fig.savefig(outfile)


def test_plot_rsl_ppi(outfile=None):
    rslfile = pyart.io._rsl_interface.RslFile(RSL_PPI)
    display = pyart.graph.RslDisplay(rslfile)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('DZ', 0)
    if outfile:
        fig.savefig(outfile)


# version of above which use the RadarDisplay to return similar images.


def test_plot_rsl_radar_rhi(outfile=None):
    radar = pyart.io.read_rsl(RSL_RHI)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('reflectivity_horizontal', 0, vmin=-16, vmax=64)
    if outfile:
        fig.savefig(outfile)


def test_plot_rsl_radar_ppi(outfile=None):
    radar = pyart.io.read_rsl(RSL_PPI)
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
