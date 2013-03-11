#! /usr/bin/env python
# nose will check that figures can be created
# execute this script to create figure_plot_mdv_*.png files.
# To compare results
# display figure_plot_mdv_ppi.png & display figure_plot_mdv_ppi_radar.png &
# display figure_plot_mdv_rhi.png & display figure_plot_mdv_rhi_radar.png &

import matplotlib.pyplot as plt
import pyart

def test_plot_mdv_rhi(outfile=None):
    mdvfile = pyart.io.mdv.MdvFile('110041.mdv')
    display = pyart.graph.MdvDisplay(mdvfile)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('DBZ_F', 0)
    if outfile:
        fig.savefig(outfile)


def test_plot_mdv_ppi(outfile=None):
    mdvfile = pyart.io.mdv.MdvFile('110635.mdv')
    display = pyart.graph.MdvDisplay(mdvfile)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('DBZ_F', 0)
    if outfile:
        fig.savefig(outfile)


# version of above which use the RadarDisplay to return similar images.


def test_plot_mdv_radar_rhi(outfile=None):
    radar = pyart.io.read_mdv('110041.mdv')
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('reflectivity_horizontal', 0, vmin=-16, vmax=64)
    if outfile:
        fig.savefig(outfile)


def test_plot_mdv_radar_ppi(outfile=None):
    radar = pyart.io.read_mdv('110635.mdv')
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('reflectivity_horizontal', 0, vmin=-16, vmax=64)
    if outfile:
        fig.savefig(outfile)


if __name__ == "__main__":
    test_plot_mdv_rhi('figure_plot_mdv_rhi.png')
    test_plot_mdv_ppi('figure_plot_mdv_ppi.png')
    test_plot_mdv_radar_rhi('figure_plot_mdv_rhi_radar.png')
    test_plot_mdv_radar_ppi('figure_plot_mdv_ppi_radar.png')
