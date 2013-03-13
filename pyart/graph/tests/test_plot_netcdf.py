#! /usr/bin/env python
# nose will check that figures can be created
# execute this script to create figure_plot_netcdf_*.png files.
# To compare results
# display figure_plot_netcdf_ppi.png& display figure_plot_netcdf_ppi_radar.png&
# display figure_plot_netcdf_rhi.png& display figure_plot_netcdf_rhi_radar.png&

import netCDF4
import matplotlib.pyplot as plt
import pyart


def test_plot_netcdf_rhi(outfile=None):
    rhi_file = 'sgpxsaprrhicmacI5.c0.20110524.015604_NC4.nc'
    dataset = netCDF4.Dataset(rhi_file)
    display = pyart.graph.NetcdfDisplay(dataset)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('reflectivity_horizontal', 0)
    if outfile:
        fig.savefig(outfile)


def test_plot_netcdf_ppi(outfile=None):
    ppi_file = 'sgpxsaprsesurcmacI4.c0.20110520.105511.nc'
    dataset = netCDF4.Dataset(ppi_file)
    display = pyart.graph.NetcdfDisplay(dataset)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('reflectivity_horizontal', 0)
    if outfile:
        fig.savefig(outfile)


# version of above which use the RadarDisplay to return similar images.


def test_plot_netcdf_radar_rhi(outfile=None):
    rhi_file = 'sgpxsaprrhicmacI5.c0.20110524.015604_NC4.nc'
    radar = pyart.io.read_netcdf(rhi_file)
    radar.metadata['instrument_name'] = ''
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('reflectivity_horizontal', 0)
    if outfile:
        fig.savefig(outfile)


def test_plot_netcdf_radar_ppi(outfile=None):
    ppi_file = 'sgpxsaprsesurcmacI4.c0.20110520.105511.nc'
    radar = pyart.io.read_netcdf(ppi_file)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('reflectivity_horizontal', 0)
    if outfile:
        fig.savefig(outfile)


if __name__ == "__main__":
    test_plot_netcdf_rhi('figure_plot_netcdf_rhi.png')
    test_plot_netcdf_ppi('figure_plot_netcdf_ppi.png')
    test_plot_netcdf_radar_rhi('figure_plot_netcdf_rhi_radar.png')
    test_plot_netcdf_radar_ppi('figure_plot_netcdf_ppi_radar.png')
