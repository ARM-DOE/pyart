""" Unit Tests for Py-ART's graph/plot_cfradial.py module. """
# execute this script to create figure_plot_cfradial_*.png files.
# To compare results
# display figure_plot_netcdf_ppi.png&
# display figure_plot_cfradial_ppi_radar.png&
# display figure_plot_netcdf_rhi.png&
# display figure_plot_cfradial_rhi_radar.png&

# TODO use matplotlib's @image_comparison decorator to compare to file
# in baseline_images directory. Current this test only determines if files can
# be created, not that they are correct.

import netCDF4
import matplotlib.pyplot as plt
import pyart


def test_plot_netcdf_rhi(outfile=None):
    dataset = netCDF4.Dataset(pyart.testing.CFRADIAL_RHI_FILE)
    display = pyart.graph.CFRadialDisplay(dataset)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('reflectivity_horizontal', 0, ax=ax)
    if outfile:
        fig.savefig(outfile)


def test_plot_netcdf_ppi(outfile=None):
    dataset = netCDF4.Dataset(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.CFRadialDisplay(dataset)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('reflectivity_horizontal', 0, ax=ax)
    if outfile:
        fig.savefig(outfile)


def test_plot_netcdf_ppi_masked(outfile=None):
    dataset = netCDF4.Dataset(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.CFRadialDisplay(dataset)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('reflectivity_horizontal', 0, ax=ax,
                     mask_tuple=('reflectivity_horizontal', 16.0))
    if outfile:
        fig.savefig(outfile)


# version of above which use the RadarDisplay to return similar images.


def test_plot_netcdf_radar_rhi(outfile=None):
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_RHI_FILE)
    radar.metadata['instrument_name'] = ''
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('reflectivity_horizontal', 0, ax=ax)
    if outfile:
        fig.savefig(outfile)


def test_plot_netcdf_radar_ppi(outfile=None):
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('reflectivity_horizontal', 0, ax=ax)
    if outfile:
        fig.savefig(outfile)


def test_plot_netcdf_radar_ppi_masked(outfile=None):
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('reflectivity_horizontal', 0, ax=ax,
                     mask_tuple=('reflectivity_horizontal', 16.0))
    if outfile:
        fig.savefig(outfile)


if __name__ == "__main__":
    test_plot_netcdf_rhi('figure_plot_cfradial_rhi.png')
    test_plot_netcdf_ppi('figure_plot_cfradial_ppi.png')
    test_plot_netcdf_ppi_masked('figure_plot_cfradial_ppi_masked.png')
    test_plot_netcdf_radar_rhi('figure_plot_cfradial_rhi_radar.png')
    test_plot_netcdf_radar_ppi('figure_plot_cfradial_ppi_radar.png')
    test_plot_netcdf_radar_ppi_masked(
        'figure_plot_cfradial_ppi_masked_radar.png')
