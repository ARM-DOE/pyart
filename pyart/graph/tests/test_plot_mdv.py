""" Unit Tests for Py-ART's graph/plot_mdv.py module. """
# execute this script to create figure_plot_mdv_*.png files.
# To compare results:
# display figure_plot_mdv_ppi.png & display figure_plot_mdv_ppi_radar.png &
# display figure_plot_mdv_rhi.png & display figure_plot_mdv_rhi_radar.png &
# display figure_plot_mdv_rhi_masked.png &
# display figure_plot_mdv_rhi_masked_radar.png &

# TODO use matplotlib's @image_comparison decorator to compare to file
# in baseline_images directory. Current this test only determines if files can
# be created, not that they are correct.

import matplotlib.pyplot as plt
import pyart


def test_plot_mdv_rhi(outfile=None):
    mdvfile = pyart.io.mdv_common.MdvFile(pyart.testing.MDV_RHI_FILE)
    display = pyart.graph.MdvDisplay(mdvfile)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('DBZ_F', 0, ax=ax)
    if outfile:
        fig.savefig(outfile)


def test_plot_mdv_rhi_masked(outfile=None):
    mdvfile = pyart.io.mdv_common.MdvFile(pyart.testing.MDV_RHI_FILE)
    display = pyart.graph.MdvDisplay(mdvfile)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('DBZ_F', 0, mask_tuple=('DBZ_F', 16.), ax=ax)
    if outfile:
        fig.savefig(outfile)


def test_plot_mdv_ppi(outfile=None):
    mdvfile = pyart.io.mdv_common.MdvFile(pyart.testing.MDV_PPI_FILE)
    display = pyart.graph.MdvDisplay(mdvfile)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('DBZ_F', 0, ax=ax)
    if outfile:
        fig.savefig(outfile)


# version of above which use the RadarDisplay to return similar images.


def test_plot_mdv_radar_rhi(outfile=None):
    radar = pyart.io.read_mdv(pyart.testing.MDV_RHI_FILE)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('reflectivity', 0, vmin=-16, vmax=64, ax=ax)
    if outfile:
        fig.savefig(outfile)


def test_plot_mdv_radar_rhi_masked(outfile=None):
    radar = pyart.io.read_mdv(pyart.testing.MDV_RHI_FILE)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('reflectivity', 0, vmin=-16, vmax=64,
                     mask_tuple=('reflectivity', 16.), ax=ax)
    if outfile:
        fig.savefig(outfile)


def test_plot_mdv_radar_ppi(outfile=None):
    radar = pyart.io.read_mdv(pyart.testing.MDV_PPI_FILE)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('reflectivity', 0, vmin=-16, vmax=64, ax=ax)
    if outfile:
        fig.savefig(outfile)


if __name__ == "__main__":
    test_plot_mdv_rhi('figure_plot_mdv_rhi.png')
    test_plot_mdv_rhi_masked('figure_plot_mdv_rhi_masked.png')
    test_plot_mdv_ppi('figure_plot_mdv_ppi.png')
    test_plot_mdv_radar_rhi('figure_plot_mdv_rhi_radar.png')
    test_plot_mdv_radar_rhi_masked('figure_plot_mdv_rhi_masked_radar.png')
    test_plot_mdv_radar_ppi('figure_plot_mdv_ppi_radar.png')
