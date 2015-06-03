""" Unit Tests for Py-ART's graph/plot_rsl.py module. """
# execute this script to create figure_plot_rsl_*.png files.
# To compare results:
# display figure_plot_rsl_ppi.png & display figure_plot_rsl_ppi_radar.png &
# display figure_plot_rsl_rhi.png & display figure_plot_rsl_rhi_radar.png &

# TODO use matplotlib's @image_comparison decorator to compare to file
# in baseline_images directory. Current this test only determines if files can
# be created, not that they are correct.

import matplotlib.pyplot as plt
from numpy.testing.decorators import skipif
import pyart


@skipif(not pyart.io._RSL_AVAILABLE)
def test_plot_rsl_rhi(outfile=None):
    rslfile = pyart.io._rsl_interface.RslFile(
        pyart.testing.SIGMET_RHI_FILE.encode('ascii'))
    display = pyart.graph.RslDisplay(rslfile)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('DZ', 0, ax=ax)
    if outfile:
        fig.savefig(outfile)


@skipif(not pyart.io._RSL_AVAILABLE)
def test_plot_rsl_ppi(outfile=None):
    rslfile = pyart.io._rsl_interface.RslFile(
        pyart.testing.SIGMET_PPI_FILE.encode('ascii'))
    display = pyart.graph.RslDisplay(rslfile)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('DZ', 0, ax=ax)
    if outfile:
        fig.savefig(outfile)


# version of above which use the RadarDisplay to return similar images.


@skipif(not pyart.io._RSL_AVAILABLE)
def test_plot_rsl_radar_rhi(outfile=None):
    radar = pyart.io.read_rsl(pyart.testing.SIGMET_RHI_FILE)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('reflectivity', 0,
                     vmin=-16, vmax=64, ax=ax)
    if outfile:
        fig.savefig(outfile)


@skipif(not pyart.io._RSL_AVAILABLE)
def test_plot_rsl_radar_ppi(outfile=None):
    radar = pyart.io.read_rsl(pyart.testing.SIGMET_PPI_FILE)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('reflectivity', 0,
                     vmin=-16, vmax=64, ax=ax)
    if outfile:
        fig.savefig(outfile)


if __name__ == "__main__":
    test_plot_rsl_rhi('figure_plot_rsl_rhi.png')
    test_plot_rsl_ppi('figure_plot_rsl_ppi.png')
    test_plot_rsl_radar_rhi('figure_plot_rsl_rhi_radar.png')
    test_plot_rsl_radar_ppi('figure_plot_rsl_ppi_radar.png')
