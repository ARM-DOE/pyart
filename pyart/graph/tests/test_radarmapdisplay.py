""" Unit Tests for Py-ART's graph/radarmapdisplay.py module. """
# execute this script to create figure_plot_radar_display_*.png files.

# TODO use matplotlib's @image_comparison decorator to compare to file
# in baseline_images directory. Current this test only determines if files can
# be created, not that they are correct.

import matplotlib.pyplot as plt
import pyart
from numpy.testing import assert_raises


# Top level Figure generating tests
def test_radarmapdisplay_ppi(outfile=None):
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarMapDisplay(radar, shift=(0.1, 0.0))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi_map(
        'reflectivity_horizontal', 0, colorbar_flag=True,
        title="Fancy PPI Map", mask_tuple=('reflectivity_horizontal', -100),
        resolution='c', auto_range=False,
        min_lon=-100, max_lon=-93, min_lat=33, max_lat=38)
    if outfile:
        fig.savefig(outfile)
    plt.close()


# Tests of methods, these tests do not generate figures
def test_radarmapdisplay_auto_range():
    # test the auto_range=True function
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarMapDisplay(radar, shift=(0.1, 0.0))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi_map('reflectivity_horizontal', resolution='c')
    assert round(display.basemap.latmax, 2) == 36.85
    assert round(display.basemap.latmin, 2) == 36.13
    assert round(display.basemap.lonmax, 2) == -97.15
    assert round(display.basemap.lonmin, 2) == -98.04
    plt.close()


if __name__ == "__main__":
    test_radarmapdisplay_ppi('figure_radarmapdisplay_ppi.png')
