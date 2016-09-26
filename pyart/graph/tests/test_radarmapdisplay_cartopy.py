""" Unit Tests for Py-ART's graph/radarmapdisplay.py module. """
# execute this script to create figure_plot_radar_display_*.png files.

# TODO use matplotlib's @image_comparison decorator to compare to file
# in baseline_images directory. Current this test only determines if files can
# be created, not that they are correct.

import matplotlib.pyplot as plt
import pyart
import numpy as np
from numpy.testing import assert_raises
from numpy.testing.decorators import skipif


# Top level Figure generating tests
@skipif(not pyart.graph.radarmapdisplay_cartopy._CARTOPY_AVAILABLE)
def test_radarmapdisplay_cartopy_ppi(outfile=None):
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarMapDisplayCartopy(radar, shift=(0.1, 0.0))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi_map(
        'reflectivity_horizontal', 0, colorbar_flag=True,
        title="Fancy PPI Map", mask_tuple=('reflectivity_horizontal', -100),
        resolution='c', min_lon=-100, max_lon=-93, min_lat=33, max_lat=38,
        mask_outside=True)
    display.plot_point(-95, 35, label_text='TEXT')
    display.plot_range_rings([15, 30])
    display.plot_line_geo(np.array([-95, -95]), np.array([33, 38]))
    if outfile:
        fig.savefig(outfile)
    plt.close()


# Tests of methods, these tests do not generate figures
@skipif(not pyart.graph.radarmapdisplay_cartopy._CARTOPY_AVAILABLE)
def test_radarmapdisplay_cartopy_auto_range():
    # test the auto_range=True function
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarMapDisplayCartopy(radar, shift=(0.1, 0.0))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi_map('reflectivity_horizontal', resolution='c')
    assert round(display.basemap.latmax, 2) == 36.85
    assert round(display.basemap.latmin, 2) == 36.13
    assert round(display.basemap.lonmax, 2) == -97.15
    assert round(display.basemap.lonmin, 2) == -98.04
    plt.close()


@skipif(not pyart.graph.radarmapdisplay_cartopy._CARTOPY_AVAILABLE)
def test_error_raising():
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarMapDisplayCartopy(radar, shift=(0.1, 0.0))
    # no basemap
    assert_raises(ValueError, display.plot_range_ring, 10)


@skipif(not pyart.graph.radarmapdisplay_cartopy._CARTOPY_AVAILABLE)
def test_radardisplay_cylindrical_proj_error():
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarMapDisplayCartopy(radar)
    assert_raises(
        ValueError,
        display.plot_ppi_map, 'reflectivity_horizontal', projection='cyl',
        resolution='c')


if __name__ == "__main__":
    test_radarmapdisplay_ppi('figure_radarmapdisplay_ppi.png')
