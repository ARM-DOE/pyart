""" Unit Tests for Py-ART's graph/radarmapdisplay.py module. """
# execute this script to create figure_plot_radar_display_*.png files.

# TODO use matplotlib's @image_comparison decorator to compare to file
# in baseline_images directory. Current this test only determines if files can
# be created, not that they are correct.

import matplotlib.pyplot as plt
import numpy as np
import pytest

import pyart


# Top level Figure generating tests
@pytest.mark.skipif(not pyart.graph.radarmapdisplay_basemap._BASEMAP_AVAILABLE,
                    reason="Basemap is not installed.")
def test_radarmapdisplay_ppi_basemap(outfile=None):
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarMapDisplayBasemap(radar, shift=(0.1, 0.0))
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
@pytest.mark.skipif(not pyart.graph.radarmapdisplay_basemap._BASEMAP_AVAILABLE,
                    reason="Basemap is not installed.")
def test_radarmapdisplay_auto_range_basemap():
    # test the auto_range=True function
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarMapDisplayBasemap(radar, shift=(0.1, 0.0))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi_map('reflectivity_horizontal', resolution='c')
    assert round(display.basemap.latmax, 2) == 36.85
    assert round(display.basemap.latmin, 2) == 36.13
    assert round(display.basemap.lonmax, 2) == -97.15
    assert round(display.basemap.lonmin, 2) == -98.04
    plt.close()


@pytest.mark.skipif(not pyart.graph.radarmapdisplay_basemap._BASEMAP_AVAILABLE,
                    reason="Basemap is not installed.")
def test_error_raising_basemap():
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarMapDisplayBasemap(radar, shift=(0.1, 0.0))
    # no basemap
    pytest.raises(ValueError, display.plot_range_ring, 10)


@pytest.mark.skipif(not pyart.graph.radarmapdisplay_basemap._BASEMAP_AVAILABLE,
                    reason="Basemap is not installed.")
def test_radardisplay_cylindrical_proj_error_basemap():
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarMapDisplayBasemap(radar)
    pytest.raises(
        ValueError,
        display.plot_ppi_map, 'reflectivity_horizontal', projection='cyl',
        resolution='c')


if __name__ == "__main__":
    test_radarmapdisplay_ppi_basemap('figure_radarmapdisplay_ppi.png')
