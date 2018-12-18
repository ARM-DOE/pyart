""" Unit Tests for Py-ART's graph/radarmapdisplay.py module. """
# execute this script to create figure_plot_radar_display_*.png files.

# TODO use matplotlib's @image_comparison decorator to compare to file
# in baseline_images directory. Currently this test only determines if files
# can be created, not that they are correct.

import matplotlib.pyplot as plt
import numpy as np
import pytest

import pyart


# Top level Figure generating tests
@pytest.mark.skipif(not pyart.graph.radarmapdisplay._CARTOPY_AVAILABLE,
                    reason="Cartopy is not installed.")
def test_radarmapdisplay_cartopy_ppi(outfile=None):
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarMapDisplay(radar, shift=(0.1, 0.0))
    display.plot_ppi_map(
        'reflectivity_horizontal', 0, colorbar_flag=True,
        title="Fancy PPI Map", mask_tuple=('reflectivity_horizontal', -100),
        min_lon=-100, max_lon=-93, min_lat=33, max_lat=38,
        mask_outside=True)
    display.plot_point(-95, 35, label_text='TEXT')
    display.plot_range_rings([15, 30])
    display.plot_line_geo(np.array([-95, -95]), np.array([33, 38]))
    if outfile:
        plt.savefig(outfile)
    plt.close()


@pytest.mark.skipif(not pyart.graph.radarmapdisplay._CARTOPY_AVAILABLE,
                    reason="Cartopy is not installed.")
def test_radarmapdisplay_cartopy_preexisting_ax(outfile=None):
    import cartopy
    from cartopy.io.img_tiles import StamenTerrain
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarMapDisplay(radar, shift=(0.1, 0.0))
    fig = plt.figure()
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_image(StamenTerrain(), 6)
    display.plot_ppi_map('reflectivity_horizontal', 0, ax=ax, embelish=False)
    if outfile:
        fig.savefig(outfile)
    plt.close()


# Tests of methods, these tests do not generate figures
""" XXX Determine why this is failing
@pytest.mark.skipif(not pyart.graph.radarmapdisplay_cartopy._CARTOPY_AVAILABLE,
                    reason="Cartopy is not installed.")
def test_radarmapdisplay_cartopy_auto_range():
    # test the auto_range=True function
    import cartopy
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarMapDisplayCartopy(radar, shift=(0.1, 0.0))
    display.plot_ppi_map('reflectivity_horizontal')
    extent = display.ax.get_extent(cartopy.crs.PlateCarree())
    assert round(extent[3], 2) == 36.85
    assert round(extent[2], 2) == 36.13
    assert round(extent[1], 2) == -97.15
    assert round(extent[0], 2) == -98.04
    plt.close()
"""


@pytest.mark.skipif(not pyart.graph.radarmapdisplay._CARTOPY_AVAILABLE,
                    reason="Cartopy is not installed.")
def test_error_raising():
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarMapDisplay(radar, shift=(0.1, 0.0))
    # no cartopy
    pytest.raises(ValueError, display.plot_range_ring, 10)

if __name__ == "__main__":
    test_radarmapdisplay_cartopy_ppi(
        'figure_radarmapdisplay_cartopy_ppi.png')
