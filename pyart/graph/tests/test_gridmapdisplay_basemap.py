""" Unit Tests for Py-ART's graph/gridmapdisplay_basemap.py module. """
# execute this script to create figure_gridmapdisplay_*.png files.

# TODO use matplotlib's @image_comparison decorator to compare to file
# in baseline_images directory. Current this test only determines if files can
# be created, not that they are correct.

import matplotlib.pyplot as plt
import pytest

import pyart

RESOLUTION = 'c'    # crude resolution to speed up tests


@pytest.mark.skipif(not pyart.graph.gridmapdisplay_basemap._BASEMAP_AVAILABLE,
                    reason="Basemap is not installed.")
def test_gridmapdisplay_simple_basemap(outfile=None):
    # test basic GridMapDisplay functionality
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplayBasemap(grid)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_basemap(ax=ax, resolution=RESOLUTION)
    display.plot_grid('reflectivity', vmin=-5., vmax=35)
    if outfile:
        fig.savefig(outfile)


@pytest.mark.skipif(not pyart.graph.gridmapdisplay_basemap._BASEMAP_AVAILABLE,
                    reason="Basemap is not installed.")
def test_gridmapdisplay_fancy_basemap(outfile=None):
    # test a bunch of GridMapDisplay functionality
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplayBasemap(grid)
    display.debug = True
    fig = plt.figure()

    ax = fig.add_subplot(331)
    display.plot_basemap(ax=ax, resolution=RESOLUTION)
    display.plot_grid(
        'reflectivity', vmin=-5., vmax=35, mask_outside=True,
        axislabels_flag=True, axislabels=('foo', 'bar'), title='Special title')
    display.plot_crosshairs(line_style='b--')

    ax = fig.add_subplot(332)
    display.plot_grid('reflectivity', axislabels_flag=True)

    ax = fig.add_subplot(333)
    display.plot_colorbar()

    ax = fig.add_subplot(334)
    display.plot_latitude_slice('reflectivity', mask_outside=True)

    ax = fig.add_subplot(335)
    display.plot_latitude_slice('reflectivity', title='Lat title')

    ax = fig.add_subplot(336)
    grid.fields['reflectivity']['valid_min'] = 0.
    grid.fields['reflectivity']['valid_max'] = 30.
    display.plot_longitude_slice('reflectivity', mask_outside=True)

    ax = fig.add_subplot(337)
    display.plot_longitude_slice('reflectivity', title='Lon title')

    ax = fig.add_subplot(338)
    del display.grid.fields['reflectivity']['long_name']
    display.plot_colorbar()

    if outfile:
        fig.savefig(outfile)


@pytest.mark.skipif(not pyart.graph.gridmapdisplay_basemap._BASEMAP_AVAILABLE,
                    reason="Basemap is not installed.")
def test_plot_basemap_not_using_corners(outfile=None):
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplayBasemap(grid)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_basemap(
        ax=ax, resolution=RESOLUTION, max_lon=None, auto_range=False)


@pytest.mark.skipif(not pyart.graph.gridmapdisplay_basemap._BASEMAP_AVAILABLE,
                    reason="Basemap is not installed.")
def test_generate_filename_basemap():
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplayBasemap(grid)
    filename = display.generate_filename('reflectivity', 0)
    assert isinstance(filename, str)


@pytest.mark.skipif(not pyart.graph.gridmapdisplay_basemap._BASEMAP_AVAILABLE,
                    reason="Basemap is not installed.")
def test_generate_titles_basemap():
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplayBasemap(grid)

    title = display.generate_longitudinal_level_title('reflectivity', 0)
    assert isinstance(title, str)

    title = display.generate_latitudinal_level_title('reflectivity', 0)
    assert isinstance(title, str)


@pytest.mark.skipif(not pyart.graph.gridmapdisplay_basemap._BASEMAP_AVAILABLE,
                    reason="Basemap is not installed.")
def test_get_basemap():
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplayBasemap(grid)
    basemap = display.get_basemap()


@pytest.mark.skipif(not pyart.graph.gridmapdisplay_basemap._BASEMAP_AVAILABLE,
                    reason="Basemap is not installed.")
def test_error_raising_basemap():
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplayBasemap(grid)

    # no mappable
    pytest.raises(ValueError, display.plot_colorbar)

    # no field
    display.mappables.append(None)  # mock the mappable
    pytest.raises(ValueError, display.plot_colorbar)


if __name__ == "__main__":
    test_gridmapdisplay_simple('figure_gridmapdisplay_simple_basemap.png')
    test_gridmapdisplay_fancy('figure_gridmapdisplay_fancy_basemap.png')
