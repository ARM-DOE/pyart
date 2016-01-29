""" Unit Tests for Py-ART's graph/gridmapdisplay.py module. """
# execute this script to create figure_gridmapdisplay_*.png files.

# TODO use matplotlib's @image_comparison decorator to compare to file
# in baseline_images directory. Current this test only determines if files can
# be created, not that they are correct.

import warnings

import matplotlib.pyplot as plt
import pyart
from pyart.exceptions import DeprecatedAttribute
from numpy.testing import assert_raises
from numpy.testing.decorators import skipif

RESOLUTION = 'c'    # crude resolution to speed up tests


@skipif(not pyart.graph.gridmapdisplay._BASEMAP_AVAILABLE)
def test_gridmapdisplay_simple(outfile=None):
    # test basic GridMapDisplay functionality
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_basemap(ax=ax, resolution=RESOLUTION)
    display.plot_grid('reflectivity', vmin=-5., vmax=35)
    if outfile:
        fig.savefig(outfile)


@skipif(not pyart.graph.gridmapdisplay._BASEMAP_AVAILABLE)
def test_gridmapdisplay_fancy(outfile=None):
    # test a bunch of GridMapDisplay functionality
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplay(grid)
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


@skipif(not pyart.graph.gridmapdisplay._BASEMAP_AVAILABLE)
def test_plot_basemap_not_using_corners(outfile=None):
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_basemap(
        ax=ax, resolution=RESOLUTION, max_lon=None, auto_range=False)


@skipif(not pyart.graph.gridmapdisplay._BASEMAP_AVAILABLE)
def test_generate_filename():
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplay(grid)
    filename = display.generate_filename('reflectivity', 0)
    assert isinstance(filename, str)


@skipif(not pyart.graph.gridmapdisplay._BASEMAP_AVAILABLE)
def test_generate_titles():
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplay(grid)

    title = display.generate_longitudinal_level_title('reflectivity', 0)
    assert isinstance(title, str)

    title = display.generate_latitudinal_level_title('reflectivity', 0)
    assert isinstance(title, str)


@skipif(not pyart.graph.gridmapdisplay._BASEMAP_AVAILABLE)
def test_get_basemap():
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplay(grid)
    basemap = display.get_basemap()


@skipif(not pyart.graph.gridmapdisplay._BASEMAP_AVAILABLE)
def test_error_raising():
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplay(grid)

    # no mappable
    assert_raises(ValueError, display.plot_colorbar)

    # no field
    display.mappables.append(None)  # mock the mappable
    assert_raises(ValueError, display.plot_colorbar)


@skipif(not pyart.graph.gridmapdisplay._BASEMAP_AVAILABLE)
def test_deprecated_attributes():

    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplay(grid)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=DeprecatedAttribute)

        assert display.grid_lons.ndim == 2
        assert display.grid_lats.ndim == 2
        assert callable(display.proj)


if __name__ == "__main__":
    test_gridmapdisplay_simple('figure_gridmapdisplay_simple.png')
    test_gridmapdisplay_fancy('figure_gridmapdisplay_fancy.png')
