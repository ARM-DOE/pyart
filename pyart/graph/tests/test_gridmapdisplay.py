""" Unit Tests for Py-ART's graph/gridmapdisplay.py module. """
# execute this script to create figure_gridmapdisplay_*.png files.

# TODO use matplotlib's @image_comparison decorator to compare to file
# in baseline_images directory. Current this test only determines if files can
# be created, not that they are correct.

import matplotlib.pyplot as plt
import pyart
from numpy.testing import assert_raises
from numpy.testing.decorators import skipif

RESOLUTION = 'c'    # crude resolution to speed up tests


@skipif('GridMapDisplay' not in dir(pyart.graph))
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


@skipif('GridMapDisplay' not in dir(pyart.graph))
def test_gridmapdisplay_fancy(outfile=None):
    # test a bunch of GridMapDisplay functionality
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure()

    ax = fig.add_subplot(321)
    display.plot_basemap(ax=ax, resolution=RESOLUTION)
    display.plot_grid('reflectivity', vmin=-5., vmax=35)
    display.plot_crosshairs(line_style='b--')

    ax = fig.add_subplot(322)
    display.plot_colorbar()

    ax = fig.add_subplot(323)
    display.plot_latitude_slice('reflectivity')  # no vmin/vmax

    ax = fig.add_subplot(325)
    grid.fields['reflectivity']['valid_min'] = 0.
    grid.fields['reflectivity']['valid_max'] = 30.
    display.plot_longitude_slice('reflectivity')

    ax = fig.add_subplot(326)
    del display.grid.fields['reflectivity']['long_name']
    display.plot_colorbar()

    if outfile:
        fig.savefig(outfile)


@skipif('GridMapDisplay' not in dir(pyart.graph))
def test_error_raising():
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplay(grid)
    # no basemap
    assert_raises(ValueError, display.plot_grid, 'reflectivity')

    # no mappable
    assert_raises(ValueError, display.plot_colorbar)

    # no field
    display.mappables.append(None)  # mock the mappable
    assert_raises(ValueError, display.plot_colorbar)


if __name__ == "__main__":
    test_gridmapdisplay_simple('figure_gridmapdisplay_simple.png')
    test_gridmapdisplay_fancy('figure_gridmapdisplay_fancy.png')
