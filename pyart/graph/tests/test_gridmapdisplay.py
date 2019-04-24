""" Unit Tests for Py-ART's graph/gridmapdisplay.py module. """
# execute this script to create figure_gridmapdisplay_*.png files.

# TODO use matplotlib's @image_comparison decorator to compare to file
# in baseline_images directory. Currently this test only determines if files
# can be created, not that they are correct.

import matplotlib.pyplot as plt
import pytest

import cartopy.crs as ccrs
import pyart

@pytest.mark.skipif(not pyart.graph.gridmapdisplay._XARRAY_AVAILABLE,
                    reason='Xarray is not installed')
def test_gridmapdisplay_simple(outfile=None):
    # test basic GridMapDisplat functionally.
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_grid('reflectivity', vmin=-5, vmax=35, ax=ax)
    if outfile:
        fig.savefig(outfile)

@pytest.mark.skipif(not pyart.graph.gridmapdisplay._XARRAY_AVAILABLE,
                    reason='Xarray is not installed')
def test_gridmapdisplay_fancy(outfile=None):
    # test a bunch of GridMapDisplay functionaliy
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure(1, figsize=[25., 17.])
    projection = ccrs.Mercator()

    ax1 = plt.subplot(331, projection=projection)
    display.plot_grid('reflectivity', vmin=-5., vmax=35., mask_outside=True,
                      axislabels=('foo', 'bar'), axislabels_flag=True,
                      title='Special title')
    display.plot_crosshairs(color='b')

    ax2 = plt.subplot(332, projection=projection)
    display.plot_grid('reflectivity', axislabels_flag=True)

    ax3 = plt.subplot(333)
    display.plot_colorbar()

    ax4 = plt.subplot(334)
    display.plot_latitude_slice('reflectivity', mask_outside=True)

    ax5 = plt.subplot(335)
    display.plot_latitude_slice('reflectivity', title='Lat title')

    ax6 = plt.subplot(336)
    grid.fields['reflectivity']['valid_min'] = 0
    grid.fields['reflectivity']['valid_max'] = 30
    display.plot_longitude_slice('reflectivity', mask_outside=True)

    ax7 = plt.subplot(337)
    display.plot_longitude_slice('reflectivity', title='Lon title')

    ax8 = plt.subplot(338)
    del display.grid.fields['reflectivity']['long_name']
    display.plot_colorbar()

    if outfile:
        fig.savefig(outfile)

@pytest.mark.skipif(not pyart.graph.gridmapdisplay._XARRAY_AVAILABLE,
                    reason='Xarray is not installed')
def test_generate_filename():
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplay(grid)
    filename = display.generate_filename('reflectivity', 0)
    assert isinstance(filename, str)

@pytest.mark.skipif(not pyart.graph.gridmapdisplay._XARRAY_AVAILABLE,
                    reason='Xarray is not installed')
def test_generate_titles():
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplay(grid)

    title = display.generate_longitudinal_level_title('reflectivity', 0)
    assert isinstance(title, str)

    title = display.generate_latitudinal_level_title('reflectivity', 0)
    assert isinstance(title, str)

@pytest.mark.skipif(not pyart.graph.gridmapdisplay._XARRAY_AVAILABLE,
                    reason='Xarray is not installed')
def test_error_raising():
    grid = pyart.testing.make_target_grid()
    display = pyart.graph.GridMapDisplay(grid)

    # no mappable
    pytest.raises(ValueError, display.plot_colorbar)

    # no field
    display.mappables.append(None) # mock the mappable
    pytest.raises(ValueError, display.plot_colorbar)

if __name__ == '__main__':
    test_gridmapdisplay_simple('figure_grid_mapdisplay_simple.png')
    test_gridmapdisplay_fancy('figure_gridmapdisplay_fancy.png')
