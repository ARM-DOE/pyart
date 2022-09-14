""" Unit Tests for Py-ART's graph/radar_display.py module. """
# execute this script to create figure_plot_radar_display_*.png files.

# TODO use matplotlib's @image_comparison decorator to compare to file
# in baseline_images directory. Current this test only determines if files can
# be created, not that they are correct.

import warnings
import datetime

import matplotlib.pyplot as plt

import numpy as np
from numpy.testing import assert_almost_equal
import pytest

import pyart


# Top level Figure generating tests
def test_radardisplay_rhi(outfile=None):
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_RHI_FILE)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot('reflectivity_horizontal', 0, ax=ax, mask_outside=True)
    if outfile:
        fig.savefig(outfile)
    plt.close()


def test_radardisplay_ppi(outfile=None):
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    radar.antenna_transition = {'data': np.zeros(radar.nrays)}
    display = pyart.graph.RadarDisplay(radar, shift=(0.1, 0.0))
    fig = plt.figure()
    ax = fig.add_subplot(111)

    gatefilter = pyart.filters.GateFilter(radar)
    display.plot('reflectivity_horizontal', 0, colorbar_flag=True,
                 title="Fancy PPI Plot", mask_outside=True,
                 mask_tuple=('reflectivity_horizontal', -100),
                 gatefilter=gatefilter)
    display.plot_colorbar()
    display.plot_range_rings([10, 20, 30, 40], ax=ax)
    display.plot_labels(['tree'], [(36.68, -97.62)], symbols='k+', ax=ax)
    display.plot_cross_hair(2)
    display.plot_grid_lines()
    display.set_aspect_ratio()
    display.set_limits(ylim=[-50, 50], xlim=[-50, 50])
    if outfile:
        fig.savefig(outfile)
    plt.close()


def test_radardisplay_ray(outfile=None):
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_RHI_FILE)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ray(
        'reflectivity_horizontal', 0, ray_min=15, ray_max=40,
        mask_outside=True, ax=ax, mask_tuple=('reflectivity_horizontal', -10))
    if outfile:
        fig.savefig(outfile)
    plt.close()


def test_radardisplay_vpt(outfile=None):
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    radar.antenna_transition = {'data': np.zeros(radar.nrays)}
    pyart.util.to_vpt(radar)      # hack to make the data a VPT
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    gatefilter = pyart.filters.GateFilter(radar)
    display.plot(
        'reflectivity_horizontal', colorbar_flag=True, mask_outside=True,
        mask_tuple=('reflectivity_horizontal', -100), ax=ax,
        gatefilter=gatefilter)
    if outfile:
        fig.savefig(outfile)
    plt.close()


def test_radardisplay_cr_raster(outfile=None):
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_CR_RASTER_FILE)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_cr_raster(target_range=478., el_limits=[-0.5, 2.5])
    if outfile:
        fig.savefig(outfile)
    plt.close()


def test_radardisplay_vpt_time(outfile=None):
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    pyart.util.to_vpt(radar)      # hack to make the data a VPT
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot('reflectivity_horizontal', colorbar_flag=True,
                 time_axis_flag=True, edges=False,
                 mask_tuple=('reflectivity_horizontal', -100), ax=ax)
    if outfile:
        fig.savefig(outfile)
    plt.close()


def test_radardisplay_azimuth_to_rhi(outfile=None):
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    radar.antenna_transition = {'data': np.zeros(radar.nrays)}
    gatefilter = pyart.filters.GateFilter(radar)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_azimuth_to_rhi(
        'reflectivity_horizontal', 45., colorbar_flag=True, ax=ax,
        mask_tuple=('reflectivity_horizontal', -100), mask_outside=True,
        gatefilter=gatefilter)
    if outfile:
        fig.savefig(outfile)
    plt.close()

# Tests of methods, these tests do not generate figures


def test_radardisplay_init():
    # test that a display object can be created with and without
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    radar.antenna_transition = {'data': np.zeros((40, ))}
    display = pyart.graph.RadarDisplay(radar)
    assert display.antenna_transition is not None
    plt.close()


def test_radardisplay_plot_rhi_reverse():
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_RHI_FILE)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot('reflectivity_horizontal', 0, ax=ax, reverse_xaxis=True)
    plt.close()


def test_radardisplay_plot_azimuth_to_rhi_reverse(outfile=None):
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_azimuth_to_rhi(
        'reflectivity_horizontal', 45., reverse_xaxis=True)
    plt.close()


def test_radardisplay_generate_filename():
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarDisplay(radar, shift=(0.1, 0.0))
    filename = display.generate_filename('test', 0)
    assert filename == 'xsapr-sgp_test_00_20110520105416.png'

    filename = display.generate_filename('test', 0, datetime_format='%Y%m%d')
    assert filename == 'xsapr-sgp_test_00_20110520.png'
    plt.close()


def test_radardisplay_plot_labels_errors():
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarDisplay(radar)
    # len(labels) != len(locations)
    pytest.raises(ValueError, display.plot_labels, ['a', 'a'], [(0, 0)])
    # len(labels) != len(symbols)
    pytest.raises(ValueError, display.plot_labels, ['a'], [(0, 0)],
                  symbols=['r+', 'r+'])
    # unknown scan_type
    display.scan_type = 'foo'
    pytest.raises(ValueError, display.plot, 'fake_field')
    plt.close()


def test_radardisplay_user_specified_labels():
    # test that labels are set when a user specifies them.
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    display._set_ray_title('field', 0, 'foo', ax)
    assert ax.get_title() == 'foo'

    display._label_axes_ppi(('foo', 'bar'), ax)
    assert ax.get_xlabel() == 'foo'
    assert ax.get_ylabel() == 'bar'

    display._label_axes_rhi(('spam', 'eggs'), ax)
    assert ax.get_xlabel() == 'spam'
    assert ax.get_ylabel() == 'eggs'

    display._label_axes_ray(('baz', 'qux'), 'field', ax)
    assert ax.get_xlabel() == 'baz'
    assert ax.get_ylabel() == 'qux'

    display._label_axes_vpt(('nick', 'nock'), False, ax)
    assert ax.get_xlabel() == 'nick'
    assert ax.get_ylabel() == 'nock'
    plt.close()


def test_radardisplay_loc_of_moving_radar():
    radar = pyart.testing.make_empty_ppi_radar(1, 1, 1)

    radar.latitude['data'] = np.array([35, 45])
    radar.longitude['data'] = np.array([75, 85])
    pytest.warns(UserWarning, pyart.graph.RadarDisplay, radar)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        display = pyart.graph.RadarDisplay(radar)
        assert_almost_equal(display.loc, (40, 80), 0)


def test_radardisplay_get_x_z():
    radar = pyart.testing.make_empty_ppi_radar(1, 1, 1)
    display = pyart.graph.RadarDisplay(radar)
    x, z = display._get_x_z(0, False, False)
    assert x.shape == (1, 1)
    assert z.shape == (1, 1)


def test_set_title():
    radar = pyart.testing.make_empty_ppi_radar(1, 1, 1)
    radar.fields['foo'] = {}
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    display._set_az_rhi_title('foo', 0, 'special title', ax)
    assert ax.get_title() == 'special title'


def test_radardisplay_misc():
    # misc methods which are not tested above
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # _set_vpt_title with a title
    display._set_vpt_title('foo_field', 'title_string', ax)
    assert ax.get_title() == 'title_string'

    # _generate_field_name method
    fn = pyart.graph.common.generate_field_name(
        radar, 'reflectivity_horizontal')
    assert fn == 'Equivalent reflectivity factor'

    display.fields['reflectivity_horizontal'].pop('standard_name')
    fn = pyart.graph.common.generate_field_name(
        radar, 'reflectivity_horizontal')
    assert fn == 'Reflectivity'

    display.fields['reflectivity_horizontal'].pop('long_name')
    fn = pyart.graph.common.generate_field_name(
        radar, 'reflectivity_horizontal')
    assert fn == 'Reflectivity horizontal'

    plt.close()


def test_radardisplay_get_colorbar_label():
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
    display = pyart.graph.RadarDisplay(radar)

    # default is to base the label on the standard_name key
    assert (display._get_colorbar_label('reflectivity_horizontal') ==
            'equivalent reflectivity factor (dBZ)')

    # next is to look at the long_name
    del display.fields['reflectivity_horizontal']['standard_name']
    assert (display._get_colorbar_label('reflectivity_horizontal') ==
            'Reflectivity (dBZ)')

    # use the field if standard_name and long_name missing
    del display.fields['reflectivity_horizontal']['long_name']
    print(display._get_colorbar_label('reflectivity_horizontal'))
    assert (display._get_colorbar_label('reflectivity_horizontal') ==
            'reflectivity horizontal (dBZ)')

    # no units if key is missing
    del display.fields['reflectivity_horizontal']['units']
    print(display._get_colorbar_label('reflectivity_horizontal'))
    assert (display._get_colorbar_label('reflectivity_horizontal') ==
            'reflectivity horizontal (?)')
    plt.close()


if __name__ == "__main__":
    test_radardisplay_rhi('figure_radar_display_rhi.png')
    test_radardisplay_ppi('figure_radar_display_ppi.png')
    test_radardisplay_ray('figure_radar_display_ray.png')
    test_radardisplay_vpt('figure_radar_display_vpt.png')
    test_radardisplay_cr_raster('figure_radar_display_cr_raster.png')
