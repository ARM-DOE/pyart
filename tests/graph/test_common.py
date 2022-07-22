""" Unit Tests for Py-ART's graph/common.py module. """

import datetime

from numpy.testing import assert_almost_equal

import matplotlib
import matplotlib.pyplot as plt

import pyart
from pyart.graph import common


def test_parse_ax():
    ax = common.parse_ax(None)
    assert isinstance(ax, matplotlib.axes.Axes)

    ax1 = plt.gca()
    ax2 = common.parse_ax(ax1)
    assert ax1 == ax2


def test_parse_ax_fig():
    ax, fig = common.parse_ax_fig(None, None)
    assert isinstance(ax, matplotlib.axes.Axes)
    assert isinstance(fig, matplotlib.figure.Figure)

    ax1 = plt.gca()
    fig1 = plt.gcf()
    ax2, fig2 = common.parse_ax_fig(ax1, fig1)
    assert ax1 == ax2
    assert fig1 == fig2


def test_parse_cmap():
    assert common.parse_cmap('jet', 'foo') == 'jet'
    assert common.parse_cmap(None, 'reflectivity') == 'pyart_HomeyerRainbow'


def test_parse_vmin_vmax():
    radar = pyart.testing.make_empty_ppi_radar(1, 1, 1)
    radar.fields['foo'] = {}

    vmin, vmax = common.parse_vmin_vmax(radar, 'foo', None, None)
    assert vmin is None
    assert vmax is None

    vmin, vmax = common.parse_vmin_vmax(radar, 'foo', 10, 20)
    assert vmin == 10
    assert vmax == 20

    radar.fields['foo']['valid_min'] = 30
    radar.fields['foo']['valid_max'] = 40
    vmin, vmax = common.parse_vmin_vmax(radar, 'foo', None, None)
    assert vmin == 30
    assert vmax == 40


def test_parse_lon_lat():
    grid_shape = (1, 1, 1)
    grid_limits = ((0, 1), (0, 1), (0, 1))
    grid = pyart.testing.make_empty_grid(grid_shape, grid_limits)

    lon, lat = common.parse_lon_lat(grid, None, None)
    assert_almost_equal(lon, -98.1, 2)
    assert_almost_equal(lat, 36.74, 2)

    lon, lat = common.parse_lon_lat(grid, -12.34, 56.78)
    assert_almost_equal(lon, -12.34, 2)
    assert_almost_equal(lat, 56.78, 2)


def test_generate_colorbar_label():
    label = common.generate_colorbar_label('special_name', 'km')
    assert label == 'special name (km)'


def test_generate_field_name():
    radar = pyart.testing.make_empty_ppi_radar(1, 1, 1)
    radar.fields['foo'] = {}

    field_name = common.generate_field_name(radar, 'foo')
    assert field_name == 'Foo'

    radar.fields['foo']['long_name'] = 'special_name'
    field_name = common.generate_field_name(radar, 'foo')
    assert field_name == 'Special name'

    radar.fields['foo']['standard_name'] = 'other_name'
    field_name = common.generate_field_name(radar, 'foo')
    assert field_name == 'Other name'


def test_generate_radar_name():
    radar = pyart.testing.make_empty_ppi_radar(1, 1, 1)

    radar.metadata['instrument_name'] = 'foobar'
    name = common.generate_radar_name(radar)
    assert name == 'foobar'

    del radar.metadata['instrument_name']
    name = common.generate_radar_name(radar)
    assert name == ''


def test_generate_grid_name():
    grid_shape = (1, 1, 1)
    grid_limits = ((0, 1), (0, 1), (0, 1))
    grid = pyart.testing.make_empty_grid(grid_shape, grid_limits)

    grid.metadata['instrument_name'] = 'foobar'
    name = common.generate_grid_name(grid)
    assert name == 'foobar'

    del grid.metadata['instrument_name']
    name = common.generate_grid_name(grid)
    assert name == ''


def test_generate_radar_time_begin():
    radar = pyart.testing.make_empty_ppi_radar(1, 1, 1)
    time = common.generate_radar_time_begin(radar)
    assert time == datetime.datetime(1989, 1, 1, 0, 0, 1)


def test_generate_radar_time_sweep():
    radar = pyart.testing.make_empty_ppi_radar(1, 1, 2)

    time = common.generate_radar_time_sweep(radar, 0)
    assert time == datetime.datetime(1989, 1, 1, 0, 0, 1)

    time = common.generate_radar_time_sweep(radar, 1)
    assert time == datetime.datetime(1989, 1, 1, 0, 0, 2)


def test_generate_grid_time_begin():
    grid_shape = (1, 1, 1)
    grid_limits = ((0, 1), (0, 1), (0, 1))
    grid = pyart.testing.make_empty_grid(grid_shape, grid_limits)

    time = common.generate_grid_time_begin(grid)
    assert time == datetime.datetime(2000, 1, 1, 0, 0, 0)

    del grid.time['calendar']
    time = common.generate_grid_time_begin(grid)
    assert time == datetime.datetime(2000, 1, 1, 0, 0, 0)


def test_generate_filename():
    radar = pyart.testing.make_empty_ppi_radar(1, 1, 2)
    filename = common.generate_filename(radar, 'foobar', 0)
    assert filename == 'fake_radar_foobar_00_19890101000001.png'

    filename = common.generate_filename(radar, 'foobar', 1)
    assert filename == 'fake_radar_foobar_01_19890101000001.png'

    filename = common.generate_filename(radar, 'foobar', 1, datetime_format='%Y%m%d')
    assert filename == 'fake_radar_foobar_01_19890101.png'

    filename = common.generate_filename(radar, 'foobar', 1, use_sweep_time=True)
    assert filename == 'fake_radar_foobar_01_19890101000002.png'


def test_generate_grid_filename():
    grid_shape = (1, 1, 1)
    grid_limits = ((0, 1), (0, 1), (0, 1))
    grid = pyart.testing.make_empty_grid(grid_shape, grid_limits)
    grid.metadata['instrument_name'] = 'foo'

    filename = common.generate_grid_filename(grid, 'bar', 0)
    assert filename == 'foo_bar_00_20000101000000.png'


def test_generate_title():
    radar = pyart.testing.make_empty_ppi_radar(1, 1, 2)
    radar.fields['foo'] = {}

    title = common.generate_title(radar, 'foo', 0)
    assert title == 'fake_radar 0.8 Deg. 1989-01-01T00:00:01Z \nFoo'

    title = common.generate_title(radar, 'foo', 1)
    assert title == 'fake_radar 0.8 Deg. 1989-01-01T00:00:02Z \nFoo'

    title = common.generate_title(radar, 'foo', 1, datetime_format='%Y-%m-%d')
    assert title == 'fake_radar 0.8 Deg. 1989-01-01 \nFoo'

    title = common.generate_title(radar, 'foo', 1, use_sweep_time=False)
    assert title == 'fake_radar 0.8 Deg. 1989-01-01T00:00:01Z \nFoo'


def test_generate_grid_title():
    grid_shape = (1, 1, 1)
    grid_limits = ((0, 1), (0, 1), (0, 1))
    grid = pyart.testing.make_empty_grid(grid_shape, grid_limits)
    grid.metadata['instrument_name'] = 'bar'
    grid.fields['foo'] = {}

    title = common.generate_grid_title(grid, 'foo', 0)
    assert title == 'bar 0.0 km 2000-01-01T00:00:00Z \nFoo'


def test_generate_longitudinal_level_title():
    grid_shape = (1, 1, 1)
    grid_limits = ((0, 1), (0, 1), (0, 1))
    grid = pyart.testing.make_empty_grid(grid_shape, grid_limits)
    grid.metadata['instrument_name'] = 'bar'
    grid.fields['foo'] = {}

    title = common.generate_longitudinal_level_title(grid, 'foo', 0)
    assert title == (
        'bar 0.0 km east of origin 2000-01-01T00:00:00Z \nFoo')

    grid.x['data'][0] = -1000.
    title = common.generate_longitudinal_level_title(grid, 'foo', 0)
    assert title == (
        'bar 1.0 km west of origin 2000-01-01T00:00:00Z \nFoo')


def test_generate_latitudinal_level_title():
    grid_shape = (1, 1, 1)
    grid_limits = ((0, 1), (0, 1), (0, 1))
    grid = pyart.testing.make_empty_grid(grid_shape, grid_limits)
    grid.metadata['instrument_name'] = 'bar'
    grid.fields['foo'] = {}

    title = common.generate_latitudinal_level_title(grid, 'foo', 0)
    assert title == (
        'bar 0.0 km north of origin 2000-01-01T00:00:00Z \nFoo')

    grid.y['data'][0] = -1000.
    title = common.generate_latitudinal_level_title(grid, 'foo', 0)
    assert title == (
        'bar 1.0 km south of origin 2000-01-01T00:00:00Z \nFoo')


def test_generate_vpt_title():
    radar = pyart.testing.make_empty_ppi_radar(1, 1, 1)
    radar.fields['foo'] = {}

    title = common.generate_vpt_title(radar, 'foo')
    assert title == 'fake_radar 1989-01-01T00:00:01Z \nFoo'


def test_generate_ray_title():
    radar = pyart.testing.make_empty_ppi_radar(1, 1, 1)
    radar.fields['foo'] = {}

    title = common.generate_ray_title(radar, 'foo', 0)
    assert title == ('fake_radar 1989-01-01T00:00:01Z\n' +
                     'Ray: 0  Elevation: 0.8 Azimuth: 0.0\nFoo')


def test_generate_az_rhi_title():
    radar = pyart.testing.make_empty_ppi_radar(1, 1, 1)
    radar.fields['foo'] = {}

    title = common.generate_az_rhi_title(radar, 'foo', 0)
    assert title == (
        'fake_radar 1989-01-01T00:00:01Z \nAzimuth: 0.0 deg\nFoo')


def test_set_limits():
    common.set_limits(xlim=(-20, 20), ylim=(-10, 10))
    ax = plt.gca()
    assert ax.get_ylim() == (-10, 10)
    assert ax.get_xlim() == (-20, 20)
