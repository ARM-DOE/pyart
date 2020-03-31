""" Unit Tests for Py-ART's core/radar_spectra.py module. """

import inspect

import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal
import pytest

try:
    import xarray as xr
    _XARRAY_AVAILABLE = True
except ImportError:
    _XARRAY_AVAILABLE = False

import pyart


@pytest.mark.skipif(not _XARRAY_AVAILABLE,
                    reason="Xarray is not installed.")
def test_gate_longitude_latitude():
    radar = pyart.testing.make_empty_spectra_radar(5, 4, 2)
    radar.ds['azimuth'] = xr.DataArray(
        np.array([0, 90, 180, 270, 0]), dims='time')
    radar.ds['elevation'] = xr.DataArray(
        np.array([0, 0, 0, 0, 10]), dims='time')
    radar.ds['range'] = xr.DataArray(
        np.array([5, 15, 25, 35]), dims='range')
    assert radar.gate_longitude.data.shape == (5, 4)
    assert radar.gate_latitude.data.shape == (5, 4)
    assert_almost_equal(radar.gate_longitude.data[0, 0], -97.5, 1)
    assert_almost_equal(radar.gate_latitude.data[0, 0], 36.5, 1)
    # reset and try again with a non-default lat_0/lon_0
    radar.ds.projection.pop('_include_lon_0_lat_0')
    radar.ds.projection['lat_0'] = 20.0
    radar.ds.projection['lon_0'] = 60.0
    radar.init_gate_longitude_latitude()
    assert_almost_equal(radar.gate_longitude.data[0, 0], 60.0, 1)
    assert_almost_equal(radar.gate_latitude.data[0, 0], 20.0, 1)


@pytest.mark.skipif(not _XARRAY_AVAILABLE,
                    reason="Xarray is not installed.")
def test_gate_altitude():
    radar = pyart.testing.make_empty_spectra_radar(5, 4, 2)
    radar.ds['azimuth'] = xr.DataArray(
        np.array([0, 90, 180, 270, 0]), dims='time')
    radar.ds['elevation'] = xr.DataArray(
        np.array([0, 0, 0, 0, 10]), dims='time')
    radar.ds['range'] = xr.DataArray(
        np.array([5, 15, 25, 35]), dims='range')

    assert radar.gate_altitude.values.shape == (5, 4)
    assert_almost_equal(radar.gate_altitude.values[0, 0], 200.0, 1)

    radar.altitude.values[0] = 150.
    radar.init_gate_altitude()
    assert_almost_equal(radar.gate_altitude.values[0, 0], 150.0, 1)


@pytest.mark.skipif(not _XARRAY_AVAILABLE,
                    reason="Xarray is not installed.")
def test_gate_x_y_z():
    radar = pyart.testing.make_empty_spectra_radar(5, 5, 2)
    radar.ds['azimuth'] = xr.DataArray(
        np.array([0, 90, 180, 270, 0]), dims='time')
    radar.ds['elevation'] = xr.DataArray(
        np.array([0, 0, 0, 0, 10]), dims='time')
    radar.ds['range'] = xr.DataArray(
        np.array([5, 15, 25, 35, 45]), dims='range')
    radar.init_gate_x_y_z()
    assert radar.gate_x.values.shape == (5, 5)
    assert_allclose(radar.gate_x.values[0], [0, 0, 0, 0, 0], atol=1e-14)
    assert_allclose(radar.gate_x.values[1], [5, 15, 25, 35, 45], atol=1e-14)
    assert_allclose(radar.gate_x.values[2], [0, 0, 0, 0, 0], atol=1e-5)
    assert_allclose(
        radar.gate_x.values[3], [-5, -15, -25, -35, -45], atol=1e-14)

    assert radar.gate_y.values.shape == (5, 5)
    assert_allclose(radar.gate_y.values[0], [5, 15, 25, 35, 45], atol=1e-14)
    assert_allclose(radar.gate_y.values[1], [0, 0, 0, 0, 0], atol=1e-5)
    assert_allclose(
        radar.gate_y.values[2], [-5, -15, -25, -35, -45], atol=1e-14)
    assert_allclose(radar.gate_y.values[3], [0, 0, 0, 0, 0], atol=1e-6)

    assert radar.gate_z.values.shape == (5, 5)
    z_sweep0 = np.array([1.47e-6, 1.324e-5, 3.679e-5, 7.210e-5, 1.1919e-4])
    assert_allclose(radar.gate_z.values[0], z_sweep0, atol=1e-3)
    assert_allclose(radar.gate_z.values[1], z_sweep0, atol=1e-3)
    assert_allclose(radar.gate_z.values[2], z_sweep0, atol=1e-3)
    assert_allclose(radar.gate_z.values[3], z_sweep0, atol=1e-3)


@pytest.mark.skipif(not _XARRAY_AVAILABLE,
                    reason="Xarray is not installed.")
def test_init_gate_x_y_z():
    radar = pyart.testing.make_empty_spectra_radar(4, 5, 1)
    radar.ds['azimuth'] = xr.DataArray(np.array([0, 90, 180, 270]), dims='time')
    radar.ds['elevation'] = xr.DataArray(np.array([0, 0, 0, 0]), dims='time')
    radar.ds['range'] = xr.DataArray(
        np.array([5, 15, 25, 35, 45]), dims='range')
    # access and check initial gate locations
    assert_allclose(
        radar.gate_x.values[1], [0., 4.3627267, 8.72545, 13.08817, 17.450884],
        atol=1e-14)
    assert_allclose(
        radar.gate_y.values[0], [0., 249.97852, 499.95685, 749.935, 999.9128],
        atol=1e-14)
    z_sweep0 = np.array([0., 3., 6., 9., 13.])
    assert_allclose(radar.gate_z.values[0], z_sweep0, atol=1e-3)

    # change range, gate_x, y, z are not updated
    radar.ds['range'] = xr.DataArray(
        np.array([15, 25, 35, 45, 55]), dims='range')
    assert_allclose(
        radar.gate_x.values[1], [0., 4.3627267, 8.72545, 13.08817, 17.450884],
        atol=1e-14)
    assert_allclose(
        radar.gate_y.values[0], [0., 249.97852, 499.95685, 749.935, 999.9128],
        atol=1e-14)
    z_sweep0 = np.array([0., 3., 6., 9., 13.])
    assert_allclose(radar.gate_z.values[0], z_sweep0, atol=1e-3)

    # call init_gate_x_y_z, now the attributes are updated
    radar.init_gate_x_y_z()
    assert_allclose(radar.gate_x.values[1], [15, 25, 35, 45, 55], atol=1e-14)
    assert_allclose(radar.gate_y.values[0], [15, 25, 35, 45, 55], atol=1e-14)
    z_sweep0 = np.array([1.324e-5, 3.679e-5, 7.210e-5, 1.1919e-4, 1.7805e-4])
    assert_allclose(radar.gate_z.values[0], z_sweep0, atol=1e-3)


@pytest.mark.skipif(not _XARRAY_AVAILABLE,
                    reason="Xarray is not installed.")
def test_rays_per_sweep_attribute():
    radar = pyart.testing.make_target_spectra_radar()
    rays_per_sweep = radar.rays_per_sweep
    assert isinstance(rays_per_sweep, xr.DataArray)
    assert rays_per_sweep.values.shape == (1, )
    assert rays_per_sweep.values[0] == 10


@pytest.mark.skipif(not _XARRAY_AVAILABLE,
                    reason="Xarray is not installed.")
def test_radar_creation():
    radar = pyart.testing.make_target_spectra_radar()
    assert isinstance(radar, pyart.core.RadarSpectra)


@pytest.mark.skipif(not _XARRAY_AVAILABLE,
                    reason="Xarray is not installed.")
def test_iterators():
    radar = pyart.testing.make_empty_spectra_radar(30, 20, 5)

    starts = [0]
    ends = [29]
    starts_ends = [(s, e) for s, e in zip(starts, ends)]

    assert inspect.isgenerator(radar.iter_start())
    assert [s for s in radar.iter_start()] == starts

    assert inspect.isgenerator(radar.iter_end())
    assert [s for s in radar.iter_end()] == ends

    assert inspect.isgenerator(radar.iter_start_end())
    assert [s for s in radar.iter_start_end()] == starts_ends

    assert inspect.isgenerator(radar.iter_slice())
    for s, start, end in zip(radar.iter_slice(), starts, ends):
        assert s.start == start
        assert s.stop == end + 1
        assert s.step is None

    assert inspect.isgenerator(radar.iter_field('spectra'))
    for d in radar.iter_field('spectra'):
        assert d.shape == (30, 20, 5)
        assert d.dtype == np.float64
    pytest.raises(KeyError, radar.iter_field, 'foobar')

    assert inspect.isgenerator(radar.iter_azimuth())
    for d in radar.iter_azimuth():
        assert d.shape == (30, )

    assert inspect.isgenerator(radar.iter_elevation())
    for d in radar.iter_elevation():
        assert d.shape == (30, )


@pytest.mark.skipif(not _XARRAY_AVAILABLE,
                    reason="Xarray is not installed.")
def test_vpt():
    radar = pyart.testing.make_target_spectra_radar()
    vpt = radar.to_vpt()
    assert 'instrument_name' in vpt.metadata
    assert vpt.fields['reflectivity']['data'].shape == (10, 20)
    assert list(vpt.fields.keys()) == [
        'reflectivity', 'velocity', 'spectrum_width',
        'skewness', 'kurtosis']
    assert vpt.range['data'].shape == (20,)
    assert vpt.time['data'].shape == (10,)
