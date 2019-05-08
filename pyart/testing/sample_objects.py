"""
pyart.testing.sample_objects
============================

Functions for creating sample Radar and Grid objects.

.. autosummary::
    :toctree: generated/

    make_empty_ppi_radar
    make_target_radar
    make_velocity_aliased_radar
    make_single_ray_radar
    make_empty_grid
    make_target_grid
    make_storm_grid
    make_normal_storm

"""

import numpy as np

from .sample_files import _EXAMPLE_RAYS_FILE
from ..config import get_metadata
from ..core.radar import Radar
from ..core.grid import Grid


def make_empty_ppi_radar(ngates, rays_per_sweep, nsweeps):
    """
    Return an Radar object, representing a PPI scan.

    Parameters
    ----------
    ngates : int
        Number of gates per ray.
    rays_per_sweep : int
        Number of rays in each PPI sweep.
    nsweeps : int
        Number of sweeps.

    Returns
    -------
    radar : Radar
        Radar object with no fields, other parameters are set to default
        values.

    """
    nrays = rays_per_sweep * nsweeps

    time = get_metadata('time')
    _range = get_metadata('range')
    latitude = get_metadata('latitude')
    longitude = get_metadata('longitude')
    altitude = get_metadata('altitude')
    sweep_number = get_metadata('sweep_number')
    sweep_mode = get_metadata('sweep_mode')
    fixed_angle = get_metadata('fixed_angle')
    sweep_start_ray_index = get_metadata('sweep_start_ray_index')
    sweep_end_ray_index = get_metadata('sweep_end_ray_index')
    azimuth = get_metadata('azimuth')
    elevation = get_metadata('elevation')

    fields = {}
    scan_type = 'ppi'
    metadata = {'instrument_name': 'fake_radar'}

    time['data'] = np.arange(nrays, dtype='float64')
    time['units'] = 'seconds since 1989-01-01T00:00:01Z'
    _range['data'] = np.linspace(0, 1000, ngates).astype('float32')

    latitude['data'] = np.array([36.5], dtype='float64')
    longitude['data'] = np.array([-97.5], dtype='float64')
    altitude['data'] = np.array([200], dtype='float64')

    sweep_number['data'] = np.arange(nsweeps, dtype='int32')
    sweep_mode['data'] = np.array(['azimuth_surveillance'] * nsweeps)
    fixed_angle['data'] = np.array([0.75] * nsweeps, dtype='float32')
    sweep_start_ray_index['data'] = np.arange(0, nrays, rays_per_sweep,
                                              dtype='int32')
    sweep_end_ray_index['data'] = np.arange(rays_per_sweep - 1, nrays,
                                            rays_per_sweep, dtype='int32')

    azimuth['data'] = np.arange(nrays, dtype='float32')
    elevation['data'] = np.array([0.75] * nrays, dtype='float32')

    return Radar(time, _range, fields, metadata, scan_type,
                 latitude, longitude, altitude,
                 sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
                 sweep_end_ray_index,
                 azimuth, elevation,
                 instrument_parameters=None)


def make_empty_rhi_radar(ngates, rays_per_sweep, nsweeps):
    """
    Return an Radar object, representing a RHI scan.

    Parameters
    ----------
    ngates : int
        Number of gates per ray.
    rays_per_sweep : int
        Number of rays in each PPI sweep.
    nsweeps : int
        Number of sweeps.

    Returns
    -------
    radar : Radar
        Radar object with no fields, other parameters are set to default
        values.

    """
    radar = make_empty_ppi_radar(ngates, rays_per_sweep, nsweeps)
    radar.scan_type = 'rhi'
    nrays = rays_per_sweep * nsweeps
    radar.sweep_mode['data'] = np.array(['rhi'] * nsweeps)
    radar.elevation['data'] = np.arange(nrays, dtype='float32')
    radar.azimuth['data'] = np.array([0.75] * nrays, dtype='float32')
    return radar


def make_target_radar():
    """
    Return a PPI radar with a target like reflectivity field.
    """
    radar = make_empty_ppi_radar(50, 360, 1)
    fields = {
        'reflectivity': get_metadata('reflectivity')}
    fdata = np.zeros((360, 50), dtype='float32')
    fdata[:, 0:10] = 0.
    fdata[:, 10:20] = 10.
    fdata[:, 20:30] = 20.
    fdata[:, 30:40] = 30.
    fdata[:, 40:50] = 40.
    fields['reflectivity']['data'] = fdata
    radar.fields = fields
    return radar


def make_velocity_aliased_radar(alias=True):
    """
    Return a PPI radar with a target like reflectivity field.

    Set alias to False to return a de-aliased radar.
    """
    radar = make_empty_ppi_radar(50, 360, 1)
    radar.range['meters_between_gates'] = 1.0
    radar.range['meters_to_center_of_first_gate'] = 1.0
    radar.instrument_parameters = {
        'nyquist_velocity': {'data': np.array([10.0] * 360)}}

    fields = {
        'reflectivity': get_metadata('reflectivity'),
        'velocity': get_metadata('velocity')}

    # fake reflectivity data, all zero reflectivity
    fdata = np.zeros((360, 50), dtype='float32')
    fields['reflectivity']['data'] = fdata

    # fake velocity data, all zeros except a wind burst on at ~13 degrees.
    # burst is partially aliased.
    vdata = np.zeros((360 * 1, 50), dtype='float32')

    for i, idx in enumerate(range(13, -1, -1)):
        vdata[i, idx:idx + i + 1] = np.arange(0.5, 0.5 + i * 1. + 0.001)
    vdata[:14, 14:27] = vdata[:14, 12::-1]  # left/right flip
    vdata[14:27] = vdata[12::-1, :]         # top/bottom flip
    aliased = np.where(vdata > 10.0)
    if alias:
        vdata[aliased] += -20.
    fields['velocity']['data'] = vdata

    radar.fields = fields
    return radar


def make_velocity_aliased_rhi_radar(alias=True):
    """
    Return a RHI radar with a target like reflectivity field.

    Set alias to False to return a de-aliased radar.
    """
    radar = make_empty_rhi_radar(50, 180, 1)
    radar.range['meters_between_gates'] = 1.0
    radar.range['meters_to_center_of_first_gate'] = 1.0
    radar.instrument_parameters = {
        'nyquist_velocity': {'data': np.array([10.0] * 360)}}

    fields = {
        'reflectivity': get_metadata('reflectivity'),
        'velocity': get_metadata('velocity')}

    # fake reflectivity data, all zero reflectivity
    fdata = np.zeros((180, 50), dtype='float32')
    fields['reflectivity']['data'] = fdata

    # fake velocity data, all zeros except a wind burst on at ~13 degrees.
    # burst is partially aliased.
    vdata = np.zeros((180 * 1, 50), dtype='float32')

    for i, idx in enumerate(range(13, -1, -1)):
        vdata[i, idx:idx + i + 1] = np.arange(0.5, 0.5 + i * 1. + 0.001)
    vdata[:14, 14:27] = vdata[:14, 12::-1]  # left/right flip
    vdata[14:27] = vdata[12::-1, :]         # top/bottom flip
    aliased = np.where(vdata > 10.0)
    if alias:
        vdata[aliased] += -20.
    fields['velocity']['data'] = vdata

    radar.fields = fields
    return radar


def make_single_ray_radar():
    """
    Return a PPI radar with a single ray taken from a ARM C-SAPR Radar.

    Radar object returned has 'reflectivity_horizontal',
    'norm_coherent_power', 'copol_coeff', 'dp_phase_shift', 'diff_phase', and
    'differential_reflectivity' fields with no metadata but a 'data' key.
    This radar is used for unit tests in correct modules.

    """
    radar = make_empty_ppi_radar(983, 1, 1)
    radar.range['data'] = 117.8784 + np.arange(983) * 119.91698
    f = np.load(_EXAMPLE_RAYS_FILE)
    for field_name in f:
        radar.fields[field_name] = {'data': f[field_name]}
    f.close()
    return radar


def make_empty_grid(grid_shape, grid_limits):
    """
    Make an empty grid object without any fields or metadata.

    Parameters
    ----------
    grid_shape : 3-tuple of floats
        Number of points in the grid (z, y, x).
    grid_limits : 3-tuple of 2-tuples
        Minimum and maximum grid location (inclusive) in meters for the
        z, y, x coordinates.

    Returns
    -------
    grid : Grid
        Empty Grid object, centered near the ARM SGP site (Oklahoma).

    """
    time = get_metadata('grid_time')
    time['data'] = np.array([0.0])
    time['units'] = 'seconds since 2000-01-01T00:00:00Z'

    # grid coordinate dictionaries
    nz, ny, nx = grid_shape
    (z0, z1), (y0, y1), (x0, x1) = grid_limits

    x = get_metadata('x')
    x['data'] = np.linspace(x0, x1, nx)

    y = get_metadata('y')
    y['data'] = np.linspace(y0, y1, ny)

    z = get_metadata('z')
    z['data'] = np.linspace(z0, z1, nz)

    origin_altitude = get_metadata('origin_altitude')
    origin_altitude['data'] = np.array([300.])

    origin_latitude = get_metadata('origin_latitude')
    origin_latitude['data'] = np.array([36.74])

    origin_longitude = get_metadata('origin_longitude')
    origin_longitude['data'] = np.array([-98.1])

    fields = {}
    metadata = {}

    radar_latitude = get_metadata('radar_latitude')
    radar_latitude['data'] = np.array([36.74])

    radar_longitude = get_metadata('radar_longitude')
    radar_longitude['data'] = np.array([-98.1])

    radar_altitude = get_metadata('radar_altitude')
    radar_altitude['data'] = np.array([300.])

    radar_time = get_metadata('radar_time')
    radar_time['data'] = np.array([0.0])
    radar_time['units'] = 'seconds since 2000-01-01T00:00:00Z'

    radar_name = get_metadata('radar_name')
    radar_name['data'] = np.array(['ExampleRadar'])

    return Grid(time, fields, metadata,
                origin_latitude, origin_longitude, origin_altitude, x, y, z,
                radar_latitude=radar_latitude, radar_longitude=radar_longitude,
                radar_altitude=radar_altitude, radar_time=radar_time,
                radar_name=radar_name)


def make_target_grid():
    """
    Make a sample Grid with a rectangular target.
    """
    grid_shape = (2, 400, 320)
    grid_limits = ((0, 500), (-400000, 400000), (-300000, 300000))
    grid = make_empty_grid(grid_shape, grid_limits)

    fdata = np.zeros((2, 400, 320), dtype='float32')
    fdata[:, 50:-50, 40:-40] = 10.
    fdata[:, 100:-100, 80:-80] = 20.
    fdata[:, 150:-150, 120:-120] = 30.
    fdata[1] += 5
    rdic = {
        'data': fdata,
        'long_name': 'reflectivity',
        'units': 'dBz'}
    grid.fields = {'reflectivity': rdic}
    return grid


def make_storm_grid():
    """
    Make a sample Grid with a rectangular storm target.
    """
    grid_shape = (2, 50, 40)
    grid_limits = ((0, 500), (-400000, 400000), (-300000, 300000))
    grid = make_empty_grid(grid_shape, grid_limits)

    fdata = np.ma.zeros((2, 50, 40), dtype='float32')
    fdata[:] = np.ma.masked
    fdata[:, 5:-5, 4:-4] = 20.
    fdata[:, 10:-10, 8:-8] = 30.
    fdata[:, 15:-15, 12:-12] = 60.
    fdata[1] += 5
    rdic = {
        'data': fdata,
        'long_name': 'reflectivity',
        'units': 'dBz'}
    grid.fields = {'reflectivity': rdic}
    return grid


def make_normal_storm(sigma, mu):
    """
    Make a sample Grid with a gaussian storm target.
    """
    test_grid = make_empty_grid(
        [1, 101, 101], [(1, 1), (-50, 50), (-50, 50)])
    x = test_grid.x['data']
    y = test_grid.y['data']
    z = test_grid.z['data']
    zg, yg, xg = np.meshgrid(z, y, x, indexing='ij')
    r = np.sqrt((xg - mu[0])**2 + (yg - mu[1])**2)
    term1 = 1.0 / (sigma * np.sqrt(2.0 * np.pi))
    term2 = -1.0 * (r**2 / (2.0 * sigma**2))
    data = term1 * np.exp(term2)
    rdic = {
        'data': data,
        'long_name': 'reflectivity',
        'units': 'dBz'}
    test_grid.fields.update({'reflectivity': rdic})
    return test_grid
