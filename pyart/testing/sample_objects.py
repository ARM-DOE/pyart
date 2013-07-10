"""
pyart.testing.sample_objects
============================

Functions for creating sample Radar and Grid objects.

.. autosummary::
    :toctree: generated/

    make_empty_ppi_radar
    make_target_radar

"""

import numpy as np

from .sample_files import _EXAMPLE_RAYS_FILE
from ..io.common import get_metadata
from ..io.radar import Radar


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


def make_target_radar():
    """
    Return a PPI radar with a target like reflectivity_horizontal field.
    """
    radar = make_empty_ppi_radar(50, 360, 1)
    fields = {
        'reflectivity_horizontal': get_metadata('reflectivity_horizontal')}
    fdata = np.zeros((360, 50), dtype='float32')
    fdata[:, 0:10] = 0.
    fdata[:, 10:20] = 10.
    fdata[:, 20:30] = 20.
    fdata[:, 30:40] = 30.
    fdata[:, 40:50] = 40.
    fields['reflectivity_horizontal']['data'] = fdata
    radar.fields = fields
    return radar


def make_single_ray_radar():
    """
    Return a PPI radar with a single ray taken from a ARM C-SAPR Radar

    Radar object returned has 'reflectivity_horizontal',
    'norm_coherent_power', 'copol_coeff', 'dp_phase_shift', and 'diff_phase'
    fields with no metadata but a 'data' key.  This radar is used for unit
    tests in correct modules.

    """
    radar = make_empty_ppi_radar(983, 1, 1)
    radar.range['data'] = 117.8784 + np.arange(983) * 119.91698
    f = np.load(_EXAMPLE_RAYS_FILE)
    for field_name in f:
        radar.fields[field_name] = {'data': f[field_name]}
    f.close()
    return radar
