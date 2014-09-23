"""
pyart.correct.unwrap
====================

Dealias using multidimensional phase unwrapping algorithms.

.. autosummary::
    :toctree: generated/

    dealias_unwrap_phase

"""

from __future__ import print_function

import numpy as np

from ..config import get_field_name, get_metadata

from ._unwrap_1d import unwrap_1d
from ._unwrap_2d import unwrap_2d
from ._unwrap_3d import unwrap_3d

# TODO
# add choice to change sweep periodicity and auto detect for ppi/rhi
# test for not allowed unwrap_unit choices?
# unit tests
# clean up correct namespace


def dealias_unwrap_phase(radar, nyquist_vel=None, unwrap_unit='sweep',
                         vel_field=None, corr_vel_field=None):
    """
    Dealias Doppler velocities using multi-dimensional phase unwrapping.

    Parameters
    ----------
    radar : Radar
        Radar object containing Doppler velocities to dealias.
    unwrap_unit : {'ray', 'sweep', 'volume', None}
        Unit to unwrap independently.  'ray' will unwrap each ray
        individually, 'sweep' each sweep, and 'volume' will unwrap the entire
        volume in a single pass.  None will determine the larges allowed
        unit from the structure of the radar.  'sweep' often gives superior
        results when the lower sweeps of the radar volume are contaminated by
        clutter.
    nyquist_velocity : float, optional
        Nyquist velocity in unit identical to those stored in the radar's
        velocity field.  None will attempt to determine this value from the
        instrument_parameters attribute.
    vel_field : str
        Field in radar to use as the Doppler velocities during dealiasing.
        None will use the default field name from the Py-ART configuration
        file.
    corr_vel_field : str
        Name to use for the dealiased Doppler velocity field metadata.  None
        will use the default field name from the Py-ART configuration file.

    Returns
    -------

    References
    ----------
    .. [1] Miguel Arevallilo Herraez, David R. Burton, Michael J. Lalor,
           and Munther A. Gdeisat, "Fast two-dimensional phase-unwrapping
           algorithm based on sorting by reliability following a noncontinuous
           path", Journal Applied Optics, Vol. 41, No. 35 (2002) 7437,
    .. [2] Abdul-Rahman, H., Gdeisat, M., Burton, D., & Lalor, M., "Fast
           three-dimensional phase-unwrapping algorithm based on sorting by
           reliability following a non-continuous path. In W. Osten,
           C. Gorecki, & E. L. Novak (Eds.), Optical Metrology (2005) 32--40,
           International Society for Optics and Photonics.

    """
    # parse the field parameters.
    if vel_field is None:
        vel_field = get_field_name('velocity')
    if corr_vel_field is None:
        corr_vel_field = get_field_name('corrected_velocity')

    # find the nyquist velocity if not specified.
    if nyquist_vel is None:
        if (radar.instrument_parameters is None) or (
                'nyquist_velocity' not in radar.instrument_parameters):
            message = ('Nyquist velocity not specified in radar object, '
                       'provide this value explicitly in the function call.')
            raise ValueError(message)
        nyquist_vel = radar.instrument_parameters[
            'nyquist_velocity']['data'][0]

    # find unwrap_unit if not specified.
    if unwrap_unit is None:
        if _is_radar_sequential(radar):
            if _is_radar_cubic(radar) and _is_radar_sweep_aligned(radar):
                unwrap_unit = 'volume'
            else:
                unwrap_unit = 'sweep'
        else:
            unwrap_unit = 'ray'
        print("Unwrapping using unwrap_unit:", unwrap_unit)

    # perform dealiasing
    if unwrap_unit == 'ray':
        data = _dealias_unwrap_1d(radar, vel_field, nyquist_vel)
    elif unwrap_unit == 'sweep':
        data = _dealias_unwrap_2d(radar, vel_field, nyquist_vel)
    elif unwrap_unit == 'volume':
        data = _dealias_unwrap_3d(radar, vel_field, nyquist_vel)
    else:
        message = ("Unknown `unwrap_unit` parameter, must be one of"
                   "'ray', 'sweep', or 'volume'")
        raise ValueError(message)

    # return field dictionary containing dealiased Doppler velocities
    corr_vel = get_metadata(corr_vel_field)
    corr_vel['data'] = data
    return corr_vel


def _dealias_unwrap_3d(radar, vel_field, nyquist_vel):
    """ Dealias using 3D phase unwrapping (full volume at once). """

    # form cube and scale to phase units
    volume = radar.fields[vel_field]['data']
    scaled_volume = np.pi * volume / nyquist_vel
    scaled_cube = scaled_volume.reshape(radar.nsweeps, -1, radar.ngates)

    # perform unwrapping
    if np.ma.isMaskedArray(scaled_cube):
        mask = np.require(scaled_cube.mask, np.uint8, ['C'])
    else:
        mask = np.zeros_like(scaled_cube, dtype=np.uint8, order='C')
    if mask.ndim == 0:  # fix when the mask is a single True/False
        mask = np.ones_like(scaled_cube, dtype=np.uint8, order='C') * mask

    mask = np.asarray(mask, dtype=np.uint8, order='C')
    cube_no_mask = np.asarray(scaled_cube, dtype=np.float64, order='C')
    cube_scaled_unwrapped = np.empty_like(cube_no_mask, dtype=np.float64,
                                          order='C')
    unwrap_3d(cube_no_mask, mask, cube_scaled_unwrapped,
              [False, True, False])

    # scale back to velocity units, prepare output
    unwrapped_cube = cube_scaled_unwrapped * nyquist_vel / np.pi
    unwrapped_volume = unwrapped_cube.reshape(-1, radar.ngates)
    unwrapped_volume = unwrapped_volume.astype(volume.dtype)

    if np.ma.isMaskedArray(volume):
        return np.ma.array(unwrapped_volume, mask=volume.mask)
    else:
        return unwrapped_volume


def _dealias_unwrap_1d(radar, vel_field, nyquist_vel):
    """ Dealias using 1D phase unwrapping (ray-by-ray) """
    # XXX Does not support masked arrays
    volume = radar.fields[vel_field]['data']
    data = np.empty_like(volume)

    for i, ray in enumerate(volume):
        scaled_ray = ray * np.pi / nyquist_vel

        no_mask = np.array(scaled_ray, dtype=np.float64, order='C')
        scaled_unwrapped = np.empty_like(no_mask, dtype=np.float64, order='C')
        unwrap_1d(no_mask, scaled_unwrapped)
        #scaled_unwrapped = skimage.restoration.unwrap_phase(scaled_ray.data)
        unwrapped = scaled_unwrapped * nyquist_vel / np.pi
        data[i] = unwrapped[:]

    if np.ma.isMaskedArray(radar.fields[vel_field]['data']):
        return np.ma.array(data, mask=radar.fields[vel_field]['data'].mask)
    else:
        return data


def _dealias_unwrap_2d(radar, vel_field, nyquist_vel):
    """ Dealias using 2D phase unwrapping (sweep-by-sweep). """
    data = np.zeros_like(radar.fields[vel_field]['data'])
    for sweep_i in range(radar.nsweeps):

        # extract the sweep
        start = radar.sweep_start_ray_index['data'][sweep_i]
        end = radar.sweep_end_ray_index['data'][sweep_i]
        sweep = radar.fields[vel_field]['data'][start: end + 1]

        # scale to phase units (-pi to pi)
        scaled_sweep = sweep * np.pi / nyquist_vel

        # perform unwrapping
        if np.ma.isMaskedArray(sweep):
            mask = np.require(sweep.mask, np.uint8, ['C'])
        else:
            mask = np.zeros_like(sweep, dtype=np.uint8, order='C')
        if mask.ndim == 0:  # fix when the mask is a single True/False
            mask = np.ones_like(sweep, dtype=np.uint8, order='C') * mask
        mask = np.asarray(mask, dtype=np.uint8, order='C')
        csweep = np.asarray(scaled_sweep, dtype=np.float64, order='C')
        scaled_unwrapped = np.empty_like(csweep, dtype=np.float64, order='C')
        unwrap_2d(csweep, mask, scaled_unwrapped, [True, True])

        # scale back into velocity units
        unwrapped = scaled_unwrapped * nyquist_vel / np.pi
        data[start:end + 1, :] = unwrapped[:]

    if np.ma.isMaskedArray(radar.fields[vel_field]['data']):
        return np.ma.array(data, mask=radar.fields[vel_field]['data'].mask)
    else:
        return data


def _is_radar_cubic(radar):
    """ Test if a radar is cubic (sweeps have the same number of rays). """
    rays_per_sweep = (radar.sweep_end_ray_index['data'] -
                      radar.sweep_start_ray_index['data']) + 1
    return np.all(rays_per_sweep == rays_per_sweep[0])


def _is_radar_sweep_aligned(radar, diff=0.1):
    """
    Test that all sweeps in the radar sample nearly the same angles.

    Test that the maximum difference in sweep sampled angles is below
    `diff` degrees. The radar should first be tested to verify that is cubic
    before calling this function using the _is_radar_cubic function.

    """
    if radar.scan_type == 'ppi':
        angles = radar.azimuth['data']
    elif radar.scan_type == 'rhi':
        angles = radar.elevation['data']
    starts = radar.sweep_start_ray_index['data']
    ends = radar.sweep_end_ray_index['data']
    ref_angles = angles[starts[0]:ends[0] + 1]
    for start, end in zip(starts, ends):
        test_angles = angles[start:end+1]
        if np.any(np.abs(test_angles - ref_angles) > diff):
            return False
    return True


def _is_radar_sequential(radar):
    """ Test if all sweeps in radar are sequentially ordered. """
    for i in xrange(radar.nsweeps):
        if not _is_sweep_sequential(radar, i):
            return False
    return True


def _is_sweep_sequential(radar, sweep_number):
    """ Test if a specific sweep is sequentially ordered. """
    start = radar.sweep_start_ray_index['data'][sweep_number]
    end = radar.sweep_end_ray_index['data'][sweep_number]
    if radar.scan_type == 'ppi':
        angles = radar.azimuth['data'][start:end+1]
    elif radar.scan_type == 'rhi':
        angles = radar.elevation['data'][start:end+1]
    else:
        raise ValueError('invalid scan_type: %s' % (radar.scan_type))
    rolled_angles = np.roll(angles, -np.argmin(angles))
    return np.all(np.diff(rolled_angles) >= 0)
