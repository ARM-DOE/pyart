"""
pyart.correct.dealias
=====================

Front end to the University of Washington 4DD code for Doppler dealiasing.

.. autosummary::
    :toctree: generated/

    dealias_fourdd
    find_time_in_interp_sonde
    _create_rsl_volume

"""
# Nothing from this module is imported to pyart.correct if RSL is not
# installed.

import numpy as np

from ..config import get_field_name, get_fillvalue, get_metadata
from ..io import _rsl_interface
from . import _fourdd_interface
from ..util import datetime_utils


def dealias_fourdd(radar, last_radar=None, sounding_heights=None,
                   sounding_wind_speeds=None, sounding_wind_direction=None,
                   prep=1, filt=1, rsl_badval=131072.0, keep_original=False,
                   extra_prep=False, ncp_min=0.3, rhv_min=0.7,
                   refl_field=None, vel_field=None, corr_vel_field=None,
                   last_vel_field=None, ncp_field=None, rhv_field=None,
                   debug=False, max_shear=0.05, sign=1, **kwargs):
    """
    Dealias Doppler velocities using the 4DD algorithm.

    Dealias the Doppler velocities field using the University of Washington
    4DD algorithm utilizing information from a previous volume scan and/or
    sounding data. Either last_radar or sounding_heights,
    sounding_wind_speeds and sounding_wind_direction must be provided.
    For best results provide both a previous volume scan and sounding data.
    Radar and last_radar must contain the same number of rays per sweep.

    Additional arguments are passed to
    :py:func:`_fourdd_interface.fourdd_dealias`.
    These can be used to fine tune the behavior of the FourDD algorithm.
    See the documentation of Other Parameters for details.

    Parameters
    ----------
    radar : Radar
        Radar object to use for dealiasing.  Must have a Nyquist defined in
        the instrument_parameters attribute and have a
        reflectivity_horizontal and mean_doppler_velocity fields.
    last_radar : Radar, optional
        The previous radar volume, which has been successfully
        dealiased. Using a previous volume as an initial condition can
        greatly improve the dealiasing, and represents the final dimension
        in the 4DD algorithm.
    sounding_heights : ndarray, optional
        Sounding heights is meters above mean sea level.  If altitude
        attribute of the radar object if reference against something other
        than mean sea level then this parameter should also be referenced in
        that manner.
    sounding_wind_speeds : ndarray, optional
        Sounding wind speeds in m/s.
    sounding_wind_direction : ndarray, optional
        Sounding wind directions in degrees.

    Other Parameters
    ----------------
    prep : int
        Flag controlling thresholding, 1 = yes, 0 = no.
    filt : int
        Flag controlling Bergen and Albers filter, 1 = yes, 0 = no.
    rsl_badval : float
        Value which represents a bad value in RSL.
    keep_original : bool
        True to keep original doppler velocity values when the dealiasing
        procedure fails, otherwises these gates will be masked.  NaN values
        are still masked.
    extra_prep : bool
        True to use extra volume preperation, which will use the normalized
        coherent power and cross correlation ratio fields to further remove
        bad gates. Set this to False (default) if those fields are not
        available.
    ncp_min : float
        Minimum normalized coherent power allowed for any gate. Only
        applicable when extra_prep is True.
    rhv_min : float
        Minimum cross correlation ratio allowed for any gate. Only
        applicable when extra_prep is True.
    refl_field : str
        Field in radar to use as the reflectivity during dealiasing.
        None will use the default field name from the Py-ART configuration
        file.
    vel_field : str
        Field in radar to use as the Doppler velocities during dealiasing.
        None will use the default field name from the Py-ART configuration
        file.
    corr_vel_field : str
        Name to use for the dealiased Doppler velocity field metadata.  None
        will use the default field name from the Py-ART configuration file.
    last_vel_field : str
        Name to use for the dealiased Doppler velocity field metadata in
        last_radar.  None will use the corr_vel_field name.
    ncp_field, rhv_field : str
        Fields in radar to use as the normalized coherent power and cross
        correlation ratio, respectively. None will use the default field name
        from the Py-ART configuration file. Only applicable when extra_prep
        is True.
    maxshear : float
        Maximum vertical shear which will be incorperated into the created
        volume from the sounding data.  Parameter not used when no
        sounding data is provided.
    sign : int
        Sign convention which the radial velocities in the volume created
        from the sounding data will will. This should match the convention
        used in the radar data. A value of 1 represents when positive values
        velocities are towards the radar, -1 represents when negative
        velocities are towards the radar.
    debug : bool
        Set True to return RSL Volume objects for debugging:
        usuccess, DZvolume, radialVelVolume, unfoldedVolume, sondVolume

    Returns
    -------
    vr_corr : dict
        Field dictionary containing dealiased Doppler velocities.  Dealiased
        array is stored under the 'data' key.

    Notes
    -----
    Due to limitations in the C code do not call with sounding arrays over
    999 elements long.

    References
    ----------
    C. N. James and R. A Houze Jr, A Real-Time Four-Dimensional Doppler
    Dealising Scheme, Journal of Atmospheric and Oceanic Technology, 2001, 18,
    1674.

    """
    # TODO test with RHI radar scan

    # verify that sounding data or last_volume is provided
    sounding_available = ((sounding_heights is not None) and
                          (sounding_wind_speeds is not None) and
                          (sounding_wind_direction is not None))
    if (not sounding_available) and (last_radar is None):
        raise ValueError('sounding data or last_radar must be provided')

    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if vel_field is None:
        vel_field = get_field_name('velocity')
    if corr_vel_field is None:
        corr_vel_field = get_field_name('corrected_velocity')
    if last_vel_field is None:
        last_vel_field = get_field_name('corrected_velocity')
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')
    if rhv_field is None:
        rhv_field = get_field_name('cross_correlation_ratio')

    # get fill value
    fill_value = get_fillvalue()

    # check extra volume preparation inputs
    if extra_prep:
        if ncp_field not in radar.fields and rhv_field not in radar.fields:
            raise KeyError(('Radar does not have necessary fields for'
                            'extra volume preparation'))

    # extra volume preparation
    # this assumes the radar has a normalized coherent power field and a
    # cross correlation ratio field
    if prep and extra_prep:
        not_coherent = np.logical_or(
            radar.fields[ncp_field]['data'] < ncp_min,
            radar.fields[rhv_field]['data'] < rhv_min)
        fdata = np.copy(radar.fields[refl_field]['data']).astype(np.float32)
        fdata[not_coherent] = fill_value
    else:
        fdata = None

    # create RSL volumes containing the reflectivity, doppler velocity, and
    # doppler velocity in the last radar (if provided)
    refl_volume = _create_rsl_volume(radar, refl_field, 0, rsl_badval, fdata)
    vel_volume = _create_rsl_volume(radar, vel_field, 1, rsl_badval)

    if last_radar is not None:
        last_vel_volume = _create_rsl_volume(
            last_radar, last_vel_field, 1, rsl_badval)
    else:
        last_vel_volume = None

    # create an RslVolume containing the sounding data if it available
    if sounding_available:
        # convert the sounding data to 1D float32 arrays
        hc = np.ascontiguousarray(sounding_heights, dtype=np.float32)
        sc = np.ascontiguousarray(sounding_wind_speeds, dtype=np.float32)
        dc = np.ascontiguousarray(sounding_wind_direction, dtype=np.float32)

        success, sound_volume = _fourdd_interface.create_soundvolume(
            refl_volume, hc, sc, dc, sign, max_shear)
        if success == 0:
            raise ValueError('Error when loading sounding data')
    else:
        sound_volume = None

    # perform dealiasing
    if debug:
        return _fourdd_interface.fourdd_dealias(
            vel_volume, last_vel_volume, sound_volume, refl_volume,
            prep, filt, debug=True, **kwargs)

    flag, data = _fourdd_interface.fourdd_dealias(
        vel_volume, last_vel_volume, sound_volume, refl_volume, prep, filt,
        debug=False, **kwargs)

    # prepare data for output, set bad values and mask data
    is_bad_data = np.logical_or(np.isnan(data), data == rsl_badval)
    if keep_original:
        vel_array = radar.fields[vel_field]['data']
        data = np.where(is_bad_data, vel_array, data)
    else:
        data[is_bad_data] = fill_value
    data = np.ma.masked_equal(data, fill_value)

    # return field dictionary containing dealiased Doppler velocities
    vr_corr = get_metadata(corr_vel_field)
    vr_corr['data'] = data
    vr_corr['_FillValue'] = data.fill_value
    return vr_corr


def _create_rsl_volume(radar, field_name, vol_num, rsl_badval, fdata=None):
    """
    Create a RSLVolume containing data from a field in radar.
    """
    fill_value = get_fillvalue()
    if fdata is None:
        fdata = np.copy(radar.fields[field_name]['data']).astype(np.float32)
    fdata = np.ma.filled(fdata, fill_value)
    is_bad = np.logical_or(fdata == fill_value, np.isnan(fdata))
    fdata[is_bad] = rsl_badval
    rays_per_sweep = (radar.sweep_end_ray_index['data'] -
                      radar.sweep_start_ray_index['data'] + 1)
    rays_per_sweep = rays_per_sweep.astype(np.int32)
    rsl_volume = _rsl_interface.create_volume(fdata, rays_per_sweep, vol_num)
    _rsl_interface._label_volume(rsl_volume, radar)
    return rsl_volume


def find_time_in_interp_sonde(interp_sonde, target, debug=False):
    """
    Find the wind parameter for a given time in a ARM interpsonde file.

    Parameters
    ----------
    interp_sonde : netCDF4.Dataset
        netCDF4 object pointing to a ARM interpsonde file.
    target : datetime
        Target datetime, the closest time in the interpsonde file will be
        used.

    Other Parameters
    ----------------
    debug : bool
        Print debugging information.

    Returns
    -------
    height : np.ndarray
        Heights above the ground for the time closest to target.
    speed : np.ndarray
        Wind speeds at given height for the time closest to taget.
    direction : np.ndarray
        Wind direction at given height for the time closest to target.

    """

    sonde_datetimes = datetime_utils.datetimes_from_dataset(interp_sonde)

    # get closest time index to target
    idx = np.abs(sonde_datetimes - target).argmin()

    if debug:
        print('Target time is %s' % (target))
        print('Interpolated sounding time is %s' % (sonde_datetimes[idx]))

    return (interp_sonde.variables['height'][:],
            interp_sonde.variables['wspd'][idx, :],
            interp_sonde.variables['wdir'][idx, :])
