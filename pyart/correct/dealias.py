"""
pyart.correct.dealias
=====================

Front end to the University of Washington 4DD code for Doppler dealiasing.

.. autosummary::
    :toctree: generated/

    dealias_fourdd
    find_time_in_interp_sonde


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
                   refl_field=None, vel_field=None, corr_vel_field=None,
                   last_vel_field=None, debug=False):
    """
    Dealias Doppler velocities using the 4DD algorithm.

    Dealias the Doppler velocities field using the University of Washington
    4DD algorithm utilizing information from a previous volume scan and/or
    sounding data. Either last_radar or sounding_heights,
    sounding_wind_speeds and sounding_wind_direction must be provided.
    For best results provide both a previous volume scan and sounding data.
    Radar and last_radar must contain the same number of rays per sweep.

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
    Due to limitations in the C code do not call with numpy arrays over 900
    elements long

    References
    ----------
    C. N. James and R. A Houze Jr, A Real-Time Four-Dimensional Doppler
    Dealising Scheme, Journal of Atmospheric and Oceanic Technology, 2001, 18,
    1674.

    """
    # TODO test with RHI radar scan

    # verify that sounding data or last_volume is provided
    sounding_available = (
        None not in [sounding_heights, sounding_wind_speeds,
                     sounding_wind_direction])
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
        last_vel_field = str(corr_vel_field)

    # get fill value
    fill_value = get_fillvalue()

    # transform the dimension of the field arrays from the CF/Radial standard
    # of (ns*nr, ng) to the RSL standard of (ns, nr, ng) where
    #
    # ns = number of sweeps
    # nr = number of rays per sweep
    # ng = number of range gates
    nshape = (radar.nsweeps, -1, radar.ngates)

    # extract the reflectivity and create a RslVolume containing it
    refl_array = np.copy(radar.fields[refl_field]['data'])
    refl_array = refl_array.reshape(nshape).astype(np.float32)
    refl_array = np.ma.filled(refl_array, fill_value)
    is_bad_refl = np.logical_or(refl_array == fill_value, np.isnan(refl_array))
    refl_array[is_bad_refl] = rsl_badval
    refl_volume = _rsl_interface.create_volume(refl_array, 0)
    _rsl_interface._label_volume(refl_volume, radar)

    # extract the Doppler velocity and create a RslVolume containing it
    vel_array = np.copy(radar.fields[vel_field]['data'])
    vel_array = vel_array.reshape(nshape).astype(np.float32)
    vel_array = np.ma.filled(vel_array, fill_value)
    is_bad_vel = np.logical_or(vel_array == fill_value, np.isnan(vel_array))
    vel_array[is_bad_vel] = rsl_badval
    vel_volume = _rsl_interface.create_volume(vel_array, 1)
    _rsl_interface._label_volume(vel_volume, radar)

    # extract the Doppler velocity from the previous volume scan and create a
    # RslVolume containing it
    if last_radar is not None:
        nshape_last = (last_radar.nsweeps, -1, last_radar.ngates)
        last_vel_array = np.copy(last_radar.fields[last_vel_field]['data'])
        last_vel_array = last_vel_array.reshape(nshape_last).astype(np.float32)
        last_vel_array = np.ma.filled(last_vel_array, fill_value)
        is_bad_last = np.logical_or(last_vel_array == fill_value,
                                    np.isnan(last_vel_array))
        last_vel_array[is_bad_last] = rsl_badval
        last_vel_volume = _rsl_interface.create_volume(last_vel_array, 1)
        _rsl_interface._label_volume(last_vel_volume, last_radar)
    else:
        last_vel_volume = None

    # create an RslVolume containing the sounding data if it available
    if sounding_available:
        # convert the sounding data to 1D float32 arrays
        hc = np.ascontiguousarray(sounding_heights, dtype=np.float32)
        sc = np.ascontiguousarray(sounding_wind_speeds, dtype=np.float32)
        dc = np.ascontiguousarray(sounding_wind_direction, dtype=np.float32)

        success, sound_volume = _fourdd_interface.create_soundvolume(
            refl_volume, hc, sc, dc, 1)
        if success == 0:
            raise ValueError('Error when loading sounding data')
    else:
        sound_volume = None

    # perform dealiasing
    if debug:
        return _fourdd_interface.fourdd_dealias(
            vel_volume, last_vel_volume, sound_volume, refl_volume,
            prep, filt, False)

    flag, data = _fourdd_interface.fourdd_dealias(
        vel_volume, last_vel_volume, sound_volume, refl_volume, prep, filt,
        False)

    # prepare data for output, which includes reshaping, setting bad values,
    # and masking data
    data.shape = (-1, radar.ngates)
    is_bad_data = np.logical_or(np.isnan(data), data == rsl_badval)

    if keep_original:
        vel_array.shape = (-1, radar.ngates)
        vel_array[vel_array == rsl_badval] = fill_value
        data = np.where(is_bad_data, vel_array, data)
    else:
        data[is_bad_data] = fill_value

    data = np.ma.masked_equal(data, fill_value)

    # create and return field dictionary containing dealiased Doppler
    # velocities
    vr_corr = get_metadata(corr_vel_field)
    vr_corr['data'] = data
    vr_corr['_FillValue'] = data.fill_value

    return vr_corr


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

    Optional parameters
    -------------------
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
        print 'Target time is %s' % (target)
        print 'Interpolated sounding time is %s' % (sonde_datetimes[idx])

    return (interp_sonde.variables['height'][:],
            interp_sonde.variables['wspd'][idx, :],
            interp_sonde.variables['wdir'][idx, :])
