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

from ..io import _rsl_interface
from . import _fourdd_interface
from ..io.common import get_metadata
from ..util import datetime_utils


def dealias_fourdd(radar, sounding_heights, sounding_wind_speeds,
                   sounding_wind_direction, datetime_sounding,
                   prep=1, filt=1, rsl_badval=131072, fill_value=-9999.0,
                   refl='reflectivity_horizontal',
                   vel='mean_doppler_velocity',
                   debug=False):
    """
    Dealias the Doppler velocities field using the University of Washington
    4DD algorithm utilizing information from sounding data.

    Parameters
    ----------
    radar : Radar
        Radar object to use for dealiasing.  Must have a Nyquist defined in
        the instrument_parameters attribute and have a
        reflectivity_horizontal and mean_doppler_velocity fields.
    sounding_heights : array
        Sounding heights is meters above mean sea level.  If altitude
        attribute of the radar object if reference against something other
        than mean sea level then this parameter should also be referenced in
        that manner.
    sounding_wind_speeds : array
        Sounding wind speeds in m/s.
    sounding_wind_direction : array
        Sounding wind directions in degrees.
    datetime_sounding : datetime
        Datetime representing the mean time of the sounding profile.

    Other Parameters
    ----------------
    prep : int
        Flag controlling thresholding, 1 = yes, 0 = no.
    filt : int
        Flag controlling Bergen and Albers filter, 1 = yes, 0 = no.
    rsl_badval : int
        Value which represented a bad, masked, points in RSL.
    fill_value : float
        Value to fill array with at masked points.
    refl : str
        Field in radar to use as the doppler velocities during dealiasing.
    vel : str
        Field in radar to use as the reflectivity during dealiasing
    debug : bool
        Set True to return RSL Volume objects for debugging:
        usuccess, DZvolume, radialVelVolume, unfoldedVolume, sondVolume

    Returns
    -------
    dealiased_fielddict : dict
        Field dictionary containing dealiased doppler velocities.  Dealiased
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
    # TODO use radar with recent correct vel instead of sond data
    # TODO option not to use refl.
    # TODO test with RHI radar

    # convert the sounding data to 1D float32 arrays
    hc = np.ascontiguousarray(sounding_heights, dtype='float32')
    sc = np.ascontiguousarray(sounding_wind_speeds, dtype='float32')
    dc = np.ascontiguousarray(sounding_wind_direction, dtype='float32')

    # interger representation of the time
    vad_time = int(datetime_sounding.strftime('%y%j%H%M'))   # YYDDDHHMM

    # extract the reflectivity and create a RslVolume containing it
    nshape = (radar.nsweeps, -1, radar.ngates)
    refl_array = radar.fields[refl]['data'].reshape(nshape).astype('float32')
    refl_array[np.where(refl_array == fill_value)] = rsl_badval
    refl_array[np.where(np.isnan(refl_array))] = rsl_badval
    refl_volume = _rsl_interface.create_volume(refl_array, 0)
    _rsl_interface._label_volume(refl_volume, radar)

    # extract the doppler velocity and create a RslVolume containing it
    vel_array = radar.fields[vel]['data'].reshape(nshape).astype('float32')
    vel_array[np.where(vel_array == fill_value)] = rsl_badval
    vel_array[np.where(np.isnan(vel_array))] = rsl_badval
    vel_volume = _rsl_interface.create_volume(vel_array, 1)
    _rsl_interface._label_volume(vel_volume, radar)

    # perform dealiasing
    if debug:
        return _fourdd_interface.fourdd_dealias(
            refl_volume, vel_volume, hc, sc, dc, vad_time, prep, filt, True)
    flag, data = _fourdd_interface.fourdd_dealias(
        refl_volume, vel_volume, hc, sc, dc, vad_time, prep, filt, False)

    # reshape and mask data
    data.shape = (-1, radar.ngates)
    data[np.where(np.isnan(data))] = fill_value
    data[np.where(data == rsl_badval)] = fill_value
    data = np.ma.masked_equal(data, fill_value)

    # create and return field dictionary containing dealiased data
    dealiased_fielddict = {'data': data}
    meta = get_metadata('VEL_COR')
    dealiased_fielddict.update(meta)
    return dealiased_fielddict


def find_time_in_interp_sonde(interp_sonde, target):
    """
    Find the wind parameter for a given time in a ARM interpsonde file.

    Parameters
    ----------
    interp_sonde : netCDF4.Dataset
        netCDF4 object pointing to a ARM interpsonde file.
    target : datetime
        Target datetime, the closest time in the interpsonde file will be
        used.

    Returns
    -------
    height : array
        Heights above the ground for the time closest to target.
    speed : array
        Wind speeds at given height for the time closest to taget.
    direction : array
        Wind direction at given height for the time closest to target.

    """
    sonde_datetimes = datetime_utils.datetimes_from_dataset(interp_sonde)

    # as of numpy v1.7 the next two lines can be replaced with
    # idx = np.abs(sonde_datetimes - radar_datetime).argmin()
    selected = sorted(sonde_datetimes,
                      key=lambda x: abs(x - target))[0]
    idx = list(sonde_datetimes).index(selected)

    return (interp_sonde.variables['height'][:],
            interp_sonde.variables['wspd'][idx, :],
            interp_sonde.variables['wdir'][idx, :])
