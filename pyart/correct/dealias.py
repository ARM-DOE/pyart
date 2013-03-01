"""
Front end to the University of Washington 4DD code for Doppler dealiasing.
"""

import numpy as np

# XXX fix these up when io is reorganized
from ..io import rsl, fourdd, rsl_utils
from ..io.common import get_metadata

from ..util import datetime_utils


def dealias(radar, sounding_heights, sounding_wind_speeds,
            sounding_wind_direction, datetime_sounding, rsl_radar=None,
            prep=1, filt=1, rsl_badval=131072, fill_value=-9999.0,
            refl='reflectivity_horizontal', vel='mean_doppler_velocity'):
    """
    Dealias the Doppler velocities field using the University of Washington
    4DD algorithm utilizing information from sounding data.

    Parameters
    ----------
    radar : Radar
        Radar object to use for dealiasing.  Must have a Nyquist defined in the
        inst_params attribute and have a reflectivity_horizontal and
        mean_doppler_velocity fields.  Data is not extracted from this object
        for dealiasing if `rsl_radar` is defined.
    sounding_heights : array
        Sounding heights is meters.
    sounding_wind_speeds : array
        Sounding wind speeds in m/s.
    sounding_wind_direction : array
        Sounding wind directions in degrees.
    datetime_sounding : datetime
        Datetime representing the mean time of the sounding profile.
    rsl_radar : RSL Radar, optional
        RSL Radar to use for dealiasing, None to extract this data from radar.

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
    # XXX add check for arrays over 900 elements (???)

    # create a RSL Radar if one is not provided
    if rsl_radar is None:
        rsl_radar = rsl_utils.radar_to_rsl(radar, {refl: 'DZ', vel: 'VR'})

    # perform dealiasing
    juldate = int(datetime_sounding.strftime('%y%j%H%M'))   # YYDDDHHMM
    new_volume, _ = fourdd.dealias_radar_array(
        rsl_radar, juldate, sounding_heights, sounding_wind_speeds,
        sounding_wind_direction, None, prep=prep, filt=filt)

    # extract data, mask and reshape
    dealiased_data = rsl.create_cube_array_lim(
        new_volume[0], new_volume.contents.h.nsweeps,
        new_volume.contents.sweeps[0].h.nrays)
    dealiased_data[np.where(np.isnan(dealiased_data))] = fill_value
    dealiased_data[np.where(dealiased_data == rsl_badval)] = fill_value
    dealiased_data = np.ma.masked_equal(dealiased_data, -9999.0)
    dealiased_data.shape = (-1, dealiased_data.shape[2])

    # build the field dictionary
    dealiased_fielddict = {'data': dealiased_data}
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
