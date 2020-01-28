"""
Utilities for ARM sonde NetCDF files.

"""

import netCDF4
import numpy as np

from ..core import HorizontalWindProfile


def read_arm_sonde(filename):
    """
    Read a ARM sonde file returning a wind profile.

    Parameters
    ----------
    filename : str
        Name of ARM sonde NetCDF file to read data from.

    Return
    ------
    launch_datetime : datetime
        Date and time corresponding to radiosonde launch time, i.e., first
        recorded time.
    profile : HorizontalWindProfile
        Profile of the horizontal winds

    """
    dset = netCDF4.Dataset(filename, 'r')

    launch_datetime = netCDF4.num2date(
        dset.variables['time'][0], dset.variables['time'].units)

    height = dset.variables['alt'][:]  # meters above mean sea level
    speed = dset.variables['wspd'][:]  # m s-1
    direction = dset.variables['deg'][:]  # degrees clockwise from north
    lat = dset.variables['lat'][:]  # degrees north
    lon = dset.variables['lon'][:]  # degrees east
    profile = HorizontalWindProfile(
        height, speed, direction, latitude=lat, longitude=lon)

    dset.close()

    return launch_datetime, profile


def read_arm_sonde_vap(filename, radar=None, target_datetime=None):
    """
    Read a ARM interpolated or merged sonde returning a wind profile.

    Parameters
    ----------
    filename : str
        Name of ARM interpolate or merged sonde NetCDF file to read data from.
    radar : Radar, optional
        If provided the profile returned is that which is closest in time to
        the first ray collected in this radar. Either radar or target_datetime
        must be provided.
    target_datetime : datetime, optional
        If specified the profile returned is that which is closest in time to
        this datetime.

    Return
    ------
    profile_datetime : datetime
        Date and time of the profile.
    profile : HorizontalWindProfile
        Profile of the horizontal winds.

    """

    # check parameters
    if radar is None and target_datetime is None:
        raise ValueError('Either radar or target_datetime must be specified.')

    if radar is not None and target_datetime is not None:
        raise ValueError(
            'Either radar or target_datetime must be specified, not both.')

    if radar is not None:
        time_0 = radar.time['data'][0]
        time_units = radar.time['units']
        target_datetime = netCDF4.num2date(time_0, time_units)

    dset = netCDF4.Dataset(filename, 'r')

    # find index of time closest to target datetime
    sonde_datetimes = netCDF4.num2date(
        dset.variables['time'][:], dset.variables['time'].units)
    idx = np.abs(sonde_datetimes - target_datetime).argmin()
    profile_datetime = sonde_datetimes[idx]

    # extract wind profile
    height = dset.variables['height'][:] * 1000   # km -> m
    speed = dset.variables['wspd'][idx, :]
    direction = dset.variables['wdir'][idx, :]
    profile = HorizontalWindProfile(height, speed, direction)

    dset.close()

    return profile_datetime, profile
