"""
pyart.graph.coord_transform
===========================

Coordinate transform routines.

.. autosummary::
    :toctree: generated/

    radar_coords_to_cart_track_relative
    radar_coords_to_cart_earth_relative
    radar_coords_to_cart_aircraft_relative

"""

import numpy as np


def radar_coords_to_cart_track_relative(rng, rot, roll, drift, tilt, pitch):
    """
    Calculate track-relative Cartesian coordinates from radar coordinates.

    Parameters
    ----------
    rng : array
        Distances to the center of the radar gates (bins) in kilometers.
    rot : array
        Rotation angle of the radar in degrees.
    roll : array
        Roll angle of the radar in degrees.
    drift : array
        Drift angle of the radar in degrees.
    tilt : array
        Tilt angle of the radar in degrees.
    pitch : array
        Pitch angle of the radar in degrees.

    Returns
    -------
    x, y, z : array
        Cartesian coordinates in meters from the radar.

    Notes
    -----
    Project native (polar) coordinate radar sweep data onto
    track-relative Cartesian coordinate grid.

    References
    ----------
    .. [1] Lee et al. (1994) Journal of Atmospheric and Oceanic Technology.

    """
    rot = np.radians(rot)               # rotation angle in radians.
    roll = np.radians(roll)             # roll angle in radians.
    drift = np.radians(drift)           # drift angle in radians.
    tilt = np.radians(tilt)             # tilt angle in radians.
    pitch = np.radians(pitch)           # pitch angle in radians.
    r = rng * 1000.0                    # distances to gates in meters.

    x = r * (np.cos(rot + roll) * np.sin(drift) * np.cos(tilt) *
             np.sin(pitch) + np.cos(drift) * np.sin(rot + roll) *
             np.cos(tilt) - np.sin(drift) * np.cos(pitch) * np.sin(tilt))
    y = r * (-1. * np.cos(rot + roll) * np.cos(drift) * np.cos(tilt) *
             np.sin(pitch) + np.sin(drift) * np.sin(rot + roll) *
             np.cos(tilt) + np.cos(drift) * np.cos(pitch) * np.sin(tilt))
    z = (r * np.cos(pitch) * np.cos(tilt) * np.cos(rot + roll) +
         np.sin(pitch) * np.sin(tilt))
    return x, y, z


def radar_coords_to_cart_earth_relative(rng, rot, roll, heading, tilt, pitch):
    """
    Calculate earth-relative Cartesian coordinates from radar coordinates

    Parameters
    ----------
    rng : array
        Distances to the center of the radar gates (bins) in kilometers.
    rot : array
        Rotation angle of the radar in degrees.
    roll : array
        Roll angle of the radar in degrees.
    heading : array
        Heading (compass) angle of the radar in degrees clockwise from north.
    tilt : array
        Tilt angle of the radar in degrees.
    pitch : array
        Pitch angle of the radar in degrees.

    Returns
    -------
    x, y, z : array
        Cartesian coordinates in meters from the radar.

    Notes
    -----
    Project native (polar) coordinate radar sweep data onto
    earth-relative Cartesian coordinate grid.

    References
    ----------
    .. [1] Lee et al. (1994) Journal of Atmospheric and Oceanic Technology.

    """
    rot = np.radians(rot)               # rotation angle in radians.
    roll = np.radians(roll)             # roll angle in radians.
    heading = np.radians(heading)       # drift angle in radians.
    tilt = np.radians(tilt)             # tilt angle in radians.
    pitch = np.radians(pitch)           # pitch angle in radians.
    r = rng * 1000.0                    # distances to gates in meters.

    x = r * (-1. * np.cos(rot + roll) * np.sin(heading) * np.cos(tilt) *
             np.sin(pitch) + np.cos(heading) * np.sin(rot + roll) *
             np.cos(tilt) + np.sin(heading) * np.cos(pitch) * np.sin(tilt))
    y = r * (-1. * np.cos(rot + roll) * np.cos(heading) * np.cos(tilt) *
             np.sin(pitch) - np.sin(heading) * np.sin(rot + roll) *
             np.cos(tilt) + np.cos(heading) * np.cos(pitch) * np.sin(tilt))
    z = (r * np.cos(pitch) * np.cos(tilt) * np.cos(rot + roll) +
         np.sin(pitch) * np.sin(tilt))
    return x, y, z


def radar_coords_to_cart_aircraft_relative(rng, rot, tilt):
    """
    Calculate aircraft-relative Cartesian coordinates from radar coordinates.

    Parameters
    ----------
    rng : array
        Distances to the center of the radar gates (bins) in kilometers.
    rot : array
        Rotation angle of the radar in degrees.
    tilt : array
        Tilt angle of the radar in degrees.

    Returns
    -------
    X, Y, Z : array
        Cartesian coordinates in meters from the radar.

    Notes
    -----
    Project native (polar) coordinate radar sweep data onto
    earth-relative Cartesian coordinate grid.

    References
    ----------
    .. [1] Lee et al. (1994) Journal of Atmospheric and Oceanic Technology.

    """
    rot = np.radians(rot)               # rotation angle in radians.
    tilt = np.radians(tilt)             # tilt angle in radians.
    r = rng * 1000.0                    # distances to gates in meters.
    x = r * np.cos(tilt) * np.sin(rot)
    y = r * np.sin(tilt)
    z = r * np.cos(rot) * np.cos(tilt)
    return x, y, z
