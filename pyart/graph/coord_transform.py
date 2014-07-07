"""
pyart.graph.coord_transform
==================

Coordinate transform routines.

.. autosummary::
    :toctree: generated/

    radar_coords_to_cart_track_relative

"""

import numpy as np

def radar_coords_to_cart_track_relative(rng, rot, roll, drift, tilt, pitch, debug=False):
    """
    Calculate track-relative Cartesian coordinates from radar coordinates

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
    X, Y, Z : array
        Cartesian coordinates in meters from the radar.

    Notes
    -----
    Project native (polar) coordinate radar sweep data onto 
    track-relative Cartesian coordinate grid.

    .. math::

        z = 

        s = R * arcsin(\\frac{r*cos(\\theta_e)}{R+z})

        x = 
        y = 

    Where r is the distance from the radar to the center of the gate,
    :math:`\\theta_a` is the azimuth angle, :math:`\\theta_e` is the
    elevation angle, s is the arc length, and R is the effective radius
    of the earth, taken to be 4/3 the mean radius of earth (6371 km).

    References
    ----------
    .. [1] Lee et al. (1994) Journal of Atmospheric and Oceanic Technology.

    """
    Rot = np.radians(rot)               # rotation angle in radians.
    Roll = np.radians(roll)             # roll angle in radians.
    Drift = np.radians(drift)           # drift angle in radians.
    Tilt = np.radians(tilt)             # tilt angle in radians.
    Pitch = np.radians(pitch)           # pitch angle in radians.
    r = rng * 1000.0                    # distances to gates in meters.
    
    X = r * (np.cos(Rot + Roll) * np.sin(Drift) * np.cos(Tilt) * np.sin(Pitch) +
        np.cos(Drift) * np.sin(Rot + Roll) * np.cos(Tilt) - 
        np.sin(Drift) * np.cos(Pitch) * np.sin(Tilt))
    Y = r * (-1.* np.cos(Rot + Roll) * np.cos(Drift) * np.cos(Tilt) * np.sin(Pitch) +
        np.sin(Drift) * np.sin(Rot + Roll) * np.cos(Tilt) + 
        np.cos(Drift) * np.cos(Pitch) * np.sin(Tilt))
    Z = r * np.cos(Pitch) * np.cos(Tilt) * np.cos(Rot + Roll) + np.sin(Pitch) * np.sin(Tilt)
    return X, Y, Z
    
def radar_coords_to_cart_earth_relative(rng, rot, roll, heading, tilt, pitch, debug=False):
    """,tilt,heading,pitch
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
    X, Y, Z : array
        Cartesian coordinates in meters from the radar.

    Notes
    -----
    Project native (polar) coordinate radar sweep data onto 
    earth-relative Cartesian coordinate grid.

    .. math::

        z = 

        s = R * arcsin(\\frac{r*cos(\\theta_e)}{R+z})

        x = 
        y = 

    Where r is the distance from the radar to the center of the gate,
    :math:`\\theta_a` is the azimuth angle, :math:`\\theta_e` is the
    elevation angle, s is the arc length, and R is the effective radius
    of the earth, taken to be 4/3 the mean radius of earth (6371 km).

    References
    ----------
    .. [1] Lee et al. (1994) Journal of Atmospheric and Oceanic Technology.

    """
    Rot = np.radians(rot)               # rotation angle in radians.
    Roll = np.radians(roll)             # roll angle in radians.
    Heading = np.radians(heading)           # drift angle in radians.
    Tilt = np.radians(tilt)             # tilt angle in radians.
    Pitch = np.radians(pitch)           # pitch angle in radians.
    r = rng * 1000.0                    # distances to gates in meters.
    
    X = r * (-1.* np.cos(Rot + Roll) * np.sin(Heading) * np.cos(Tilt) * np.sin(Pitch) +
        np.cos(Heading) * np.sin(Rot + Roll) * np.cos(Tilt) + 
        np.sin(Heading) * np.cos(Pitch) * np.sin(Tilt))
    Y = r * (-1.* np.cos(Rot + Roll) * np.cos(Heading) * np.cos(Tilt) * np.sin(Pitch) -
        np.sin(Heading) * np.sin(Rot + Roll) * np.cos(Tilt) + 
        np.cos(Heading) * np.cos(Pitch) * np.sin(Tilt))
    Z = r * np.cos(Pitch) * np.cos(Tilt) * np.cos(Rot + Roll) + np.sin(Pitch) * np.sin(Tilt)
    return X, Y, Z
    
def radar_coords_to_cart_aircraft_relative(rng, rot, tilt, debug=False):
    """
    Calculate aircraft-relative Cartesian coordinates from radar coordinates

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

    .. math::

        z = 

        s = R * arcsin(\\frac{r*cos(\\theta_e)}{R+z})

        x = 
        y = 

    Where r is the distance from the radar to the center of the gate,
    :math:`\\theta_a` is the azimuth angle, :math:`\\theta_e` is the
    elevation angle, s is the arc length, and R is the effective radius
    of the earth, taken to be 4/3 the mean radius of earth (6371 km).

    References
    ----------
    .. [1] Lee et al. (1994) Journal of Atmospheric and Oceanic Technology.

    """
    Rot = np.radians(rot)               # rotation angle in radians.
    Tilt = np.radians(tilt)             # tilt angle in radians.
    r = rng * 1000.0                    # distances to gates in meters.
    
    X = r * np.cos(Tilt) * np.sin(Rot)
    Y = r * np.sin(Tilt)
    Z = r * np.cos(Rot) * np.cos(Tilt)
    return X, Y, Z
    
def latlon2xy(latitude,longitude,altitude,): ####NOT IN WORKING ORDER######

    """/* These calculations are from the book
    * "Aerospace Coordinate Systems and Transformations"
    * by G. Minkler/J. Minkler
    * these are the ECEF/ENU point transformations
    """

    earth_rad = 6370000.
    h = earth_rad + altitude[0]
    delta_o = np.radians(latitude[0])
    lambda_o = np.radians(longitude[0])

    sinLambda = np.sin(lambda_o)
    cosLambda = np.cos(lambda_o)
    sinDelta = np.sin(delta_o)
    cosDelta = np.cos(delta_o)
    
    R_p = earth_rad + altitude
    delta_p = np.radians(latitude)
    lambda_p = np.radians(longitude)
    R_p_pr = R_p * np.cos(delta_p)

    xe = R_p * np.sin(delta_p)
    ye = -R_p_pr * np.sin(lambda_p)
    ze = R_p_pr * np. cos(lambda_p)

# transform to ENU coordinates */

    a = -h * sinDelta + xe
    b =  h * cosDelta * sinLambda + ye
    c = -h * cosDelta * cosLambda + ze

    X = -cosLambda * b - (sinLambda * c);
    Y = (cosDelta * a)  +  (sinLambda * sinDelta * b) - (cosLambda * sinDelta * c)
    Z = (sinDelta * a)  - (sinLambda * cosDelta * b) + (cosLambda * cosDelta * c)
    return X, Y, Z