"""
pyart.graph.common
==================

Common graphing routines.

.. autosummary::
    :toctree: generated/

    dms_to_d
    radar_coords_to_cart
    corner_to_point
    ax_radius

"""

import numpy as np

PI = 3.141592653589793

from ..io.common import dms_to_d, radar_coords_to_cart


def corner_to_point(corner, point):
    """
    Return the x, y distances in meters from a corner to a point.

    Assumes a spherical earth model.

    Parameters
    ----------
    corner : (float, float)
        Latitude and longitude in degrees of the corner.
    point : (float, float)
        Latitude and longtide in degrees of the point.

    Returns
    -------
    x, y : floats
        Distances from the corner to the point in meters.

    """
    Re = 6371.0 * 1000.0
    Rc = ax_radius(point[0], units='degrees')
    #print Rc/Re
    y = ((point[0] - corner[0]) / 360.0) * PI * 2.0 * Re
    x = ((point[1] - corner[1]) / 360.0) * PI * 2.0 * Rc
    return x, y


def ax_radius(lat, units='radians'):
    """
    Return the radius of a constant longitude circle at a given latitude.

    Parameters
    ----------
    lat : float
        Latitude to calculate constant longitude circle radius.
    units : 'radians' or 'degrees'
        Units of lat, either 'radians' or 'degrees'.

    Returns
    -------
    R : float
        Radius in meters of a constant longitude circle at latitude.

    """
    Re = 6371.0 * 1000.0
    if units == 'degrees':
        const = PI / 180.0
    else:
        const = 1.0
    R = Re * np.sin(PI / 2.0 - abs(lat * const))
    return R
