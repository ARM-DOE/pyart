
import numpy as np

PI = 3.141592653589793


def dms_to_d(dms):
    return dms[0] + (dms[1] + dms[2] / 60.0) / 60.0


def dt_to_dict(dt, **kwargs):
    """
    Returns a dictionary from a datetime object perfect for use in a
    formatted string usage:
        dict = (datetime.datetime object, pref = prefix
                string for the dictionary keys)

    """
    pref = kwargs.get('pref', '')
    return dict([(pref+key, getattr(dt, key)) for key in
                ['year', 'month', 'day', 'hour', 'minute', 'second']])


def corner_to_point(corner, point):
    """
    Distance from a corner to a point in lat-lons given a spherical earth
    usage: x_distance, y_distance=corner_to_point([lat1, lon1], [lat2, lon2])
    """
    Re = 6371.0 * 1000.0
    Rc = ax_radius(point[0], units='degrees')
    #print Rc/Re
    y = ((point[0] - corner[0]) / 360.0) * PI * 2.0 * Re
    x = ((point[1] - corner[1]) / 360.0) * PI * 2.0 * Rc
    return x, y


def ax_radius(lat, units='radians'):
    """Determine the radius of a circle of constant longitude at a certain
    Latitude usage radius=(latitude, units='degrees' or 'radians')
    """
    Re = 6371.0 * 1000.0
    if units == 'degrees':
        const = PI / 180.0
    else:
        const = 1.0
    R = Re * np.sin(PI / 2.0 - abs(lat * const))
    return R


def radar_coords_to_cart(rng, az, ele, debug=False):
    """
    Asumes standard atmosphere, ie R=4Re/3
    rng in km
    """
    Re = 6371.0 * 1000.0
    p_r = 4.0 * Re / 3.0
    rm = rng * 1000.0
    z = (rm ** 2 + p_r ** 2 + 2.0 * rm * p_r *
         np.sin(ele * PI / 180.0)) ** 0.5 - p_r
    #arc length
    s = p_r * np.arcsin(rm * np.cos(ele * PI / 180.) / (p_r + z))
    if debug:
        print "Z=", z, "s=", s
    y = s * np.cos(az * PI / 180.0)
    x = s * np.sin(az * PI / 180.0)
    return x, y, z
