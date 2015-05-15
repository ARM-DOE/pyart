"""
pyart.io.common
===============

Input/output routines common to many file formats.

.. autosummary::
    :toctree: generated/

    dms_to_d
    stringarray_to_chararray
    radar_coords_to_cart
    make_time_unit_str
    add_2d_latlon_axis

"""

import numpy as np
import netCDF4


def dms_to_d(dms):
    """ Degrees, minutes, seconds to degrees """
    return dms[0] + (dms[1] + dms[2] / 60.0) / 60.0


def stringarray_to_chararray(arr, numchars=None):
    """
    Convert an string array to a character array with one extra dimension.

    Parameters
    ----------
    arr : array
        Array with numpy dtype 'SN', where N is the number of characters
        in the string.

    numchars : int
        Number of characters used to represent the string.  If numchar > N
        the results will be padded on the right with blanks.  The default,
        None will use N.

    Returns
    -------
    chararr : array
        Array with dtype 'S1' and shape = arr.shape + (numchars, ).

    """
    carr = netCDF4.stringtochar(arr)
    if numchars is None:
        return carr

    arr_numchars = carr.shape[-1]
    if numchars <= arr_numchars:
        raise ValueError('numchars must be >= %i' % (arr_numchars))
    chararr = np.zeros(arr.shape + (numchars, ), dtype='S1')
    chararr[..., :arr_numchars] = carr[:]
    return chararr


# XXX move this to another module
def radar_coords_to_cart(rng, az, ele, debug=False):
    """
    Calculate Cartesian coordinate from radar coordinates

    Parameters
    ----------
    rng : array
        Distances to the center of the radar gates (bins) in kilometers.
    az : array
        Azimuth angle of the radar in degrees.
    ele : array
        Elevation angle of the radar in degrees.

    Returns
    -------
    x, y, z : array
        Cartesian coordinates in meters from the radar.

    Notes
    -----
    The calculation for Cartesian coordinate is adapted from equations
    2.28(b) and 2.28(c) of Doviak and Zrnic [1]_ assuming a
    standard atmosphere (4/3 Earth's radius model).

    .. math::

        z = \\sqrt{r^2+R^2+r*R*sin(\\theta_e)} - R

        s = R * arcsin(\\frac{r*cos(\\theta_e)}{R+z})

        x = s * sin(\\theta_a)

        y = s * cos(\\theta_a)

    Where r is the distance from the radar to the center of the gate,
    :math:`\\theta_a` is the azimuth angle, :math:`\\theta_e` is the
    elevation angle, s is the arc length, and R is the effective radius
    of the earth, taken to be 4/3 the mean radius of earth (6371 km).

    References
    ----------
    .. [1] Doviak and Zrnic, Doppler Radar and Weather Observations, Second
        Edition, 1993, p. 21.

    """
    theta_e = ele * np.pi / 180.0       # elevation angle in radians.
    theta_a = az * np.pi / 180.0        # azimuth angle in radians.
    R = 6371.0 * 1000.0 * 4.0 / 3.0     # effective radius of earth in meters.
    r = rng * 1000.0                    # distances to gates in meters.

    z = (r ** 2 + R ** 2 + 2.0 * r * R * np.sin(theta_e)) ** 0.5 - R
    s = R * np.arcsin(r * np.cos(theta_e) / (R + z))  # arc length in m.
    x = s * np.sin(theta_a)
    y = s * np.cos(theta_a)
    return x, y, z


def make_time_unit_str(dtobj):
    """ Return a time unit string from a datetime object. """
    return "seconds since " + dtobj.strftime("%Y-%m-%dT%H:%M:%SZ")


def add_2d_latlon_axis(grid, **kwargs):
    """
    Add to Grid (in place) a 2-dimensional axes for latitude and longitude
    of every point in the y,x plane. If available conversion is done using
    basemap.pyproj, extra arguments in kwargs are passed to pyproj.Proj
    function. If not available internal implementation is used.

    Parameters
    ----------
    grid: grid object
        Cartesian grid object containing the 1d axes "x_disp", "y_disp" and
        scalar axes 'lat', 'lon'.
    kwargs: Pyproj options
        Options to be passed to Proj. If projection is not specified here it
        uses proj='aeqd' (azimuthal equidistant)

    Returns
    -------
    grid: grid object
        Cartesian grid with new axes "longitude", "latitude"

    Notes
    -----
    If Basemap is not available calculation of latitude, longitude is done
    by converting spherical azimuthal equidistant projection to latlon
    projection [1].
    It uses the mean radius of earth (6371 km)

    .. math::

        c = \\sqrt(x^2 + y^2)/R

        azi = \\arctan2(y,x) \\text{  # from east to north}

        lat = \\arcsin(\\cos(c)*\\sin(lat0)+\\sin(azi)*\\sin(c)*\\cos(lat0))

        lon = \\arctan2(\\cos(azi)*\\sin(c),\\cos(c)*\\cos(lat0)-\\sin(azi)*\\sin(c)*\\sin(lat0)) + lon0

    Where x, y are the cartesian position from the center of projection;
    lat, lon the corresponding latitude and longitude; lat0, lon0 the latitude
    and longitude of the center of the projection; R the mean radius of the
    earth (6371 km)

    References
    ----------
    .. [1] Snyder, J. P. Map Projections--A Working Manual. U. S. Geological
        Survey Professional Paper 1395, 1987, pp. 191-202.
    """
    import sys
    try:
        from mpl_toolkits.basemap import pyproj
        if 'proj' not in kwargs:
            kwargs['proj'] = 'aeqd'
        x, y = np.meshgrid(
            grid.axes["x_disp"]['data'], grid.axes["y_disp"]['data'])
        b = pyproj.Proj(lat_0=grid.axes["lat"]['data'][0],
                        lon_0=grid.axes["lon"]['data'][0], **kwargs)
        lon, lat = b(x, y, inverse=True)
    except ImportError:
        import warnings
        warnings.warn('No basemap found, using internal implementation '
                      'for converting azimuthal equidistant to latlon')
        # azimutal equidistant projetion to latlon
        R = 6371.0 * 1000.0     # radius of earth in meters.

        x, y = np.meshgrid(grid.axes["x_disp"]['data'],
                           grid.axes["y_disp"]['data'])

        c = np.sqrt(x*x + y*y) / R
        phi_0 = grid.axes["lat"]['data'] * np.pi / 180
        azi = np.arctan2(y, x)  # from east to north

        lat = np.arcsin(np.cos(c) * np.sin(phi_0) +
                        np.sin(azi) * np.sin(c) * np.cos(phi_0)) * 180 / np.pi
        lon = (np.arctan2(np.cos(azi) * np.sin(c), np.cos(c) * np.cos(phi_0) -
               np.sin(azi) * np.sin(c) * np.sin(phi_0)) * 180 /
               np.pi + grid.axes["lon"]['data'])
        lon = np.fmod(lon + 180, 360) - 180

    lat_axis = {
        'data':  lat,
        'long_name': 'Latitude for points in Cartesian system',
        'axis': 'YX',
        'units': 'degree_N',
        'standard_name': 'latitude',
    }

    lon_axis = {
        'data': lon,
        'long_name': 'Longitude for points in Cartesian system',
        'axis': 'YX',
        'units': 'degree_E',
        'standard_name': 'longitude',
    }

    grid.axes["latitude"] = lat_axis
    grid.axes["longitude"] = lon_axis
