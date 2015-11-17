"""
pyart.io.common
===============

Input/output routines common to many file formats.

.. autosummary::
    :toctree: generated/

    prepare_for_read
    dms_to_d
    stringarray_to_chararray
    _test_arguments
    make_time_unit_str
    add_2d_latlon_axis

"""

import bz2
import gzip

import numpy as np
import netCDF4


def prepare_for_read(filename):
    """
    Return a file like object read for reading.

    Open a file for reading in binary mode with transparent decompression of
    Gzip and BZip2 files.  The resulting file-like object should be closed.

    Parameters
    ----------
    filename : str or file-like object
        Filename or file-like object which will be opened.  File-like objects
        will not be examined for compressed data.

    Returns
    -------
    file_like : file-like object
        File like object from which data can be read.

    """
    # if a file-like object was provided, return
    if hasattr(filename, 'read'):   # file-like object
        return filename

    # look for compressed data by examining the first few bytes
    fh = open(filename, 'rb')
    magic = fh.read(3)
    fh.close()

    if magic.startswith(b'\x1f\x8b'):
        return gzip.GzipFile(filename, 'rb')

    if magic.startswith(b'BZh'):
        return bz2.BZ2File(filename, 'rb')

    return open(filename, 'rb')


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


def _test_arguments(dic):
    """ Issue a warning if receive non-empty argument dict """
    if dic:
        import warnings
        warnings.warn('Unexpected arguments: %s' % dic.keys())


def make_time_unit_str(dtobj):
    """ Return a time unit string from a datetime object. """
    return "seconds since " + dtobj.strftime("%Y-%m-%dT%H:%M:%SZ")


def add_2d_latlon_axis(grid, **kwargs):
    """
    Add the latitude and longitude for grid points in the y, x plane.

    Adds a **latitude** and **longitude** dictionary to the axes attribute
    of a provided grid.  Addition is done in-place, nothing is returned from
    this function.  These dictionaries contain 2D arrays which specify the
    latitude and longitude of every point in the y, x plane.

    If available, the conversion is done using basemap.pyproj, extra arguments
    are passed to pyproj.Proj. If not available, an internal spherical
    azimuthal equidistant transformation is is used.

    Parameters
    ----------
    grid: grid object
        Cartesian grid object containing the 1d axes "x_disp", "y_disp" and
        scalar axes 'lat', 'lon'.
    kwargs: Pyproj options
        Options to be passed to Proj. If projection is not specified here it
        uses proj='aeqd' (azimuthal equidistant)

    Notes
    -----
    If Basemap is not available, calculation of the latitude, longitude is
    done using a azimuthal equidistant projection projection [1].
    It uses the mean radius of earth (6371 km)

    .. math::

        c = \\sqrt(x^2 + y^2)/R

        azi = \\arctan2(y,x) \\text{  # from east to north}

        lat = \\arcsin(\\cos(c)*\\sin(lat0)+\\sin(azi)*\\sin(c)*\\cos(lat0))

        lon = \\arctan2(\\cos(azi)*\\sin(c),\\cos(c)*\\cos(lat0)-
                        \\sin(azi)*\\sin(c)*\\sin(lat0)) + lon0

    Where x, y are the Cartesian position from the center of projection;
    lat, lon the corresponding latitude and longitude; lat0, lon0 the latitude
    and longitude of the center of the projection; R the mean radius of the
    earth (6371 km)

    References
    ----------
    .. [1] Snyder, J. P. Map Projections--A Working Manual. U. S. Geological
        Survey Professional Paper 1395, 1987, pp. 191-202.

    """
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
