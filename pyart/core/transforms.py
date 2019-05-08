"""
pyart.core.transforms
=====================

Transformations between coordinate systems. Routines for converting between
Cartesian/Cartographic (x, y, z), Geographic (latitude, longitude, altitude)
and antenna (azimuth, elevation, range) coordinate systems.

.. autosummary::
    :toctree: generated/

    antenna_to_cartesian
    antenna_vectors_to_cartesian
    antenna_to_cartesian_track_relative
    antenna_to_cartesian_earth_relative
    antenna_to_cartesian_aircraft_relative

    cartesian_to_geographic
    cartesian_vectors_to_geographic
    geographic_to_cartesian
    cartesian_to_geographic_aeqd
    geographic_to_cartesian_aeqd

    _interpolate_axes_edges
    _interpolate_azimuth_edges
    _interpolate_elevation_edges
    _interpolate_range_edges
    _half_angle_complex


"""

import warnings

import numpy as np
try:
    import pyproj
    _PYPROJ_AVAILABLE = True
except ImportError:
    _PYPROJ_AVAILABLE = False

from ..exceptions import MissingOptionalDependency

PI = np.pi


def antenna_to_cartesian(ranges, azimuths, elevations):
    """
    Return Cartesian coordinates from antenna coordinates.

    Parameters
    ----------
    ranges : array
        Distances to the center of the radar gates (bins) in kilometers.
    azimuths : array
        Azimuth angle of the radar in degrees.
    elevations : array
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

        z = \\sqrt{r^2+R^2+2*r*R*sin(\\theta_e)} - R

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
    theta_e = elevations * np.pi / 180.0    # elevation angle in radians.
    theta_a = azimuths * np.pi / 180.0      # azimuth angle in radians.
    R = 6371.0 * 1000.0 * 4.0 / 3.0     # effective radius of earth in meters.
    r = ranges * 1000.0                 # distances to gates in meters.

    z = (r ** 2 + R ** 2 + 2.0 * r * R * np.sin(theta_e)) ** 0.5 - R
    s = R * np.arcsin(r * np.cos(theta_e) / (R + z))  # arc length in m.
    x = s * np.sin(theta_a)
    y = s * np.cos(theta_a)
    return x, y, z


def antenna_vectors_to_cartesian(ranges, azimuths, elevations, edges=False):
    """
    Calculate Cartesian coordinate for gates from antenna coordinate vectors.

    Calculates the Cartesian coordinates for the gate centers or edges for
    all gates from antenna coordinate vectors assuming a standard atmosphere
    (4/3 Earth's radius model). See :py:func:`pyart.util.antenna_to_cartesian`
    for details.

    Parameters
    ----------
    ranges : array, 1D.
        Distances to the center of the radar gates (bins) in meters.
    azimuths : array, 1D.
        Azimuth angles of the rays in degrees.
    elevations : array, 1D.
        Elevation angles of the rays in degrees.
    edges : bool, optional
        True to calculate the coordinates of the gate edges by interpolating
        between gates and extrapolating at the boundaries. False to
        calculate the gate centers.

    Returns
    -------
    x, y, z : array, 2D
        Cartesian coordinates in meters from the center of the radar to the
        gate centers or edges.

    """
    if edges:
        if len(ranges) != 1:
            ranges = _interpolate_range_edges(ranges)
        if len(elevations) != 1:
            elevations = _interpolate_elevation_edges(elevations)
        if len(azimuths) != 1:
            azimuths = _interpolate_azimuth_edges(azimuths)
    rg, azg = np.meshgrid(ranges, azimuths)
    rg, eleg = np.meshgrid(ranges, elevations)
    return antenna_to_cartesian(rg / 1000., azg, eleg)


def _interpolate_range_edges(ranges):
    """ Interpolate the edges of the range gates from their centers. """
    edges = np.empty((ranges.shape[0] + 1, ), dtype=ranges.dtype)
    edges[1:-1] = (ranges[:-1] + ranges[1:]) / 2.
    edges[0] = ranges[0] - (ranges[1] - ranges[0]) / 2.
    edges[-1] = ranges[-1] - (ranges[-2] - ranges[-1]) / 2.
    edges[edges < 0] = 0    # do not allow range to become negative
    return edges


def _interpolate_elevation_edges(elevations):
    """ Interpolate the edges of the elevation angles from their centers. """
    edges = np.empty((elevations.shape[0]+1, ), dtype=elevations.dtype)
    edges[1:-1] = (elevations[:-1] + elevations[1:]) / 2.
    edges[0] = elevations[0] - (elevations[1] - elevations[0]) / 2.
    edges[-1] = elevations[-1] - (elevations[-2] - elevations[-1]) / 2.
    edges[edges > 180] = 180.   # prevent angles from going below horizon
    edges[edges < 0] = 0.
    return edges


def _interpolate_azimuth_edges(azimuths):
    """ Interpolate the edges of the azimuth angles from their centers. """
    edges = np.empty((azimuths.shape[0]+1, ), dtype=azimuths.dtype)
    # perform interpolation and extrapolation in complex plane to
    # account for periodic nature of azimuth angle.
    azimuths = np.exp(1.j*np.deg2rad(azimuths))

    edges[1:-1] = np.angle(azimuths[1:] + azimuths[:-1], deg=True)

    half_angle = _half_angle_complex(azimuths[0], azimuths[1])
    edges[0] = (np.angle(azimuths[0], deg=True) - half_angle) % 360.

    half_angle = _half_angle_complex(azimuths[-1], azimuths[-2])
    edges[-1] = (np.angle(azimuths[-1], deg=True) + half_angle) % 360.

    edges[edges < 0] += 360     # range from [-180, 180] to [0, 360]
    return edges


def _half_angle_complex(complex_angle1, complex_angle2):
    """
    Return half the angle between complex numbers on the unit circle.

    Parameters
    ----------
    complex_angle1, complex_angle2 : complex
        Complex numbers representing unit vectors on the unit circle

    Returns
    -------
    half_angle : float
        Half the angle between the unit vectors in degrees.

    """
    dot_product = np.real(complex_angle1 * np.conj(complex_angle2))
    if dot_product > 1:
        warnings.warn("dot_product is larger than one.")
        dot_product = 1.
    full_angle_rad = np.arccos(dot_product)
    half_angle_rad = full_angle_rad / 2.
    half_angle_deg = np.rad2deg(half_angle_rad)
    return half_angle_deg


def _interpolate_axes_edges(axes):
    """ Interpolate the edges of the axes gates from their centers. """
    edges = np.empty((axes.shape[0] + 1, ), dtype=axes.dtype)
    edges[1:-1] = (axes[:-1] + axes[1:]) / 2.
    edges[0] = axes[0] - (axes[1] - axes[0]) / 2.
    edges[-1] = axes[-1] - (axes[-2] - axes[-1]) / 2.
    return edges


def antenna_to_cartesian_track_relative(ranges, rot, roll, drift, tilt, pitch):
    """
    Calculate track-relative Cartesian coordinates from radar coordinates.

    Parameters
    ----------
    ranges : array
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
    r = ranges * 1000.0                 # distances to gates in meters.

    x = r * (np.cos(rot + roll) * np.sin(drift) * np.cos(tilt) *
             np.sin(pitch) + np.cos(drift) * np.sin(rot + roll) *
             np.cos(tilt) - np.sin(drift) * np.cos(pitch) * np.sin(tilt))
    y = r * (-1. * np.cos(rot + roll) * np.cos(drift) * np.cos(tilt) *
             np.sin(pitch) + np.sin(drift) * np.sin(rot + roll) *
             np.cos(tilt) + np.cos(drift) * np.cos(pitch) * np.sin(tilt))
    z = (r * np.cos(pitch) * np.cos(tilt) * np.cos(rot + roll) +
         np.sin(pitch) * np.sin(tilt))
    return x, y, z


def antenna_to_cartesian_earth_relative(
        ranges, rot, roll, heading, tilt, pitch):
    """
    Calculate earth-relative Cartesian coordinates from radar coordinates

    Parameters
    ----------
    ranges : array
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
    r = ranges * 1000.0                 # distances to gates in meters.

    x = r * (-1. * np.cos(rot + roll) * np.sin(heading) * np.cos(tilt) *
             np.sin(pitch) + np.cos(heading) * np.sin(rot + roll) *
             np.cos(tilt) + np.sin(heading) * np.cos(pitch) * np.sin(tilt))
    y = r * (-1. * np.cos(rot + roll) * np.cos(heading) * np.cos(tilt) *
             np.sin(pitch) - np.sin(heading) * np.sin(rot + roll) *
             np.cos(tilt) + np.cos(heading) * np.cos(pitch) * np.sin(tilt))
    z = (r * np.cos(pitch) * np.cos(tilt) * np.cos(rot + roll) +
         np.sin(pitch) * np.sin(tilt))
    return x, y, z


def antenna_to_cartesian_aircraft_relative(ranges, rot, tilt):
    """
    Calculate aircraft-relative Cartesian coordinates from radar coordinates.

    Parameters
    ----------
    ranges : array
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
    r = ranges * 1000.0                 # distances to gates in meters.
    x = r * np.cos(tilt) * np.sin(rot)
    y = r * np.sin(tilt)
    z = r * np.cos(rot) * np.cos(tilt)
    return x, y, z


def geographic_to_cartesian(lon, lat, projparams):
    """
    Geographic to Cartesian coordinate transform.

    Transform a set of Geographic coordinate (lat, lon) to a
    Cartesian/Cartographic coordinate (x, y) using pyproj or a build in
    Azimuthal equidistant projection.

    Parameters
    ----------
    lon, lat : array-like
        Geographic coordinates in degrees.
    projparams : dict or str
        Projection parameters passed to pyproj.Proj. If this parameter is a
        dictionary with a 'proj' key equal to 'pyart_aeqd' then a azimuthal
        equidistant projection will be used that is native to Py-ART and
        does not require pyproj to be installed. In this case a non-default
        value of R can be specified by setting the 'R' key to the desired
        value.

    Returns
    -------
    x, y : array-like
        Cartesian coordinates in meters unless projparams defines a value for R
        in different units.

    """
    if isinstance(projparams, dict) and projparams.get('proj') == 'pyart_aeqd':
        # Use Py-ART's Azimuthal equidistance projection
        lon_0 = projparams['lon_0']
        lat_0 = projparams['lat_0']
        if 'R' in projparams:
            R = projparams['R']
            x, y = geographic_to_cartesian_aeqd(lon, lat, lon_0, lat_0, R)
        else:
            x, y = geographic_to_cartesian_aeqd(lon, lat, lon_0, lat_0)
    else:
        # Use pyproj for the projection
        # check that pyproj is available
        if not _PYPROJ_AVAILABLE:
            raise MissingOptionalDependency(
                "PyProj is required to use geographic_to_cartesian "
                "with a projection other than pyart_aeqd but it is not "
                "installed")
        proj = pyproj.Proj(projparams)
        x, y = proj(lon, lat, inverse=False)
    return x, y


def geographic_to_cartesian_aeqd(lon, lat, lon_0, lat_0, R=6370997.):
    """
    Azimuthal equidistant geographic to Cartesian coordinate transform.

    Transform a set of geographic coordinates (lat, lon) to
    Cartesian/Cartographic coordinates (x, y) using a azimuthal equidistant
    map projection [1].

    .. math::

        x = R * k * \\cos(lat) * \\sin(lon - lon_0)

        y = R * k * [\\cos(lat_0) * \\sin(lat) -
                     \\sin(lat_0) * \\cos(lat) * \\cos(lon - lon_0)]

        k = c / \\sin(c)

        c = \\arccos(\\sin(lat_0) * \\sin(lat) +
                     \\cos(lat_0) * \\cos(lat) * \\cos(lon - lon_0))

    Where x, y are the Cartesian position from the center of projection;
    lat, lon the corresponding latitude and longitude; lat_0, lon_0 are the
    latitude and longitude of the center of the projection; R is the radius of
    the earth (defaults to ~6371 km).

    Parameters
    ----------
    lon, lat : array-like
        Longitude and latitude coordinates in degrees.
    lon_0, lat_0 : float
        Longitude and latitude, in degrees, of the center of the projection.
    R : float, optional
        Earth radius in the same units as x and y. The default value is in
        units of meters.

    Returns
    -------
    x, y : array
        Cartesian coordinates in the same units as R, typically meters.

    References
    ----------
    .. [1] Snyder, J. P. Map Projections--A Working Manual. U. S. Geological
        Survey Professional Paper 1395, 1987, pp. 191-202.

    """
    lon = np.atleast_1d(np.asarray(lon))
    lat = np.atleast_1d(np.asarray(lat))

    lon_rad = np.deg2rad(lon)
    lat_rad = np.deg2rad(lat)

    lat_0_rad = np.deg2rad(lat_0)
    lon_0_rad = np.deg2rad(lon_0)

    lon_diff_rad = lon_rad - lon_0_rad

    # calculate the arccos after ensuring all values in valid domain, [-1, 1]
    arg_arccos = (np.sin(lat_0_rad) * np.sin(lat_rad) +
                  np.cos(lat_0_rad) * np.cos(lat_rad) * np.cos(lon_diff_rad))
    arg_arccos[arg_arccos > 1] = 1
    arg_arccos[arg_arccos < -1] = -1
    c = np.arccos(arg_arccos)
    with warnings.catch_warnings():
        # division by zero may occur here but is properly addressed below so
        # the warnings can be ignored
        warnings.simplefilter("ignore", RuntimeWarning)
        k = c / np.sin(c)
    # fix cases where k is undefined (c is zero), k should be 1
    k[c == 0] = 1

    x = R * k * np.cos(lat_rad) * np.sin(lon_diff_rad)
    y = R * k * (np.cos(lat_0_rad) * np.sin(lat_rad) -
                 np.sin(lat_0_rad) * np.cos(lat_rad) * np.cos(lon_diff_rad))
    return x, y


def cartesian_to_geographic(x, y, projparams):
    """
    Cartesian to Geographic coordinate transform.

    Transform a set of Cartesian/Cartographic coordinates (x, y) to a
    geographic coordinate system (lat, lon) using pyproj or a build in
    Azimuthal equidistant projection.

    Parameters
    ----------
    x, y : array-like
        Cartesian coordinates in meters unless R is defined in different units
        in the projparams parameter.
    projparams : dict or str
        Projection parameters passed to pyproj.Proj. If this parameter is a
        dictionary with a 'proj' key equal to 'pyart_aeqd' then a azimuthal
        equidistant projection will be used that is native to Py-ART and
        does not require pyproj to be installed. In this case a non-default
        value of R can be specified by setting the 'R' key to the desired
        value.

    Returns
    -------
    lon, lat : array
        Longitude and latitude of the Cartesian coordinates in degrees.

    """
    if isinstance(projparams, dict) and projparams.get('proj') == 'pyart_aeqd':
        # Use Py-ART's Azimuthal equidistance projection
        lon_0 = projparams['lon_0']
        lat_0 = projparams['lat_0']
        if 'R' in projparams:
            R = projparams['R']
            lon, lat = cartesian_to_geographic_aeqd(x, y, lon_0, lat_0, R)
        else:
            lon, lat = cartesian_to_geographic_aeqd(x, y, lon_0, lat_0)
    else:
        # Use pyproj for the projection
        # check that pyproj is available
        if not _PYPROJ_AVAILABLE:
            raise MissingOptionalDependency(
                "PyProj is required to use cartesian_to_geographic "
                "with a projection other than pyart_aeqd but it is not "
                "installed")
        proj = pyproj.Proj(projparams)
        lon, lat = proj(x, y, inverse=True)
    return lon, lat


def cartesian_vectors_to_geographic(x, y, projparams, edges=False):
    """
    Cartesian vectors to Geographic coordinate transform.

    Transform a set of Cartesian/Cartographic coordinate vectors (x, y) to a
    geographic coordinate system (lat, lon) using pyproj or a build in
    Azimuthal equidistant projection finding the coordinates edges in
    Cartesian space if requested.

    Parameters
    ----------
    x, y : array 1D.
        Cartesian coordinate vectors in meters unless R is defined in
        different units in the projparams parameter.
    projparams : dict or str
        Projection parameters passed to pyproj.Proj. If this parameter is a
        dictionary with a 'proj' key equal to 'pyart_aeqd' then a azimuthal
        equidistant projection will be used that is native to Py-ART and
        does not require pyproj to be installed. In this case a
        non-default value of R can be specified by setting the 'R' key to the
        desired value.
    edges : bool, optional
        True to calculate the coordinates of the geographic edges by
        interpolating between Cartesian points and extrapolating at the
        boundaries. False to calculate the coordinate centers.

    Returns
    -------
    lon, lat : array
        Longitude and latitude of the Cartesian coordinates in degrees.

    """
    if edges:
        if len(x) > 1:
            x = _interpolate_axes_edges(x)
        if len(y) > 1:
            y = _interpolate_axes_edges(y)
    x, y = np.meshgrid(x, y)
    return cartesian_to_geographic(x, y, projparams)


def cartesian_to_geographic_aeqd(x, y, lon_0, lat_0, R=6370997.):
    """
    Azimuthal equidistant Cartesian to geographic coordinate transform.

    Transform a set of Cartesian/Cartographic coordinates (x, y) to
    geographic coordinate system (lat, lon) using a azimuthal equidistant
    map projection [1].

    .. math::

        lat = \\arcsin(\\cos(c) * \\sin(lat_0) +
                       (y * \\sin(c) * \\cos(lat_0) / \\rho))

        lon = lon_0 + \\arctan2(
            x * \\sin(c),
            \\rho * \\cos(lat_0) * \\cos(c) - y * \\sin(lat_0) * \\sin(c))

        \\rho = \\sqrt(x^2 + y^2)

        c = \\rho / R

    Where x, y are the Cartesian position from the center of projection;
    lat, lon the corresponding latitude and longitude; lat_0, lon_0 are the
    latitude and longitude of the center of the projection; R is the radius of
    the earth (defaults to ~6371 km). lon is adjusted to be between -180 and
    180.

    Parameters
    ----------
    x, y : array-like
        Cartesian coordinates in the same units as R, typically meters.
    lon_0, lat_0 : float
        Longitude and latitude, in degrees, of the center of the projection.
    R : float, optional
        Earth radius in the same units as x and y. The default value is in
        units of meters.

    Returns
    -------
    lon, lat : array
        Longitude and latitude of Cartesian coordinates in degrees.

    References
    ----------
    .. [1] Snyder, J. P. Map Projections--A Working Manual. U. S. Geological
        Survey Professional Paper 1395, 1987, pp. 191-202.

    """
    x = np.atleast_1d(np.asarray(x))
    y = np.atleast_1d(np.asarray(y))

    lat_0_rad = np.deg2rad(lat_0)
    lon_0_rad = np.deg2rad(lon_0)

    rho = np.sqrt(x*x + y*y)
    c = rho / R

    with warnings.catch_warnings():
        # division by zero may occur here but is properly addressed below so
        # the warnings can be ignored
        warnings.simplefilter("ignore", RuntimeWarning)
        lat_rad = np.arcsin(np.cos(c) * np.sin(lat_0_rad) +
                            y * np.sin(c) * np.cos(lat_0_rad) / rho)
    lat_deg = np.rad2deg(lat_rad)
    # fix cases where the distance from the center of the projection is zero
    lat_deg[rho == 0] = lat_0

    x1 = x * np.sin(c)
    x2 = rho*np.cos(lat_0_rad)*np.cos(c) - y*np.sin(lat_0_rad)*np.sin(c)
    lon_rad = lon_0_rad + np.arctan2(x1, x2)
    lon_deg = np.rad2deg(lon_rad)
    # Longitudes should be from -180 to 180 degrees
    lon_deg[lon_deg > 180] -= 360.
    lon_deg[lon_deg < -180] += 360.

    return lon_deg, lat_deg
