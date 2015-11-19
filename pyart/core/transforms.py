"""
pyart.core.transforms
=====================

Transformations between coordinate systems. Routines for converting between
Cartesian (x, y, z), cartographic (latitude, longitude, altitude) and antenna
(azimuth, elevation, range) coordinate systems.

.. autosummary::
    :toctree: generated/

    antenna_to_cartesian
    antenna_vectors_to_cartesian
    antenna_to_cartesian_track_relative
    antenna_to_cartesian_earth_relative
    antenna_to_cartesian_aircraft_relative
    add_2d_latlon_axis
    corner_to_point
    _interpolate_axes_edges
    _interpolate_azimuth_edges
    _interpolate_elevation_edges
    _interpolate_range_edges
    _half_angle_complex
    _ax_radius


"""

import numpy as np

PI = np.pi


def antenna_to_cartesian(ranges, azimuths, elevations, debug=False):
    """
    Return cartesian coordinates from antenna coordinates.

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
        between gates and extrapolating at the boundaries.  False to
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
        Complex numbers represeting unit vectors on the unit circle

    Returns
    -------
    half_angle : float
        Half the angle between the unit vectors in degrees.

    """
    dot_product = np.real(complex_angle1 * np.conj(complex_angle2))
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


def corner_to_point(corner, point):
    """
    Return the x, y distances in meters from a corner to a point.

    Assumes a spherical earth model.

    Parameters
    ----------
    corner : (float, float)
        Latitude and longitude in degrees of the corner.
    point : (float, float)
        Latitude and longitude in degrees of the point.

    Returns
    -------
    x, y : floats
        Distances from the corner to the point in meters.

    """
    Re = 6371.0 * 1000.0
    Rc = _ax_radius(point[0], units='degrees')
    # print Rc/Re
    y = ((point[0] - corner[0]) / 360.0) * PI * 2.0 * Re
    x = ((point[1] - corner[1]) / 360.0) * PI * 2.0 * Rc
    return x, y


def _ax_radius(lat, units='radians'):
    """
    Return the radius of a constant latitude circle for a given latitude.

    Parameters
    ----------
    lat : float
        Latitude at which to calculate constant latitude circle (parallel)
        radius.
    units : 'radians' or 'degrees'
        Units of lat, either 'radians' or 'degrees'.

    Returns
    -------
    R : float
        Radius in meters of a constant latitude circle (parallel).

    """
    Re = 6371.0 * 1000.0
    if units == 'degrees':
        const = PI / 180.0
    else:
        const = 1.0
    R = Re * np.sin(PI / 2.0 - abs(lat * const))
    return R
