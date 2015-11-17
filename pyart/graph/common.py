"""
pyart.graph.common
==================

Common graphing routines.

.. autosummary::
    :toctree: generated/

    dms_to_d
    corner_to_point
    ax_radius
    parse_ax
    parse_ax_fig
    parse_vmin_vmax
    parse_lon_lat
    generate_colorbar_label
    generate_field_name
    generate_radar_name
    generate_grid_name
    generate_radar_time_begin
    generate_grid_time_begin
    generate_filename
    generate_grid_filename
    generate_title
    generate_grid_title
    generate_longitudinal_level_title
    generate_latitudinal_level_title
    generate_vpt_title
    generate_ray_title
    set_limits

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import num2date

from ..io.common import dms_to_d

PI = 3.141592653589793


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
    Rc = ax_radius(point[0], units='degrees')
    # print Rc/Re
    y = ((point[0] - corner[0]) / 360.0) * PI * 2.0 * Re
    x = ((point[1] - corner[1]) / 360.0) * PI * 2.0 * Rc
    return x, y


def ax_radius(lat, units='radians'):
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


########################
# Common radar methods #
########################


def parse_ax(ax):
    """ Parse and return ax parameter. """
    if ax is None:
        ax = plt.gca()
    return ax


def parse_ax_fig(ax, fig):
    """ Parse and return ax and fig parameters. """
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()
    return ax, fig


def parse_vmin_vmax(container, field, vmin, vmax):
    """ Parse and return vmin and vmax parameters. """
    field_dict = container.fields[field]
    if vmin is None:
        if 'valid_min' in field_dict:
            vmin = field_dict['valid_min']
        else:
            vmin = -6   # default value
    if vmax is None:
        if 'valid_max' in field_dict:
            vmax = field_dict['valid_max']
        else:
            vmax = 100
    return vmin, vmax


def parse_lon_lat(grid, lon, lat):
    """ Parse lat and lon parameters """
    if lat is None:
        lat = grid.axes['lat']['data'][0]
    if lon is None:
        lon = grid.axes['lon']['data'][0]
    return lon, lat


def generate_colorbar_label(standard_name, units):
    """ Generate and return a label for a colorbar. """
    return standard_name.replace('_', ' ') + ' (' + units + ')'


def generate_field_name(container, field):
    """ Return a nice field name for a particular field. """
    if 'standard_name' in container.fields[field]:
        field_name = container.fields[field]['standard_name']
    elif 'long_name' in container.fields[field]:
        field_name = container.fields[field]['long_name']
    else:
        field_name = str(field)
    field_name = field_name.replace('_', ' ')
    field_name = field_name[0].upper() + field_name[1:]
    return field_name


def generate_radar_name(radar):
    """ Return radar name. """
    if 'instrument_name' in radar.metadata:
        return radar.metadata['instrument_name']
    else:
        return ''


def generate_grid_name(grid):
    """ Return grid name. """
    if 'instrument_name' in grid.metadata:
        return grid.metadata['instrument_name']
    else:
        return ''


def generate_radar_time_begin(radar):
    """ Return time begin in datetime instance. """
    # datetime object describing first sweep time
    times = radar.time['data'][0]
    units = radar.time['units']
    calendar = radar.time['calendar']
    return num2date(times, units, calendar)


def generate_grid_time_begin(grid):
    """ Return time begin in datetime instance. """
    # datetime object describing time
    if "time_start" in grid.axes:
        time = "time_start"
    elif 'time' in grid.axes:
        time = 'time'
    elif 'time_end' in grid.axes:
        time = 'time_end'
    times = grid.axes[time]['data'][0]
    units = grid.axes[time]['units']
    calendar = grid.axes[time]['calendar']
    return num2date(times, units, calendar)


def generate_filename(radar, field, sweep, ext='png'):
    """
    Generate a filename for a plot.

    Generated filename has form:
        radar_name_field_sweep_time.ext

    Parameters
    ----------
    radar : Radar
        Radar structure.
    field : str
        Field plotted.
    sweep : int
        Sweep plotted.
    ext : str
        Filename extension.

    Returns
    -------
    filename : str
        Filename suitable for saving a plot.

    """
    name_s = generate_radar_name(radar).replace(' ', '_')
    field_s = field.replace(' ', '_')
    time_s = generate_radar_time_begin(radar).strftime('%Y%m%d%H%M%S')
    sweep_s = str(sweep).zfill(2)
    return '%s_%s_%s_%s.%s' % (name_s, field_s, sweep_s, time_s, ext)


def generate_grid_filename(grid, field, level, ext='png'):
    """
    Generate a filename for a plot.

    Generated filename has form:
        grid_name_field_level_time.ext

    Parameters
    ----------
    grid : Grid
        Grid structure.
    field : str
        Field plotted.
    level : int
        Level plotted.
    ext : str
        Filename extension.

    Returns
    -------
    filename : str
        Filename suitable for saving a plot.

    """
    name_s = generate_grid_name(grid).replace(' ', '_')
    field_s = field.replace(' ', '_')
    time_s = generate_grid_time_begin(grid).strftime('%Y%m%d%H%M%S')
    level_s = str(level).zfill(2)
    return '%s_%s_%s_%s.%s' % (name_s, field_s, level_s, time_s, ext)


def generate_title(radar, field, sweep):
    """
    Generate a title for a plot.

    Parameters
    ----------
    radar : Radar
        Radar structure.
    field : str
        Field plotted.
    sweep : int
        Sweep plotted.

    Returns
    -------
    title : str
        Plot title.

    """
    time_str = generate_radar_time_begin(radar).isoformat() + 'Z'
    fixed_angle = radar.fixed_angle['data'][sweep]
    l1 = "%s %.1f Deg. %s " % (generate_radar_name(radar), fixed_angle,
                               time_str)
    field_name = generate_field_name(radar, field)
    return l1 + '\n' + field_name


def generate_grid_title(grid, field, level):
    """
    Generate a title for a plot.

    Parameters
    ----------
    grid : Grid
        Radar structure.
    field : str
        Field plotted.
    level : int
        Verical level plotted.


    Returns
    -------
    title : str
        Plot title.

    """
    time_str = generate_grid_time_begin(grid).isoformat() + 'Z'
    height = grid.axes["z_disp"]['data'][level] / 1000.
    l1 = "%s %.1f km %s " % (generate_grid_name(grid), height,
                             time_str)
    field_name = generate_field_name(grid, field)
    return l1 + '\n' + field_name


def generate_longitudinal_level_title(grid, field, level):
    """
    Generate a title for a plot.

    Parameters
    ----------
    grid : Grid
        Radar structure.
    field : str
        Field plotted.
    level : int
        Longitudinal level plotted.


    Returns
    -------
    title : str
        Plot title.

    """
    time_str = generate_grid_time_begin(grid).isoformat() + 'Z'
    disp = grid.axes["x_disp"]['data'][level] / 1000.
    if disp >= 0:
        direction = "east"
    else:
        direction = "west"
        disp = - disp
    l1 = "%s %.1f km %s of origin %s " % (generate_grid_name(grid), disp,
                                          direction, time_str)
    field_name = generate_field_name(grid, field)
    return l1 + '\n' + field_name


def generate_latitudinal_level_title(grid, field, level):
    """
    Generate a title for a plot.

    Parameters
    ----------
    grid : Grid
        Radar structure.
    field : str
        Field plotted.
    level : int
        Latitudinal level plotted.


    Returns
    -------
    title : str
        Plot title.

    """
    time_str = generate_grid_time_begin(grid).isoformat() + 'Z'
    disp = grid.axes["y_disp"]['data'][level] / 1000.
    if disp >= 0:
        direction = "north"
    else:
        direction = "south"
        disp = - disp
    l1 = "%s %.1f km %s of origin %s " % (generate_grid_name(grid), disp,
                                          direction, time_str)
    field_name = generate_field_name(grid, field)
    return l1 + '\n' + field_name


def generate_vpt_title(radar, field):
    """
    Generate a title for a VPT plot.

    Parameters
    ----------
    radar : Radar
        Radar structure.
    field : str
        Field plotted.

    Returns
    -------
    title : str
        Plot title.

    """
    time_str = generate_radar_time_begin(radar).isoformat() + 'Z'
    l1 = "%s %s " % (generate_radar_name(radar), time_str)
    field_name = generate_field_name(radar, field)
    return l1 + '\n' + field_name


def generate_ray_title(radar, field, ray):
    """
    Generate a title for a ray plot.

    Parameters
    ----------
    radar : Radar
        Radar structure.
    field : str
        Field plotted.
    ray : int
        Ray plotted.

    Returns
    -------
    title : str
        Plot title.

    """
    time_str = generate_radar_time_begin(radar).isoformat() + 'Z'
    l1 = "%s %s" % (generate_radar_name(radar), time_str)
    azim = radar.azimuth['data'][ray]
    elev = radar.elevation['data'][ray]
    l2 = "Ray: %i  Elevation: %.1f Azimuth: %.1f" % (ray, azim, elev)
    field_name = generate_field_name(radar, field)
    return l1 + '\n' + l2 + '\n' + field_name


def generate_az_rhi_title(radar, field, azimuth):
    """
    Generate a title for a pseudo-RHI from PPI azimuth plot.

    Parameters
    ----------
    radar : Radar
        Radar structure.
    field : str
        Field plotted.
    azimuth : float
        Azimuth plotted.

    Returns
    -------
    title : str
        Plot title.

    """
    time_str = generate_radar_time_begin(radar).isoformat() + 'Z'
    l1 = "%s %s " % (generate_radar_name(radar), time_str)
    l2 = "Azimuth: %.1f deg" % azimuth
    field_name = generate_field_name(radar, field)
    return l1 + '\n' + l2 + '\n' + field_name


def set_limits(xlim=None, ylim=None, ax=None):
    """
    Set the display limits.

    Parameters
    ----------
    xlim : tuple, optional
        2-Tuple containing y-axis limits in km. None uses default limits.
    ylim : tuple, optional
        2-Tuple containing x-axis limits in km. None uses default limits.
    ax : Axis
        Axis to adjust.  None will adjust the current axis.

    """
    if ax is None:
        ax = plt.gca()
    if ylim is not None:
        ax.set_ylim(ylim)
    if xlim is not None:
        ax.set_xlim(xlim)
