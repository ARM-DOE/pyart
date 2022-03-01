"""
Common graphing routines.

"""

import matplotlib.pyplot as plt
from netCDF4 import num2date

from ..config import get_field_colormap, get_field_limits


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


def parse_cmap(cmap, field=None):
    """ Parse and return the cmap parameter. """
    if cmap is None:
        cmap = get_field_colormap(field)
    return cmap


def parse_vmin_vmax(container, field, vmin, vmax):
    """ Parse and return vmin and vmax parameters. """
    field_dict = container.fields[field]
    field_default_vmin, field_default_vmax = get_field_limits(field)
    if vmin is None:
        if 'valid_min' in field_dict:
            vmin = field_dict['valid_min']
        else:
            vmin = field_default_vmin
    if vmax is None:
        if 'valid_max' in field_dict:
            vmax = field_dict['valid_max']
        else:
            vmax = field_default_vmax
    return vmin, vmax


def parse_lon_lat(grid, lon, lat):
    """ Parse lat and lon parameters """
    if lat is None:
        lat = grid.origin_latitude['data'][0]
    if lon is None:
        lon = grid.origin_longitude['data'][0]
    return lon, lat


def generate_colorbar_label(standard_name, units):
    """ Generate and return a label for a colorbar. """
    return str(standard_name).replace('_', ' ') + ' (' + units + ')'


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
    return num2date(times, units, calendar, only_use_cftime_datetimes=False,
                    only_use_python_datetimes=True)


def generate_radar_time_sweep(radar, sweep):
    """ Return time that a specific sweep began in a datetime instance. """
    first_ray = radar.sweep_start_ray_index['data'][sweep]
    times = radar.time['data'][first_ray]
    units = radar.time['units']
    calendar = radar.time['calendar']
    return num2date(times, units, calendar, only_use_cftime_datetimes=False,
                    only_use_python_datetimes=True)


def generate_grid_time_begin(grid):
    """ Return time begin in datetime instance. """
    times = grid.time['data'][0]
    units = grid.time['units']
    if 'calendar' in grid.time:
        calendar = grid.time['calendar']
    else:
        calendar = 'standard'
    return num2date(times, units, calendar, only_use_cftime_datetimes=False,
                    only_use_python_datetimes=True)


def generate_filename(radar, field, sweep, ext='png',
                      datetime_format='%Y%m%d%H%M%S', use_sweep_time=False):
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
    datetime_format : str
        Format of datetime (using strftime format).
    use_sweep_time : bool
        If true, the current sweep's beginning time is used.

    Returns
    -------
    filename : str
        Filename suitable for saving a plot.

    """
    name_s = generate_radar_name(radar).replace(' ', '_')
    field_s = field.replace(' ', '_')
    if use_sweep_time:
        time_s = generate_radar_time_sweep(radar, sweep).strftime(datetime_format)
    else:
        time_s = generate_radar_time_begin(radar).strftime(datetime_format)
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


def generate_title(radar, field, sweep, datetime_format=None, use_sweep_time=True):
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
    datetime_format : str
        Format of datetime (using strftime format).
    use_sweep_time : bool
        If true, the current sweep's beginning time is used.

    Returns
    -------
    title : str
        Plot title.

    """
    if use_sweep_time:
        begin_time = generate_radar_time_sweep(radar, sweep)
    else:
        begin_time = generate_radar_time_begin(radar)
    if datetime_format:
        time_str = begin_time.strftime(datetime_format)
    else:
        time_str = begin_time.isoformat() + 'Z'
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
    height = grid.z['data'][level] / 1000.
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
    time_str = generate_grid_time_begin(grid).strftime('%Y-%m-%dT%H:%M:%SZ')
    disp = grid.x['data'][level] / 1000.
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
    disp = grid.y['data'][level] / 1000.
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
    l2 = "Ray: %i  Elevation: %.1f Azimuth: %.1f" % (ray, elev, azim)
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
