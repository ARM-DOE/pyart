""" Functions for converting date and time between various forms """

from netCDF4 import num2date


def datetime_from_radar(radar):
    """ Return a datetime for the first ray in a Radar. """
    return num2date(radar.time['data'][0], radar.time['units'])


def datetimes_from_radar(radar):
    """ Return an array of datetimes for the rays in a Radar. """
    return num2date(radar.time['data'][:], radar.time['units'])


def datetime_from_dataset(dataset):
    """ Return a datetime for the first time in a netCDF Dataset. """
    return num2date(dataset.variables['time'][0],
                    dataset.variables['time'].units)


def datetimes_from_dataset(dataset):
    """ Return an array of datetimes for the times in a netCDF Dataset. """
    return num2date(dataset.variables['time'][:],
                    dataset.variables['time'].units)


def datetime_from_grid(grid):
    """ Return a datetime for the volume start in a Grid. """
    return num2date(grid.axes['time_start']['data'][0],
                    grid.axes['time_start']['units'])
