"""
Functions for converting date and time between various forms.

"""

try:
    from cftime import num2date, date2num
except ImportError:
    from netCDF4 import num2date, date2num

EPOCH_UNITS = "seconds since 1970-01-01T00:00:00Z"


def datetime_from_radar(radar, epoch=False, **kwargs):
    """ Return a datetime for the first ray in a Radar. """
    if epoch:
        dtrad = num2date(radar.time['data'][0], radar.time['units'])
        epnum = date2num(dtrad, EPOCH_UNITS)
        return num2date(epnum, EPOCH_UNITS, **kwargs)
    else:
        return num2date(radar.time['data'][0], radar.time['units'],
                        **kwargs)


def datetimes_from_radar(radar, epoch=False, **kwargs):
    """ Return an array of datetimes for the rays in a Radar. """
    if epoch:
        dtrad = num2date(radar.time['data'][:], radar.time['units'])
        epnum = date2num(dtrad, EPOCH_UNITS)
        return num2date(epnum, EPOCH_UNITS, **kwargs)
    else:
        return num2date(radar.time['data'][:], radar.time['units'],
                        **kwargs)


def datetime_from_dataset(dataset, epoch=False, **kwargs):
    """ Return a datetime for the first time in a netCDF Dataset. """
    if epoch:
        dtdata = num2date(dataset.variables['time'][0],
                          dataset.variables['time'].units)
        epnum = date2num(dtdata, EPOCH_UNITS)
        return num2date(epnum, EPOCH_UNITS, **kwargs)
    else:
        return num2date(dataset.variables['time'][0],
                        dataset.variables['time'].units, **kwargs)


def datetimes_from_dataset(dataset, epoch=False, **kwargs):
    """ Return an array of datetimes for the times in a netCDF Dataset. """
    if epoch:
        dtdata = num2date(dataset.variables['time'][:],
                          dataset.variables['time'].units)
        epnum = date2num(dtdata, EPOCH_UNITS)
        return num2date(epnum, EPOCH_UNITS, **kwargs)
    else:
        return num2date(dataset.variables['time'][:],
                        dataset.variables['time'].units, **kwargs)


def datetime_from_grid(grid, epoch=False, **kwargs):
    """ Return a datetime for the volume start in a Grid. """
    if epoch:
        dtrad = num2date(grid.time['data'][0], grid.time['units'])
        epnum = date2num(dtrad, EPOCH_UNITS)
        return num2date(epnum, EPOCH_UNITS, **kwargs)
    else:
        return num2date(grid.time['data'][0], grid.time['units'],
                        **kwargs)
