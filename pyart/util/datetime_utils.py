""" Functions for converting date and time between various forms """

from netCDF4 import num2date, date2num


def datetime_from_radar(radar, epoch=False):
    """ Return a datetime for the first ray in a Radar. """
    if epoch:
        dtrad = num2date(radar.time['data'][0], radar.time['units'])
        epnum = date2num(dtrad, epoch_units())
        return num2date(epnum, epoch_units())
    else:
        return num2date(radar.time['data'][0], radar.time['units'])


def datetimes_from_radar(radar, epoch=False):
    """ Return an array of datetimes for the rays in a Radar. """
    if epoch:
        dtrad = num2date(radar.time['data'][:], radar.time['units'])
        epnum = date2num(dtrad, epoch_units())
        return num2date(epnum, epoch_units())
    else:
        return num2date(radar.time['data'][:], radar.time['units'])


def datetime_from_dataset(dataset, epoch=False):
    """ Return a datetime for the first time in a netCDF Dataset. """
    if epoch:
        dtdata = num2date(dataset.variables['time'][0],
                         dataset.variables['time'].units)
        epnum = date2num(dtdata, epoch_units())
        return num2date(epnum, epoch_units())
    else:
        return num2date(dataset.variables['time'][0],
                    dataset.variables['time'].units)


def datetimes_from_dataset(dataset, epoch=False):
    """ Return an array of datetimes for the times in a netCDF Dataset. """
    if epoch:
        dtdata = num2date(dataset.variables['time'][:],
                         dataset.variables['time'].units)
        epnum = date2num(dtdata, epoch_units())
        return num2date(epnum, epoch_units())
    else:
        return num2date(dataset.variables['time'][:],
                    dataset.variables['time'].units)


def datetime_from_grid(grid, epoch=False):
    """ Return a datetime for the volume start in a Grid. """
    if epoch:
        dtrad = num2date(grid.time['data'][0], grid.time['units'])
        epnum = date2num(dtrad, epoch_units())
        return num2date(epnum, epoch_units())
    else:
        return num2date(grid.time['data'][0], grid.time['units'])

def epoch_units():
    """ Return the units used in converting to datetime in the
    standard Epoch units. """
    return "seconds since 1970-01-01T00:00:00Z"
