"""
pyart.io.grid
=============

Reading and writing Grid objects.

.. autosummary::
    :toctree: generated/

    read_grid
    write_grid

    _read_grid_cf
    _read_grid_wrf

"""

from warnings import warn

import numpy as np
import netCDF4
import datetime

from ..core.grid import Grid
from .cfradial import _ncvar_to_dict, _create_ncvar


def read_grid(filename, exclude_fields=None):
    """
    Read a netCDF grid file

    Parameters
    ----------
    filename : str
        Filename of NetCDF grid file to read.

    Other Parameters
    ----------------
    exclude_fields : list
        A list of fields to exclude from the grid object.

    Returns
    -------
    grid : Grid
        Grid object containing gridded data.

    """

    if exclude_fields is None:
        exclude_fields = []

    ncobj = netCDF4.Dataset(filename, mode='r')

    # metadata
    metadata = dict([(k, getattr(ncobj, k)) for k in ncobj.ncattrs()])

    # fields are any thing with full shape
    # axes are anything with parcial shape
    # first search for fields
    # determine the correct shape of the fields
    # ARM standard requires the left-most dimension to be time, so the shape
    # of the fields in the file is (1, nz, ny, nx) but the field data should
    # be shaped (nz, ny, nx) in the Grid object
    dim_keys = ['nz', 'ny', 'nx', 'z', 'y', 'x']
    field_shape = tuple([len(ncobj.dimensions[k]) for k in dim_keys
                         if k in ncobj.dimensions])
    field_shape_with_time = (1, ) + field_shape
    dimensions_allowed = set(dim_keys + ['time'])  # dimensions allowed
    dimension_to_axis = {'nz': 'Z', 'ny': 'Y', 'nx': 'X',
                         'z': 'Z', 'y': 'Z', 'x': 'X','time': 'T'}

    # check all variables, those with the correct shape
    # are added to the field dictionary, if a wrong sized add to axes
    var_keys = [k for k in ncobj.variables if k not in exclude_fields]
    axes_key = []
    fields = {}
    for field in var_keys:
        field_dic = _ncvar_to_dict(ncobj.variables[field])
        if field_dic['data'].shape == field_shape_with_time:
            field_dic['data'].shape = field_shape
            fields[field] = field_dic
        else:
            axes_key.append(field)

    #search for axes with correct dimensions, if wrong warning is raised
    axes = {}
    for axis in axes_key:
        dimensions = ncobj.variables[axis].dimensions
        if dimensions_allowed.issuperset(dimensions):
            axis_dic = _ncvar_to_dict(ncobj.variables[axis])
            if 'axis' not in axis_dic:
                form = "".join([dimension_to_axis[k] for k in dimensions])
                axis_dic['axis'] = form
            axes[axis] = axis_dic
        else:
            warn('Field %s skipped due to incorrect shape' % (axis))

    return Grid(fields, axes, metadata)


def write_grid(filename, grid, format='NETCDF4', arm_time_variables=False):
    """
    Write a Grid object to a CF-1.5 and ARM standard netcdf file

    Parameters
    ----------
    filename : str
        Filename to save grid to.
    grid : Grid
        Grid object to write.
    format : str, optional
        NetCDF format, one of 'NETCDF4', 'NETCDF4_CLASSIC',
        'NETCDF3_CLASSIC' or 'NETCDF3_64BIT'. See netCDF4 documentation for
        details.
    arm_time_variables : bool
        True to write the ARM standard time variables base_time and
        time_offset. False will not write these variables.

    """
    ncobj = netCDF4.Dataset(filename, mode='w', format=format)

    # create the time dimension
    ncobj.createDimension('time', None)

    # create additional dimensions
    grid_shape = grid.fields[grid.fields.keys()[0]]['data'].shape
    nz, ny, nx = grid_shape
    ncobj.createDimension('nz', nz)
    ncobj.createDimension('ny', ny)
    ncobj.createDimension('nx', nx)

    # axes variables
    # allow writing all elements of axes to variables
    # so that the user may change the output netcdf by changing direct
    # the grid object at the same time as Grid created by Pyart are still
    # written in ARM format.
    # We uses the 'axis' attribute to determine the dimension, one however
    # should discuss adding this information direct in the grid object
    for key in grid.axes.keys():
        if "axis" in grid.axes[key].keys():
            if 'X' == grid.axes[key]["axis"]:
                dimension = ('nx',)
            elif 'Y' == grid.axes[key]["axis"]:
                dimension = ('ny',)
            elif 'Z' == grid.axes[key]["axis"]:
                dimension = ('nz',)
            elif 'YX' == grid.axes[key]["axis"]:
                dimension = ('ny', 'nx')
            else:
                dimension = ()
        else:
            dimension = ('time',)
        _create_ncvar(grid.axes[key], ncobj, key, dimension)

    # create ARM time variables base_time and time_offset, if requested
    if arm_time_variables:
        time = grid.axes['time']
        dt = netCDF4.num2date(time['data'][0], time['units'])
        td = dt - datetime.datetime.utcfromtimestamp(0)
        td = td.seconds + td.days * 24 * 3600

        base_time = {
            'data': np.array([td], dtype=np.int32),
            'string': dt.strftime('%d-%b-%Y,%H:%M:%S GMT'),
            'units': 'seconds since 1970-1-1 0:00:00 0:00',
            'ancillary_variables': 'time_offset',
            'long_name': 'Base time in Epoch',
        }
        _create_ncvar(base_time, ncobj, 'base_time', ())

        time_offset = {
            'data': np.array(time['data'], dtype=np.float64),
            'long_name': 'Time offset from base_time',
            'units': time['units'].replace('T', ' ').replace('Z', ''),
            'ancillary_variables': 'time_offset',
            'calendar': 'gregorian',
        }
        _create_ncvar(time_offset, ncobj, 'time_offset', ('time', ))

    # field variables
    for field, field_dic in grid.fields.iteritems():
        # append 1, to the shape of all data to indicate the time var.
        field_dic['data'].shape = (1, ) + field_dic['data'].shape
        _create_ncvar(field_dic, ncobj, field, ('time', 'nz', 'ny', 'nx'))
        field_dic['data'].shape = field_dic['data'].shape[1:]

    # metadata
    for k, v in grid.metadata.iteritems():
        setattr(ncobj, k, v)

    ncobj.close()

    return


def _read_grid_cf(filename):
    """
    Read a CF compliant netCDF file containing a grid.

    Parameters
    ----------
    filename : str
        Filename of the netCDF file.

    Returns
    -------
    grid : Grid
        Grid object containing gridded data.

    Notes
    -----
    This function does only the most basic variable checking.  The resulting
    Grid object is most likely not writable.

    """
    ncobj = netCDF4.Dataset(filename)
    ncvars = ncobj.variables
    fields = {}
    axes = {}
    for var in ncvars:
        if len(ncvars[var].shape) > 1:
            # dimensionality of 2+ are fields variables
            fields[var] = _ncvar_to_dict(ncvars[var])
        else:
            # dimensionality of 1 are axes variables
            axes[var] = _ncvar_to_dict(ncvars[var])
    return Grid(fields, axes, {})


def _read_grid_wrf(filename):
    """
    Read a WRF netCDF file containing a grid.

    Parameters
    ----------
    filename : str
        Filename of the WRF netCDF file.

    Returns
    -------
    grid : Grid
        Grid object containing data.

    Notes
    -----
    This function does only the most basic variable checking.  The resulting
    Grid object is most likely not writable.

    """
    ncobj = netCDF4.Dataset(filename)
    ncvars = ncobj.variables
    fields = {}
    axes = {}
    for var in ncvars:
        if len(ncvars[var].shape) > 1:
            # dimensionality of 2+ are fields variables
            fields[var] = _ncvar_to_dict(ncvars[var])
        else:
            # dimensionality of 1 are axes variables
            axes[var] = _ncvar_to_dict(ncvars[var])
    return Grid(fields, axes, {})
