"""
pyart.io.grid
=============

Reading and writing PyGrid objects.

.. autosummary::
    :toctree: generated/

    read_grid
    write_grid

    _read_grid_cf
    _read_grid_wrf

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    PyGrid

"""

from warnings import warn

import netCDF4

from .netcdf import _ncvar_to_dict, _create_ncvar


def read_grid(filename):
    """
    Read a netCDF grid file

    Parameters
    ----------
    filename : str
        Filename of netCDF grid file to read

    Returns
    -------
    pygrid : PyGrid
        PyGrid object containing Grid data.

    """
    ncobj = netCDF4.Dataset(filename, 'r')
    ncvars = ncobj.variables

    # metadata
    metadata = dict([(k, getattr(ncobj, k)) for k in ncobj.ncattrs()])

    # axes
    axes_keys = ['time', 'time_start', 'time_end',
                 'z_disp', 'y_disp', 'x_disp',
                 'alt', 'lat', 'lon']
    axes = dict((k, _ncvar_to_dict(ncvars[k])) for k in axes_keys)

    # read in the fields
    # determine the correct shape of the fields
    # ARM standard required a time dimension, so the shape of the fields
    # in the file are (1, nz, ny, nx) but the field data should be shaped
    # (nz, ny, nx) in the PyGrid object
    ncdims = ncobj.dimensions
    field_shape = tuple([len(ncdims[i]) for i in ['nz', 'ny', 'nx']])
    field_shape_with_time = (1, ) + field_shape  # 1, nz, ny, nx on disk

    # check all non-axes variables, those with the correct shape
    # are added to the field dictionary, if a wrong sized field is
    # detected a warning is raised
    field_keys = [k for k in ncvars.keys() if k not in axes_keys]
    fields = {}
    for field in field_keys:
        field_dic = _ncvar_to_dict(ncvars[field])
        if field_dic['data'].shape == field_shape_with_time:
            field_dic['data'].shape = field_shape
            fields[field] = field_dic
        else:
            bad_shape = field_dic['data'].shape
            warn('Field %s skipped, incorrect shape %s', (field, bad_shape))

    return PyGrid(fields, axes, metadata)


def write_grid(filename, pygrid, format='NETCDF4'):
    """
    Write a PyGrid object to a CF-1.5 and ARM standard netcdf file

    Parameters
    ----------
    filename : str
        Filename to save pygrid to.
    pygrid : PyGrid
        PyGrid object to write.
    format : str, optional
        NetCDF format, one of 'NETCDF4', 'NETCDF4_CLASSIC',
        'NETCDF3_CLASSIC' or 'NETCDF3_64BIT'. See netCDF4 documentation for
        details.

    """
    ncobj = netCDF4.Dataset(filename, 'w')

    # create the time dimension
    ncobj.createDimension('time', None)

    # create additional dimensions
    grid_shape = pygrid.fields[pygrid.fields.keys()[0]]['data'].shape
    nz, ny, nx = grid_shape
    ncobj.createDimension('nz', nz)
    ncobj.createDimension('ny', ny)
    ncobj.createDimension('nx', nx)

    # field variables
    for field, field_dic in pygrid.fields.iteritems():
        # append 1, to the shape of all data to indicate the time var.
        field_dic['data'].shape = (1, ) + field_dic['data'].shape
        _create_ncvar(field_dic, ncobj, field, ('time', 'nz', 'ny', 'nx'))
        field_dic['data'].shape = field_dic['data'].shape[1:]

    # axes variables
    _create_ncvar(pygrid.axes['time'], ncobj, 'time', ('time', ))
    _create_ncvar(pygrid.axes['time_end'], ncobj, 'time_end', ('time', ))
    _create_ncvar(pygrid.axes['time_start'], ncobj, 'time_start', ('time', ))
    _create_ncvar(pygrid.axes['x_disp'], ncobj, 'x_disp', ('nx', ))
    _create_ncvar(pygrid.axes['y_disp'], ncobj, 'y_disp', ('ny', ))
    _create_ncvar(pygrid.axes['z_disp'], ncobj, 'z_disp', ('nz', ))
    _create_ncvar(pygrid.axes['lat'], ncobj, 'lat', ('time', ))
    _create_ncvar(pygrid.axes['lon'], ncobj, 'lon', ('time', ))
    _create_ncvar(pygrid.axes['alt'], ncobj, 'alt', ('time', ))

    # metadata
    for k, v in pygrid.metadata.iteritems():
        setattr(ncobj, k, v)
    if 'Conventions' not in pygrid.metadata.keys():
        ncobj.Conventions = 'CF-1.5'

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
    pygrid : PyGrid
        PyGrid object containing data.

    Notes
    -----
    This function does only the most basic variable checking.  The resulting
    PyGrid object is most likely not writable.

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
    return PyGrid(fields, axes, {})


def _read_grid_wrf(filename):
    """
    Read a WRF netCDF file containing a grid.

    Parameters
    ----------
    filename : str
        Filename of the WRF netCDF file.

    Returns
    -------
    pygrid : PyGrid
        PyGrid object containing data.

    Notes
    -----
    This function does only the most basic variable checking.  The resulting
    PyGrid object is most likely not writable.

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
    return PyGrid(fields, axes, {})


class PyGrid:
    """
    An object for holding gridded Radar data.

    Parameters
    ----------
    fields : dict

    axes : dict

    metadata : dict


    Attributes
    ----------
    fields: dict

    axes: dict

    metadata: dict

    """
    def __init__(self, fields, axes, metadata):
        """ initalize """
        self.fields = fields
        self.metadata = metadata
        self.axes = axes
        return

    def write(self, filename, format='NETCDF4'):
        """
        Write the the PyGrid object to a netcdf file.

        Parameters
        ----------
        filename : str
            Filename to save to.
        format : str, optional
            NetCDF format, one of 'NETCDF4', 'NETCDF4_CLASSIC',
            'NETCDF3_CLASSIC' or 'NETCDF3_64BIT'. See netCDF4 documentation
            fordetails.

        """
        write_grid(filename, self, format=format)
