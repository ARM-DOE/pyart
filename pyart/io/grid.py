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

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    Grid

"""

import numpy as np
import netCDF4

from warnings import warn

from .cfradial import _ncvar_to_dict, _create_ncvar
from ..config import get_fillvalue
from ..io import grid_qc


def read_grid(filename, exclude_fields=[]):
    """
    Read a netCDF grid file

    Parameters
    ----------
    filename : str
        Filename of NetCDF grid file to read.
        
    Optional parameters
    -------------------
    exclude_fields : list
        A list of field names to exclude from Grid.

    Returns
    -------
    grid : Grid
        Grid object containing gridded data.

    """
    
    ncobj = netCDF4.Dataset(filename, mode='r')

    # Metadata
    metadata = dict([(k, getattr(ncobj, k)) for k in ncobj.ncattrs()])

    # Axes
    axes_keys = ['time', 'time_start', 'time_end', 'base_time',
                 'time_offset', 'z_disp', 'y_disp', 'x_disp',
                 'alt', 'lat', 'lon', 'z', 'y', 'x', 'lev']
    axes = dict((k, _ncvar_to_dict(ncobj.variables[k])) for k in axes_keys if
                k in ncobj.variables)
    
    # Determine the correct shape of the fields. ARM standard required a time
    # dimension, so the shape of the fields in the file are (1, nz, ny, nx)
    # but the field data should be shaped (nz, ny, nx) in the Grid object
    dim_keys = ['nz', 'ny', 'nx', 'z', 'y', 'x']
    field_shape = tuple([len(ncobj.dimensions[i]) for i in dim_keys if
                         i in ncobj.dimensions])
    field_shape_with_time = (1, ) + field_shape

    # Check all non-axes variables, those with the correct shape
    # are added to the field dictionary, if a wrong sized field is
    # detected a warning is raised
    field_keys = [k for k in ncobj.variables if k not in axes_keys and
                  k not in exclude_fields]
    fields = {}
    for field in field_keys:
        field_dic = _ncvar_to_dict(ncobj.variables[field])
        if field_dic['data'].shape == field_shape_with_time:
            field_dic['data'].shape = field_shape
            fields[field] = field_dic
        else:
            bad_shape = field_dic['data'].shape
            warn('Field %s skipped due to incorrect shape' %field)

    return Grid(fields, axes, metadata)


def write_grid(filename, grid, format='NETCDF4'):
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
    if 'base_time' in grid.axes:
        _create_ncvar(grid.axes['base_time'], ncobj, 'base_time', ())
    if 'time_offset' in grid.axes:
        _create_ncvar(grid.axes['time_offset'], ncobj, 'time_offset',
                      ('time',))
    _create_ncvar(grid.axes['time'], ncobj, 'time', ('time', ))
    _create_ncvar(grid.axes['time_end'], ncobj, 'time_end', ('time', ))
    _create_ncvar(grid.axes['time_start'], ncobj, 'time_start', ('time', ))
    _create_ncvar(grid.axes['x_disp'], ncobj, 'x_disp', ('nx', ))
    _create_ncvar(grid.axes['y_disp'], ncobj, 'y_disp', ('ny', ))
    _create_ncvar(grid.axes['z_disp'], ncobj, 'z_disp', ('nz', ))
    _create_ncvar(grid.axes['lat'], ncobj, 'lat', ('time', ))
    _create_ncvar(grid.axes['lon'], ncobj, 'lon', ('time', ))
    _create_ncvar(grid.axes['alt'], ncobj, 'alt', ('time', ))

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
        Frid object containing data.

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


class Grid:
    """
    An object for holding gridded Radar data.

    Parameters
    ----------
    fields : dict
        Dictionary of field dictionaries.
    axes : dict
        Dictionary of axes dictionaries.
    metadata : dict
        Dictionary of metadata.

    Attributes
    ----------
    fields: dict
        Dictionary of field dictionaries.
    axes: dict
        Dictionary of axes dictionaries.
    metadata: dict
        Dictionary of metadata.

    """
    
    def __init__(self, fields, axes, metadata):
        """ Initalize object. """
        
        self.fields = fields
        self.metadata = metadata
        self.axes = axes
        
        return

    def write(self, filename, format='NETCDF4'):
        """
        Write the the Grid object to a netcdf file.

        Parameters
        ----------
        filename : str
            Filename to save to.
            
        Optional parameters
        -------------------
        format : str
            Format can be one of 'NETCDF4', 'NETCDF4_CLASSIC',
            'NETCDF3_CLASSIC' or 'NETCDF3_64BIT'
            
        """
        write_grid(filename, self, format=format)
        
    def add_field(self, field_name, field_dic):
        """
        Add field to object.
        
        Parameters
        ----------
        field_name : str
            Name of the field to add.
        field_dic : dic
            Dictionary containing the field data and metadata.
        """
        
        nz, ny, nx = self.fields[self.fields.keys()[0]]['data'].shape
        
        if 'data' not in field_dic:
            raise KeyError('Field dictionary must contain a "data" key')
        if field_name in self.fields:
            raise ValueError('A field with name %s already exists' %field_name)
        if field_dic['data'].shape != (nz, ny, nx):
            raise ValueError('Field has invalid shape')
        
        self.fields[field_name] = field_dic
        
        return
    
    def despeckle_field(self, field, window_size=6, noise_ratio=85.0,
                        proc=1, fill_value=None):
        """
        Remove outliers in gridded data.
        
        Parameters
        ----------
        field : str
            Name of field to despeckle.
        
        Optional parameters
        -------------------
        window_size : int
            The size (length) of the window in the x-, y-, and z-dimension.
            This will be a square box for interior grid points and a
            rectangular box for boundary grid points.
        noise_ratio : float
            The ratio (as a percent) for when a grid point can be classified
            as noise. For example, if the noise ratio is 85, then when more
            than 85% of the grid points in the box are noise, then that grid
            point will also be classified as noise.
        proc : int
            Number of processors requested.
        fill_value : float
            Missing value used to signify bad data points. A value of None
            will use the default fill value as defined in the Py-ART
            configuration file
        """
        
        # Get fill value
        if fill_value is None:
            fill_value = get_fillvalue()
            
        # Get data
        data = self.fields[field]['data']
        data = np.ma.filled(data, fill_value).astype(np.float64)
        
        data_qc = grid_qc.despeckle(data, window_size=window_size,
                                    noise_ratio=noise_ratio, proc=proc,
                                    fill_value=fill_value)
        
        data_qc = np.ma.masked_equal(data_qc, fill_value)
        
        self.fields[field]['data'] = data_qc
        
        return
