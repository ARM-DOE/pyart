"""
pyart.io.grid
=============

Reading and writing PyGrid objects.

.. autosummary::
    :toctree: generated/

    read_grid
    write_grid

    read_grid_cf
    read_grid_wrf
    grid_from_radars
    pyGrid

"""

import getpass
import socket
import datetime
from warnings import warn

import numpy as np
import netCDF4

from .netcdf import _ncvar_to_dict
from ..map.grid_mapper import map_to_grid


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
    # in the file are (1, nz, ny, nx) but the actual data should be shaped
    # (nz, ny, nx)
    ncdims = ncobj.dimensions
    field_shape = tuple([len(ncdims[i]) for i in ['nz', 'ny', 'nx']])
    field_shape_with_time = (1, ) + field_shape  # 1, nz, ny, nx on disk 
    
    # check all non-axes variables, those with the correct shape
    # are added to the field dictionary, if a wrong size field is 
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
    
    return pyGrid(fields, axes, metadata)


def write_grid(filename, pygrid, format='NETCDF4', zlib=False, **kwargs):
    """
    Write a pyGrid object to a CF and ARM standard netcdf file

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
    zlib : bool 
        True to use zlib compression in the created netCDF file.
    
    """
    ncfobj = netCDF4.Dataset(filename, 'w')

    # create the time dimension
    ncfobj.createDimension('time', None)

    # Axes dimensions
    # grab the dimensions from the first moment field
    nz, ny, nx = pygrid.fields[pygrid.fields.keys()[0]]['data'].shape
    ncfobj.createDimension('nz', nz)
    ncfobj.createDimension('ny', ny)
    ncfobj.createDimension('nx', nx)

    #axes
    #Generate the variables
    vvars = [ncfobj.createVariable(
        key, np.float, ('time', 'nz', 'ny', 'nx'),
        fill_value=pygrid.fields[key]['_FillValue'], zlib=zlib) 
            for key in pygrid.fields.keys()]

    # loop and populate attributes

    for field in pygrid.fields.keys():
        for meta in pygrid.fields[field].keys():
            if meta != 'data' and meta != '_FillValue':
                setattr(ncfobj.variables[field], meta,
                        pygrid.fields[field][meta])
    akeys = pygrid.axes.keys()
    akeys.sort()  # makes sure time comes first

    #For ARM compliance we want alt, lat and lon to be at the end

    for mkey in ['lat', 'lon', 'alt']:
        try:
            akeys.remove(mkey)
            akeys.append(mkey)
        except ValueError:
            print(mkey, " not existing")

    dims_lookup = {'time': 'time', 'x_disp': 'nx', 'y_disp': 'ny',
                   'z_disp': 'nz', 'time_end': 'time', 'time_start': 'time',
                   'lat': 'time', 'lon': 'time', 'alt': 'time'}
    avars = [ncfobj.createVariable(key, np.float, (dims_lookup[key], ))
             for key in akeys]

    # loop and populate attributes

    for axis in akeys:
        metakeys = pygrid.axes[axis].keys()

        #again, reorder to meet ARM standards..
        for mkey in ['units', 'long_name']:
            try:
                metakeys.remove(mkey)
                metakeys.insert(0, mkey)
            except ValueError:
                print(mkey, " not existing")
        for meta in metakeys:
            if meta != 'data':
                setattr(ncfobj.variables[axis],
                        meta, pygrid.axes[axis][meta])

    # global metadata
    if 'Conventions' in pygrid.metadata.keys():
        ncfobj.Conventions = pygrid.metadata['conventions']
    else:
        ncfobj.Conventions = 'CF-1.5'
    if 'process_version' in pygrid.metadata.keys():
        ncfobj.process_version = pygrid.metadata['process_version']
    for meta in pygrid.metadata.keys():
        if meta != 'history' or meta != 'process_version':
            setattr(ncfobj, meta, pygrid.metadata[meta])
    ncfobj.history = pygrid.metadata['history']

    #now populate data.. we leave this until last to speed up..
    for i in range(len(pygrid.fields.keys())):
        vvars[i][0, :, :, :] = pygrid.fields[
            pygrid.fields.keys()[i]]['data'][:, :, :]

    for i in range(len(akeys)):
        if 'shape' in dir(pygrid.axes[akeys[i]]['data']):
            avars[i][:] = pygrid.axes[akeys[i]]['data']
            print akeys[i], "is array"
        else:
            avars[i][:] = np.array([pygrid.axes[akeys[i]]['data']])
            print np.array([pygrid.axes[akeys[i]]['data']])
            print akeys[i], "is not array"

    return


def read_grid_cf(filename):
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
    return pyGrid(fields, axes, {})


def read_grid_wrf(filename):
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
    return pyGrid(fields, axes, {})


def grid_from_radars(radars, grid_shape, grid_limits, **kwargs):
    """

    Parameters
    ----------
    radars : tuple

    grid_shape :

    grid_limits :

    **kwargs :


    Returns
    -------
    pygrid : PyGrid
        A PyGrid object containing the gridded radar data.

    """
    # map the radar(s) to a cartesian grid
    grids = map_to_grid(radars, grid_shape=grid_shape,
                        grid_limits=grid_limits, **kwargs)

    # create and populate the field dictionary
    fields = {}
    first_radar = radars[0]

    for field in grids.keys():
        if field == 'ROI':
            fields['ROI'] = {
                'data': grids['ROI'],
                'standard_name': 'radius_of_influence',
                'long_name': 'Radius of influence for mapping',
                'units': 'm',
                'least_significant_digit': 1,
                'valid_min': 0.,
                'valid_max': 100000.,
                '_FillValue': 9999.0}
        else:
            fields[field] = {'data': grids[field]}
            # copy the metadata from the radar to the grid
            for key in first_radar.fields[field].keys():
                if key == 'data':
                    continue
                fields[field][key] = first_radar.fields[field][key]

    # time dictionaries
    time = {
        'data': first_radar.time['data'][0],
        'units': first_radar.time['units'],
        'calendar': first_radar.time['calendar'],
        'standard_name': first_radar.time['standard_name'],
        'long_name': 'time in seconds of volume start'}

    time_start = {
        'data': first_radar.time['data'][0],
        'units': first_radar.time['units'],
        'calendar': first_radar.time['calendar'],
        'standard_name': first_radar.time['standard_name'],
        'long_name': 'time in seconds of volume start'}

    time_end = {
        'data': first_radar.time['data'][-1],
        'units': first_radar.time['units'],
        'calendar': first_radar.time['calendar'],
        'standard_name': first_radar.time['standard_name'],
        'long_name': 'time in seconds of volume end'}

    # grid coordinate dictionaries
    nx, ny, nz = grid_shape
    (x0, x1), (y0, y1), (z0, z1) = grid_limits

    xaxis = {'data':  np.linspace(x0, x1, nx),
             'long_name': 'x-coordinate in Cartesian system',
             'axis': 'X',
             'units': 'm'}

    yaxis = {'data': np.linspace(y0, y1, ny),
             'long_name': 'y-coordinate in Cartesian system',
             'axis': 'Y',
             'units': 'm'}

    zaxis = {'data': np.linspace(z0, z1, nz),
             'long_name': 'z-coordinate in Cartesian system',
             'axis': 'Z',
             'units': 'm',
             'positive': 'up'}

    # grid origin location dictionaries
    if 'origin' in kwargs:
        lat, lon, alt = kwargs['origin']
    else:
        location = first_radar.location
        lat = location['latitude']['data']
        lon = location['longitude']['data']
        alt = location['altitude']['data']

    altorigin = {'data': alt,
                 'long_name': 'Altitude at grid origin',
                 'units': 'm'}

    latorigin = {'data': lat,
                 'long_name': 'latitude at grid origin',
                 'units': 'degrees_north'}

    lonorigin = {'data': lon,
                 'long_name': 'longitude at grid origin',
                 'units': 'degrees_east'}

    # axes dictionary
    axes = {'time': time,
            'time_start': time_start,
            'time_end': time_end,
            'z_disp': zaxis,
            'y_disp': yaxis,
            'x_disp': xaxis,
            'alt': altorigin,
            'lat': latorigin,
            'lon': lonorigin}

    # metadata dictionary
    metadata = dict(first_radar.metadata)

    # update history
    time_str = datetime.datetime.now().strftime('%d-%b-%Y, %X')
    t = (getpass.getuser(), socket.gethostname(), time_str)
    text = "Gridded by user %s on %s at %s using Py-ART's map_to_grid" % t

    if 'history' in metadata.keys():
        metadata['history'] = (metadata['history'] + '\n' + text)
    else:
        metadata['history'] = text

    return pyGrid(fields, axes, metadata)


class pyGrid:
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
