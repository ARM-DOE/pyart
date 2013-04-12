"""
pyart.io.grid
=============

A general gridded data object for radar data.

.. autosummary::
    :toctree: generated/

    grid2
    ncvar_to_field
    pyGrid

"""

import numpy as np
import pyart.map.ballsy as ballsy
from pyart.graph import radar_display
import getpass
import socket
import datetime as dt


def radar_coords_to_cart(rng, az, ele, debug=False):
    """
    Assumes standard atmosphere, ie R=4Re/3
    """
    Re = 6371.0 * 1000.0
    #h=(r^2 + (4Re/3)^2 + 2r(4Re/3)sin(ele))^1/2 -4Re/3
    #s=4Re/3arcsin(rcos(ele)/(4Re/3+h))
    p_r = 4.0 * Re / 3.0
    rm = rng
    z = (rm ** 2 + p_r ** 2 + 2.0 * rm * p_r *
         np.sin(ele * np.pi / 180.0)) ** 0.5 - p_r
    #arc length
    s = p_r * np.arcsin(rm * np.cos(ele * np.pi / 180.) / (p_r + z))
    if debug:
        print "Z=", z, "s=", s
    y = s * np.cos(az * np.pi / 180.0)
    x = s * np.sin(az * np.pi / 180.0)
    return x, y, z


def grid2(radars, **kwargs):
    """
    map tuples of Py-Radar objects to a cartesian grid..

    Usage: (xr,yr,zr), (nx,ny,nz), grids2 = grid2(
        (radar1, radar2, radarn,  nx=81, ny=81, nz=69,
         zr=(0.,17000.), yr=(-20000., 20000.), xr=(-30000., 20000.),
         origin=(lat, lon, height_m))

    """
    print(kwargs.get('h_factor', 'h_factor not set'))
    nx, ny, nz = kwargs.get('nxyz', (81, 81, 69))
    xr, yr, zr = kwargs.get(
        'xyzr', ((-30000., 20000), (-20000., 20000.), (0., 17000.)))
    toa = kwargs.get('toa', 17000.0)
    # need to just set to the radar..
    first_radar_location = radars[0].location
    cf_lat, cf_lon, cf_alt = kwargs.get(
        'origin', (first_radar_location['latitude']['data'],
                   first_radar_location['longitude']['data'],
                   first_radar_location['altitude']['data']))

    # initialize blank arrays to be filled by the radar gate coordinates

    x = np.array([])
    y = np.array([])
    z = np.array([])

    # parameters to be mapped, reflectivity should always be first!
    parms = kwargs.get(
        'params', ['corrected_reflectivity_horizontal',
                   'reflectivity_horizontal',
                   'rain_rate_A',
                   'recalculated_diff_phase'])
    data = dict([(parms[i], np.ma.array([])) for i in range(len(parms))])

    # Loop over radars and append them to the lists to be used to call the
    # mapping object

    for myradar in radars:

        # calculate radar offset to the origin

        displacement = radar_display.corner_to_point(
            [cf_lat, cf_lon], [myradar.location['latitude']['data'],
                               myradar.location['longitude']['data']])

        # meshgrid up azimuths, ranges and elevations
        rg, azg = np.meshgrid(myradar.range['data'], myradar.azimuth['data'])
        rg, eleg = np.meshgrid(myradar.range['data'],
                               myradar.elevation['data'])

        # Calculate cartesian locations of gates

        xdash, ydash, zdash = radar_coords_to_cart(rg, azg, eleg)
        del rg, azg  # housekeeping

        #only use gates below toa

        within_sensible = np.where(zdash.flatten() < toa)[0]
        x_disp = displacement[0]
        y_disp = displacement[1]
        # append geolocation..
        x = np.append(x, (xdash + x_disp).flatten()[within_sensible])
        y = np.append(y, (ydash + y_disp).flatten()[within_sensible])
        z = np.append(z, zdash.flatten()[within_sensible])

        #append gate variables..

        for var in parms:
            data[var] = np.append(data[var], myradar.fields[
                var]['data'].flatten()[within_sensible])

    # find NaNs and crazy reflectivities
    is_finite = np.isfinite(data[parms[0]])
    is_not_crazy = data[parms[0]] < 190.0
    is_not_masked = np.logical_not(data[parms[0]].mask)

    where_the_data_is_good = np.where(np.logical_and(np.logical_and(
        is_not_masked, is_finite), is_not_crazy))[0]
    print(len(where_the_data_is_good))
    print(len(data[parms[0]]))
    # Create a meshgrid(cube) to allow calculation of radii of influence

    zg, yg, xg = np.mgrid[zr[0]:zr[1]:nz * 1j, yr[0]:yr[1]:ny * 1j,
                          xr[0]:xr[1]:nx * 1j]

    # Virtual beam width and beam spacing

    nb = kwargs.get('nb', 1.5)
    bsp = kwargs.get('bsp', 1.0)

    # Query radius of influence, flattened

    if 'qrf' in kwargs.keys():
        qrf = kwargs['qrf'](xg, yg, zg).flatten()
    else:
        qrf = (kwargs.get('h_factor', 1.0) * (zg / 20.0) +
               np.sqrt(yg ** 2 + xg ** 2) *
               np.tan(nb * bsp * np.pi / 180.0) + 500.0).flatten()

    # flattened query points

    ask = np.column_stack((xg.flatten(), yg.flatten(), zg.flatten()))

    # fill values

    badval = -9999.0

    # mask and flatten the data

    masked_data = np.ma.masked_array(
        np.column_stack([data[key][where_the_data_is_good]
                        for key in parms]),
        np.column_stack([data[key][where_the_data_is_good]
                        for key in parms]) == badval)

    # populate the ball tree

    mapping_obj = ballsy.BallsyMapper(np.column_stack((
        x[where_the_data_is_good], y[where_the_data_is_good],
        z[where_the_data_is_good])), masked_data, debug=False)

    # house keeping

    del data

    # query the tree and get the flattened interpolation
    #break it into parts for memory management
    asplit = np.split(ask, nx)
    qsplit = np.split(qrf, nx)
    interpols = []
    for i in range(nx):
        print(i)
        interpols.append(mapping_obj(asplit[i], qsplit[i], debug=True,
                         func='Barnes'))
    interpol = np.concatenate(interpols)
    grids = {}

    #reshape and store the grids in a dictionary

    for i in range(len(parms)):
        grids.update({parms[i]: interpol[:, i].reshape((nz, ny, nx))})
    grids.update({'ROI': qrf.reshape((nz, ny, nx))})
    return (xr, yr, zr), (nx, ny, nz), grids


def ncvar_to_field(ncvar):
    outdict = {'data': ncvar[:]}
    outdict.update(dict([(key, getattr(ncvar, key)) for key in
                         ncvar.ncattrs()]))
    return outdict


class pyGrid:

    def __init__(self, *args, **kwargs):
        if len(args) == 0:

            #initialize an empty pyGrid object
            self.fields = {}
            self.metadata = {}
            self.axes = {}
        elif 'variables' in dir(args[0]):
            if 'x_disp' in args[0].variables.keys():
                self.PyGridCF_to_PyGrid(args[0])
            elif 'TITLE' in dir(args[0]):
                self.WRFGridCF_to_PyGrid(args[0])
        elif 'count' in dir(args[0]):

            # a tuple of radar objects
            # grid the data

            (xr, yr, zr), (nx, ny, nz), grids = grid2(args[0], **kwargs)

            #create the fields

            self.fields = {}

            #move metadata from the radar to the grid

            for fld in grids.keys():
                if fld != 'ROI':
                    self.fields.update({fld: {'data': grids[fld]}})
                    for meta in args[0][0].fields[fld].keys():
                        if meta != 'data':
                            self.fields[fld].update(
                                {meta: args[0][0].fields[fld][meta]})

            self.fields.update({
                'ROI': {'data': grids['ROI'],
                'standard_name': 'radius_of_influence',
                'long_name': 'Radius of influence for mapping',
                'units': 'm',
                'least_significant_digit': 1,
                'valid_min': 0.,
                'valid_max': 100000.,
                '_FillValue': 9999.}})

            #create some axes

            #get location of the first radar...
            #in the absence of an origin kwarg we will default
            #to that

            first_radar_location = args[0][0].location
            lat, lon, alt = kwargs.get(
                'origin', (first_radar_location['latitude']['data'],
                           first_radar_location['longitude']['data'],
                           first_radar_location['altitude']['data']))

            x_array = np.linspace(xr[0], xr[1], nx)
            y_array = np.linspace(yr[0], yr[1], ny)
            z_array = np.linspace(zr[0], zr[1], nz)

            time = {
                'data': args[0][0].time['data'][0],
                'units': args[0][0].time['units'],
                'calendar': args[0][0].time['calendar'],
                'standard_name': args[0][0].time['standard_name'],
                'long_name': 'time in seconds of volume start'}

            time_start = {
                'data': args[0][0].time['data'][0],
                'units': args[0][0].time['units'],
                'calendar': args[0][0].time['calendar'],
                'standard_name': args[0][0].time['standard_name'],
                'long_name': 'time in seconds of volume start'}

            time_end = {
                'data': args[0][0].time['data'][-1],
                'units': args[0][0].time['units'],
                'calendar': args[0][0].time['calendar'],
                'standard_name': args[0][0].time['standard_name'],
                'long_name': 'time in seconds of volume end'}

            xaxis = {
                'data':  x_array,
                'long_name': 'x-coordinate in Cartesian system',
                'axis': 'X',
                'units': 'm'}

            yaxis = {
                'data': y_array,
                'long_name': 'y-coordinate in Cartesian system',
                'axis': 'Y',
                'units': 'm'}

            zaxis = {
                'data': z_array,
                'long_name': 'z-coordinate in Cartesian system',
                'axis': 'Z',
                'units': 'm',
                'positive': 'up'}

            altorigin = {
                'data': alt,
                'long_name': 'Altitude at grid origin',
                'units': 'm'}

            latorigin = {
                'data': lat,
                'long_name': 'latitude at grid origin',
                'units': 'degrees_north'}

            lonorigin = {
                'data': lon,
                'long_name': 'longitude at grid origin',
                'units': 'degrees_east'}

            self.axes = {'time': time, 'time_start': time_start,
                         'time_end': time_end,
                         'z_disp': zaxis, 'y_disp': yaxis, 'x_disp': xaxis,
                         'alt': altorigin, 'lat': latorigin, 'lon': lonorigin}

            self.metadata = args[0][0].metadata
            runtime = dict([(key, getattr(dt.datetime.now(), key)) for key in
                           ['year', 'month', 'day', 'hour', 'minute',
                            'second']])
            runtime.update({'strmon': dt.datetime.now().strftime('%b')})
            runtime.update({'user': getpass.getuser(),
                            'machine': socket.gethostname()})
            history_text = ("Gridded by user %(user)s on %(machine)s at"
                            "%(day)d-%(strmon)s-%(year)d,%(hour)d:%(minute)"
                            "02d:%(second)02d using grids2") % runtime
            if 'history' in self.metadata.keys():
                self.metadata['history'] = (self.metadata['history'] + '\n' +
                                            history_text)
            else:
                self.metadata['history'] = history_text
        else:
            print("foo")
            #TBI from various grid sources, WRF etc..

    def PyGridCF_to_PyGrid(self, netcdfobj):
            #netcdf file of grids
            #lets assume it is nicely formatted

            #grab the variable names

            all_variables=netcdfobj.variables.keys()
            fields=[]

            #anything that has more than 2 axes is a field

            for var in all_variables:
                if len(netcdfobj.variables[var].shape) > 1:
                    fields.append(var)
            self.fields={}

            for field in fields:
                self.fields.update({field:ncvar_to_field(netcdfobj.variables[field])})

            #now for axes, anything that is one of shape tuples
            self.axes={}
            axes=[]
            for var in all_variables:
                if len(netcdfobj.variables[var].shape) == 1:
                    axes.append(var)
            for axis in axes:
                self.axes.update({axis:ncvar_to_field(netcdfobj.variables[axis])})

    def WRFGridCF_to_PyGrid(self, netcdfobj):
            #netcdf file WRF of grids
            #lets assume it is nicely formatted

            #grab the variable names

            all_variables=netcdfobj.variables.keys()
            fields=[]

            #anything that has more than 2 axes is a field

            for var in all_variables:
                if len(netcdfobj.variables[var].shape) > 1:
                    fields.append(var)
            self.fields={}

            for field in fields:
                self.fields.update({field:ncvar_to_field(netcdfobj.variables[field])})

            #now for axes, anything that is one of shape tuples
            self.axes={}
            axes=[]
            for var in all_variables:
                if len(netcdfobj.variables[var].shape) == 1:
                    axes.append(var)
            for axis in axes:
                self.axes.update({axis:ncvar_to_field(netcdfobj.variables[axis])})


