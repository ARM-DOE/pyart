"""
A general gridded data object for radar data

"""

import numpy as np
import pyart.map.ballsy as ballsy
from pyart.graph import radar_display
import getpass
import socket
import datetime as dt


# now to grid
def dms_to_d(dms):
    """ Convert degrees minutes seconds tuple to decimal degrees """
    return dms[0] + (dms[1] + dms[2] / 60.0) / 60.0


def grid2(radars, **kwargs):
    """
    map tuples of Py-Radar objects to a cartesian grid..

    Usage: (xr,yr,zr), (nx,ny,nz), grids2 = grid2(
        (radar1, radar2, radarn,  nx=81, ny=81, nz=69,
         zr=(0.,17000.), yr=(-20000., 20000.), xr=(-30000., 20000.),
         origin=(lat, lon, height_m))

    """
    nx, ny, nz = kwargs.get('nxyz', (81, 81, 69))
    xr, yr, zr = kwargs.get(
        'xyzr', ((-30000., 20000), (-20000., 20000.), (0., 17000.)))
    toa = kwargs.get('toa', 17000.0)
    # need to just set to the radar..
    first_radar_location=radars[0].location
    cf_lat, cf_lon, cf_alt=kwargs.get(
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
    data = dict([(parms[i], np.array([])) for i in range(len(parms))])
    
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
        
        xdash, ydash, zdash = radar_display.radar_coords_to_cart(rg, azg, eleg)
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
    
    where_the_data_is_good = np.where(np.logical_and(
        np.isfinite(data[parms[0]]), data[parms[0]] < 100.0))[0]
    
    # Create a meshgrid(cube) to allow calculation of radii of influence
    
    zg, yg, xg = np.mgrid[zr[0]:zr[1]:nz * 1j, yr[0]:yr[1]:ny * 1j,
                          xr[0]:xr[1]:nx * 1j]
    
    # Virtual beam width and beam spacing
    
    nb = 1.5
    bsp = 1.0
    
    # Query radius of influence, flattened
    
    qrf = (zg / 20.0 + np.sqrt(yg ** 2 + xg ** 2) *
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
        z[where_the_data_is_good])), masked_data, debug=True)
    
    # house keeping
    
    del data
    
    # query the tree and get the flattened interpolation
    
    interpol = mapping_obj(ask, qrf, debug=True, func='Barnes')
    grids = {}
    
    #reshape and store the grids in a dictionary
    
    for i in range(len(parms)):
        grids.update({parms[i]: interpol[:, i].reshape((nz, ny, nx))})
    return (xr, yr, zr), (nx, ny, nz), grids


class pyGrid:

    def __init__(self, *args, **kwargs):
        if len(args) == 0:

            #initialize an empty pyGrid object
            self.fields = {}
            self.metadata = {}
            self.axes = {}
        elif 'count' in dir(args[0]):
            
            # a tuple of radar objects
            # grid the data
            
            (xr, yr, zr), (nx, ny, nz), grids = grid2(args[0], **kwargs)
            
            #create the fields
            
            self.fields = {}
            
            #move metadata from the radar to the grid
            
            for fld in grids.keys():
                self.fields.update({fld: {'data': grids[fld]}})
                for meta in args[0][0].fields[fld].keys():
                    if meta != 'data':
                        self.fields[fld].update(
                            {meta: args[0][0].fields[fld][meta]})
            
            #create some axes
            
            #get location of the first radar... 
            #in the absence of an origin kwarg we will default
            #to that
            
            first_radar_location=args[0][0].location
            lat, lon, alt=kwargs.get(
                'origin', (first_radar_location['latitude']['data'],
                           first_radar_location['longitude']['data'],
                           first_radar_location['altitude']['data']))
            
            x_array = np.linspace(xr[0], xr[1], nx)
            y_array = np.linspace(yr[0], yr[1], ny)
            z_array = np.linspace(zr[0], zr[1], nz)
            
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
                'data':alt,
                'long_name': 'Altitude at grid origin',
                'units': 'm'}
            
            latorigin = {
                'data':lat,
                'long_name': 'latitude at grid origin',
                'units': 'degrees_north'}
            
            lonorigin = {
                'data':lon,
                'long_name': 'longitude at grid origin',
                'units': 'degrees_east'}

            
            self.axes = {'time_start': time_start, 'time_end': time_end,
                         'z_disp': zaxis, 'y_disp': yaxis, 'x_disp': xaxis,
                         'alt':altorigin, 'lat':latorigin, 'lon':lonorigin}
            
            self.metadata = args[0][0].metadata
            runtime = dict([(key, getattr(dt.datetime.now(), key)) for key in
                           ['year', 'month', 'day', 'hour', 'minute',
                            'second']])
            runtime.update({'strmon': dt.datetime.now().strftime('%b')})
            runtime.update({'user': getpass.getuser(),
                            'machine': socket.gethostname()})
            history_text = "Gridded by user %(user)s on %(machine)s at %(day)d-%(strmon)s-%(year)d,%(hour)d:%(minute)02d:%(second)02d using grids2" % runtime
            if 'history' in self.metadata.keys():
                self.metadata['history'] = (self.metadata['history'] + '\n' +
                                            history_text)
            else:
                self.metadata['history'] = history_text
        else:
            print("foo")
            #TBI from various grid sources, WRF etc..
