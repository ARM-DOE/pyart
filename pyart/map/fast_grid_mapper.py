"""
pyart.io.fast_grid_mapper
===============

Highly optimized map_to_grid

.. autosummary::
    :toctree: generated/

    fast_map_to_grid


"""


import numpy as np
from ..core.grid import Grid
from ..graph.common import corner_to_point
from ..config import get_fillvalue
from .cython_fast_grid_mapper import from_cartesian_get_value,cython_fast_map


def fast_map_to_grid(radars, grid_shape, grid_limits, grid_origin=None,
                fields=None, refl_filter_flag=True, refl_field=None,
                max_refl=None, interpolation='nearest', max_range=None,
                round_resolution=(0.1,1,100)):
    """
    Map one or more radars, with same latlon, to a Cartesian grid.

    Generate a Cartesian grid of points for the requested fields from the collected 
    points from one or more radars. The field value for a grid point is found calculating 
    its corresponding elevation,azimuth and range,  and then rounded to a regular grid 
    that indicates the coresponding bin. Collected points are filtered according to a 
    number of criteria so that undesired points are not included in the interpolation.

    Parameters
    ----------
    radars : tuple of Radar objects.
        Radar objects which will be mapped to the Cartesian grid. However all radars must
        have the same latitude and longitude of the first one. They also need to be PPI scans, 
        have full sweeps. Partial or fail sweep will create undefined behaviour. 
    grid_shape : 3-tuple of floats
        Number of points in the grid (z, y, x).
    grid_limits : 3-tuple of 2-tuples or numpy 1d-array
        if 2-tuples: Minimum and maximum grid location (inclusive) in meters for the
        z, y, x coordinates.
        if numpy 1d-array: grid locations in meters for the z, y, x coordinates, must
        agree in size with grid_shape
    grid_origin : (float, float) or None
        Latitude and longitude of grid origin.  None sets the origin
        to the location of the first radar.
    fields : list or None
        List of fields within the radar objects which will be mapped to
        the cartesian grid. None, the default, will map the fields which are
        present in all the radar objects.
    refl_filter_flag : bool
        True to filter the collected points based on the reflectivity field.
        False to perform no filtering.  Gates where the reflectivity field,
        specified by the `refl_field` parameter, is not-finited, masked or
        has a value above the `max_refl` parameter are excluded from the
        grid interpolation. 
    refl_field : str
        Name of the field which will be used to filter the collected points.
        A value of None will use the default field name as defined in the
        Py-ART configuration file.
    max_refl : float
        Maximum allowable reflectivity.  Points in the `refl_field` which are
        above is value are not included in the interpolation. None will
        include skip this filtering. 
    interpolation: "nearest" or "linear"
        Interpolation to be user, between 'nearest' (no interpolation) and 
        'linear' (linear in the elevation coordinate). 'linear' takes fast double
        the time.
    max_range: float
        maximum range, bins with greater radius will be excluded. If `None`, it will
        be set to be the diagonal of the grid region. 
    
    Other Parameters
    ----------------
    round_resolution: 3-tuple of float
        Resolution for which data will be rounded in (elev,azi,range) grid, 
        values given in (angles,angles,meters).
    
    Returns
    -------
    grids : dict
        Dictionary of mapped fields.  The keys of the dictionary are given by
        parameter fields.  Each elements is a `grid_size` float64 array
        containing the interpolated grid for that field.
    
    See Also
    --------
    map_to_grid : Similar function with more complet interpolation options and
    support to multiple radars
    
    """
    # unpack the grid parameters
    nz, ny, nx = grid_shape
    zr, yr, xr = grid_limits
    #let us allow irregular grids, this may be useful in doing PPI's

    if len(zr)==nz:
        z=zr
    elif len(zr)==2:
        z_start, z_stop = zr
        z=np.linspace(z_start, z_stop, nz)
    else:
        raise ValueError('Z variable does not agree in size')
    if len(yr)==ny:
        y=yr
    elif len(yr)==2:
        y_start, y_stop = yr
        y=np.linspace(y_start, y_stop, ny)
    else:
        raise ValueError('Y variable does not agree in size')
    if len(xr)==nx:
        x=xr
    elif len(xr)==2:
        x_start, x_stop = xr
        x=np.linspace(x_start, x_stop, nx)
    else:
        raise ValueError('X variable does not agree in size')
    
    #for fast access, convert interpolation array into int 0=='nearest', 1=='linear'
    if interpolation=='nearest':
        interp_opt=0
    elif interpolation=='linear':
        interp_opt=1
    else:
        raise ValueError("Known interpolation option: %s"%interpolation)
    
    badval = get_fillvalue()
# find the grid origin if not given
    if grid_origin is None:
        lat = float(radars[0].latitude['data'])
        lon = float(radars[0].longitude['data'])
        grid_origin = (lat, lon)
        offset=(radars[0].altitude['data'],0,0) #all radar have the same location
    else:
        radar_lat = float(radar.latitude['data'])
        radar_lon = float(radar.longitude['data'])
        x_disp, y_disp = corner_to_point(grid_origin, (radar_lat, radar_lon))
        offset=(radars[0].altitude['data'], y_disp, x_disp) 
        
        
# find max range if not given
    if max_range==None:
        max_range=np.sqrt((grid_limits[0][1]-grid_limits[0][0])**2+(grid_limits[1][1]-grid_limits[1][0])**2+(grid_limits[2][1]-grid_limits[2][0])**2)+1

# fields which should be mapped, None for fields which are in all radars
    if fields is None:
        fields = set(radars[0].fields.keys())
        for radar in radars[1:]:
            fields = fields.intersection(radar.fields.keys())
        fields = list(fields)

    # put reflectivity as last field, this will be used in cython side for refl_filter_flag
    if refl_field==None:
        refl_field='reflectivity'
    if refl_field in fields:
        ref=fields.pop(fields.index(refl_field))
        fields.append(ref)
    elif refl_filter_flag:
        raise ValueError('None reflectivity field found')
    nfields = len(fields)
    
    if max_refl==None:
        max_refl = 120 # cython will need a number to compare
    
    # unpack round_resolution and force float
    elev_res = float(round_resolution[0])
    azi_res = float(round_resolution[1])
    r_res = float(round_resolution[2])
    
    # map elevation to sweep/radar
    sweeps = [] #[(elevation,sweep,radar)]
    for iradar, radar in enumerate(radars):
        elev = [(e,i,iradar) for i,e in enumerate(radar.fixed_angle["data"])]
        sweeps = sweeps + elev
    sweeps = sorted(sweeps, key=lambda sweep: sweep[0])
    sweeps_elev = np.array([item[0] for item in sweeps])
    sweeps_index_down = np.empty(int(180/elev_res+1), dtype=np.int) #-90 to +90 with 'elev_res' elevation precision
    sweeps_index_down.fill(-1)
    sweeps_index_up = np.empty(int(180/elev_res+1), dtype=np.int) 
    sweeps_index_up.fill(-1)
    for i in range(int(180/elev_res+1)):
        nominal_elev = (i-90./elev_res)*elev_res
        if interp_opt==0:
            sweeps_index_down[i]=np.abs(sweeps_elev-nominal_elev).argmin() #it says down, but it is nearest
        else:
            index = len([0 for sweep in sweeps if sweep[0]<= nominal_elev])-1 
            if index>=0 and sweeps[index][0]+5. > nominal_elev:
                sweeps_index_down[i] = index #map elevation down
            if index< len(sweeps)-1:
                sweeps_index_up[i] = index+1 #map elevation up
    
    # map sweep,azimuth to ray/radar
    ray_index = np.empty((len(sweeps),int(360/azi_res+1)), dtype=np.int) #0 to 360 with 'azi_res' azimuth precision
    ray_index.fill(-1)
    for isweep ,sweep in enumerate(sweeps):
        radar = radars[sweep[2]]
        ini   = radar.sweep_start_ray_index['data'][sweep[1]]
        end   = radar.sweep_end_ray_index['data'][sweep[1]]
        rays = sorted([(azi,ray+ini) for ray,azi in enumerate(radar.azimuth['data'][ini:end+1])], key=lambda ray: ray[0])
        index = 0
        for azi in range(int(360/azi_res+1)):
            nominal_azi = azi*azi_res
            if index==-1:
                index = 0
            for i in range(len(rays)-index):
                if rays[i+index][0]>nominal_azi: break #map azimuth down
            index = i+index-1
            ray_index[isweep,azi]=rays[index][1]
    # map radar,range to bin
    gate_index = np.empty((len(radars),int(max_range/r_res)), dtype=np.int) #0 to max_range with 'r_res' range precision 
    gate_index.fill(-1)
    for iradar, radar in enumerate(radars):
        index = 0
        for i in range(int(max_range/r_res)):
            if index==-1:
                index = 0
            nominal_range = i*r_res
            for k in range(len(radar.range["data"])-index):
                if radar.range["data"][k+index]>nominal_range: break  #map range down,
            if index+k+1==len(radar.range["data"]):
                index=k+index
                gate_index[iradar,i] = index
                break
            else:
                index = k+index-1
                gate_index[iradar,i] = index
    
    grid_data_simple = cython_fast_map(z,y,x, offset[0],offset[1] , offset[2], max_range,badval, interp_opt, sweeps,sweeps_index_down, sweeps_index_up, ray_index,gate_index, radars, fields,refl_filter_flag,max_refl, elev_res, azi_res, r_res)

    # create and return the grid dictionary
    grid_data = np.ma.empty((nz, ny, nx, nfields), dtype=np.float64)
    grid_data.data[:] = grid_data_simple[:]
    grid_data = np.ma.masked_where( grid_data.data == badval, grid_data)
    
    # create and return the grid dictionary
    grids = dict([(f, grid_data[..., i]) for i, f in enumerate(fields)])
    return grids
    
#XXX The rest of this file is deprecated
#XXX 
#XXX 
#XXX 
#XXX 
#XXX 
#XXX 

# for every point in grid found the coresponding bin
# could do this with numpy array but would need to much memory, so make one index at time and parallelize it by hand, also cython if not enough.
    #shared memory, must be simple array, not masked
    grid_data_shared_base = multiprocessing.Array(ctypes.c_double, nz*ny*nx*nfields)
    grid_data_shared      = np.ctypeslib.as_array(grid_data_shared_base.get_obj())
    grid_data_shared      = grid_data_shared.reshape(nz, ny, nx, nfields)
    
    # list of processes
    processes = []
    # twice as many threads as cpu's, this shall increase usuability of the actual cpu's
    threads = 2*multiprocessing.cpu_count()
    
    array_indexes = np.indices((nz,ny,nx)).T.reshape(-1,len((nz,ny,nx)))
    
    
    for thread in range(threads):
        processes.append( multiprocessing.Process( target=thread_loop, args=(thread, threads, array_indexes, grid_data_shared ,x,y,z,offset,max_range,badval,radars,fields,nfields,sweeps,sweeps_index,ray_index,gate_index) ) )
    #rolling processes
    for process in processes:
       process.start()
    #waiting to finish
    for process in processes:
       process.join() 
    #from_cartesian_get_value(0, 1, array_indexes, grid_data_shared ,x,y,z,offset,max_range,badval,radars,fields,nfields,sweeps,sweeps_index,ray_index,gate_index) 
    #masking values
    grid_data = np.ma.empty((nz, ny, nx, nfields), dtype=np.float64)
    grid_data.data[:] = grid_data_shared[:]
    grid_data = np.ma.masked_where( grid_data.data == badval, grid_data)

    # create and return the grid dictionary
    grids = dict([(f, grid_data[..., i]) for i, f in enumerate(fields)])
    return grids



def thread_loop(thread, threads, array_indexes, grid_data_shared ,x,y,z,offset,max_range,badval,radars,fields,nfields,sweeps,sweeps_index,ray_index,gate_index):

    total_points = array_indexes.shape[0]
    
    for i in range(thread, total_points, threads):#threads processes intercalated indexes
        (iz,iy,ix) = array_indexes[i]
        value = from_cartesian_get_value(z[iz],y[iy],x[ix], offset[0],offset[1] , offset[2], max_range,badval, sweeps,sweeps_index, ray_index,gate_index, radars, fields, nfields)
        grid_data_shared[iz, iy, ix] = value


# XXX move this to another module
# XXX not being used, cython version used instead
def cart_to_radar_coords(z,y,x,z_offset,y_offset,x_offset):
    """
    
    ADD FUNCTION DESCRIPTION HERE
    copy radar_coords_to_cart
    
    verify math
    
    """
    R = 6371.0 * 1000.0 * 4.0 / 3.0     # effective radius of earth in meters.
    
    h = z-z_offset
    s = np.sqrt((x-x_offset)*(x-x_offset)+(y-y_offset)*(y-y_offset)) #radius to radar
    azi = np.arctan2((x-x_offset),(y-y_offset))*180/np.pi # azi in degree
    azi = np.remainder(azi,360.) # 0 to 360
    sigma = s / R
    r0 = R * np.tan(sigma)
    betha = np.pi/2 + sigma
    h1 = h - (np.sqrt(r0*r0+R*R)-R)
    rng = np.sqrt(r0*r0+h1*h1-2*r0*h1*np.cos(betha)) 
    elev = np.arcsin(h1*np.sin(betha)/rng)*180/np.pi #in degrees
    return (rng,azi,elev)

