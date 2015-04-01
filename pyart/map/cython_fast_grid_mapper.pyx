#!python
#cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True
from libc.math cimport *
cimport numpy as np
import numpy as np
cimport cython



def  cython_fast_map(z,y,x, double z_offset, double y_offset, double x_offset , double max_range,badval, sweeps,np.ndarray[np.int_t, ndim=1] sweeps_index, np.ndarray[np.int_t, ndim=2] ray_index, np.ndarray[np.int_t, ndim=2] gate_index, radars, fields,int refl_filter_flag, double max_refl):
    """
    cython_fast_map(z,y,x, z_offset, y_offset, x_offset , max_range,badval, sweeps, sweeps_index, ray_index, gate_index, radars, fields, refl_filter_flag, max_refl):
    
    cython part of function fast_map_to_grid
    
    """
    cdef size_t nz = len(z)
    cdef size_t ny = len(y)
    cdef size_t nx = len(x)
    cdef size_t nfields = len(fields)
    #convert x,y,z to memory view for better acess
    cdef double[:] cz = z
    cdef double[:] cy = y
    cdef double[:] cx = x
    cdef double rng,azi,elev
    cdef int temp_int
    cdef int i,j,k,l,
    cdef double[:,:,:,:] grid_data_simple = np.empty((nz,ny,nx,nfields), dtype='double')
    grid_data_simple[:,:,:,:] = badval
    cdef int isweep,iray,iradar,igate
    size = radars[0].fields['DBZH']['data'].data.shape

    # remove python overhead of getting the mask and array
    list_array=[range(nfields) for i in range(len(radars))] #list of lists
    list_mask=[range(nfields) for i in range(len(radars))] #list of lists
    
    for iradar in range(len(radars)):
        for l in range(nfields):
            list_array[iradar][l] = radars[iradar].fields[fields[l]]['data'].data
            list_mask[iradar][l] = radars[iradar].fields[fields[l]]['data'].mask

    for i in range(nz):
        for j in range(ny):
            for k in range(nx):
                cart_to_radar_coords (cz[i],cy[j],cx[k],z_offset, y_offset, x_offset, &rng, &azi, &elev)
                 
                if rng<0 or rng>=max_range:
#                    grid_data_simple[i,j,k,:] = badval
                    continue
                
                temp_int = <int> ceil((elev+90)*2)
                isweep = sweeps_index[temp_int]
                if isweep<0: #no elevation, this offen happen
#                    grid_data_simple[i,j,k,:] = badval
                    continue
                
                temp_int = <int> ceil(azi)
                iray = ray_index[isweep, temp_int]
                if iray<0: #no ray, this almost does't happen
#                    grid_data_simple[i,j,k,:] = badval
                    continue
                
                temp_int = <int> ceil(rng/100)
                iradar = sweeps[isweep][2]
                igate = gate_index[iradar, temp_int ]
                if igate<0:#no gate, this almost does't happen
#                    grid_data_simple[i,j,k,:] = badval
                    continue
                
                # we have data
                for l in range(nfields-1):#last field is reflectivity, process separated
                    if list_mask[iradar][l][iray,igate]:
                        grid_data_simple[i,j,k,l] = badval
                    elif refl_filter_flag:
                        if list_mask[iradar][nfields-1][iray,igate] or list_array[iradar][nfields-1][iray,igate]>max_refl:
                            grid_data_simple[i,j,k,l] = badval
                        else:
                            grid_data_simple[i,j,k,l] = list_array[iradar][l][iray,igate] # radars[iradar].fields[fields[l]]['data'].data[iray,igate]
                # new process reflectivity
                l=nfields-1
                if list_mask[iradar][l][iray,igate]:
                    grid_data_simple[i,j,k,l] = badval
                else:
                    grid_data_simple[i,j,k,l] = list_array[iradar][l][iray,igate] # radars[iradar].fields[fields[l]]['data'].data[iray,igate]
                    if grid_data_simple[i,j,k,l]>max_refl:
                        grid_data_simple[i,j,k,l] = badval
                        
    return np.asarray(grid_data_simple)


cdef void cart_to_radar_coords( double z, double y, double x, double z_offset, double y_offset,double x_offset,double *rng,double *azi, double *elev):
    """
    Calculate Radar coordinate from Cartesian coordinates

    Parameters
    ----------
    x, y, z : double
        Cartesian coordinates in meters.
    x_offset, y_offset, z_offset : double
        Posicion of the radar in Cartesian coordinates in meters.
    rng : array
        Distances to the center of the radar gates (bins) in kilometers.
    az : array
        Azimuth angle of the radar in degrees.
    ele : array
        Elevation angle of the radar in degrees.

    Returns
    -------
    rng : double
        Distances to the center of the radar gates (bins) in kilometers.
    azi : double
        Azimuth angle of the radar in degrees.
    elev : double
        Elevation angle of the radar in degrees.

    Notes
    -----
    The calculation for radar coordinate is derivide from the model of a 
    straight ray in a 4/3 Earth's radius model.

    .. math::

        azi = \\atan2(x,y)

        s = \\sqrt(x^2+y^2)

        r0 = R * \\tan(s/R)

        h1 = z - \\sqrt(r0^2+R^2)-R

        rng = \\sqrt(r0^2+h1^2-2*r0*h1*\\cos(pi/2. + s/R)) 

        elev = \\asin(h1*\\sin(pi/2. + s/R)/rng)

    Where s is the arc length, r0 is range of the 0° elevation ray at this 
    arc length, h1 is altitude in respect to the 0° elevation ray, rng
    is the distance from the radar to the center of the gate and R is the 
    effective radius of the earth, taken to be 4/3 of the mean radius of 
    the earth (6371 km).

    """
    cdef double R = 6371.0 * 1000.0 * 4.0 / 3.0     # effective radius of earth in meters.
    cdef double pi = 3.14159265358979323846
    cdef double h = z-z_offset
    cdef double s = sqrt((x-x_offset)*(x-x_offset)+(y-y_offset)*(y-y_offset)) #radius to radar
    azi[0] = atan2((x-x_offset),(y-y_offset))*180./pi # azi in degree
    if azi[0]<0:
        azi[0] = azi[0] + 360# 0 to 360
    cdef double sigma = s / R
    cdef double r0 = R * tan(sigma)
    cdef double betha = pi/2. + sigma
    cdef double h1 = h - (sqrt(r0*r0+R*R)-R)
    rng[0] = sqrt(r0*r0+h1*h1-2*r0*h1*cos(betha)) 
    elev[0] = asin(h1*sin(betha)/(rng[0]))*180./pi #in degrees


### this function is deprecated, use full cython version above: cython_fast_map 
def from_cartesian_get_value( double z, double y,double x, double z_offset, double y_offset, double x_offset , double max_range , double badval , sweeps, np.ndarray[np.int_t, ndim=1] sweeps_index, np.ndarray[np.int_t, ndim=2] ray_index, np.ndarray[np.int_t, ndim=2] gate_index, radars, fields, int nfields):
    cdef double rng,azi,elev
    cdef int temp_int
    cart_to_radar_coords (z,y,x,z_offset, y_offset, x_offset, &rng, &azi, &elev) 
    if rng<0 or rng>=max_range:
        return badval
    
    temp_int = <int> ceil((elev+90)*2)
    cdef int isweep = sweeps_index[temp_int]
    if isweep<0: #no elevation, this offen happen
        return badval

    temp_int = <int> ceil(azi)
    cdef int iray = ray_index[isweep, temp_int]
    if iray<0: #no ray, this almost does't happen
        return badval

    temp_int = <int> ceil(rng/100)
    cdef int iradar = sweeps[isweep][2]
    cdef int igate = gate_index[iradar, temp_int ]
    if igate<0:#no gate, this almost does't happen
        return badval
        
    # we have data
    value = np.empty((nfields), dtype=np.float64)
    for ifield in range(nfields):
        if radars[iradar].fields[fields[ifield]]['data'].mask[iray,igate]:
            value[ifield] = badval
        else:
            value[ifield] = radars[iradar].fields[fields[ifield]]['data'].data[iray,igate]
    return value


