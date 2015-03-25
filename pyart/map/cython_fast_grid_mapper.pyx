#!python
#cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True
from libc.math cimport *
cimport numpy as np
import numpy as np
cimport cython


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






cdef void cart_to_radar_coords( double z, double y, double x, double z_offset, double y_offset,double x_offset,double *rng,double *azi, double *elev):
    """
    
    ADD FUNCTION DESCRIPTION HERE
    copy radar_coords_to_cart
    
    verify math
    
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


