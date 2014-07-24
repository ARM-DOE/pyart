"""
pyart.graph.common
==================

Common graphing routines.

.. autosummary::
    :toctree: generated/

    dms_to_d
    radar_coords_to_cart
    corner_to_point
    ax_radius

"""

import numpy as np

PI = 3.141592653589793

from ..io.common import dms_to_d, radar_coords_to_cart

def centers_to_edges_2d(x):
    """ 
    Create a (N+1, M+1) array of edge locations from a
    (N, M) array of grid center locations.
    
    In the interior, the edge positions set to the midpoints
    of the values in x. For the outermost edges, half the 
    closest dx is assumed to apply. This matters for polar
    meshes, where one edge of the grid becomes a point at the
    polar coordinate origin; dx/2 is a half-hearted way of
    trying to prevent negative ranges.
    
    Useful when plotting with pcolor, which requires
    X, Y of shape (N+1) and grid center values of shape (N).
    Otherwise, pcolor silently discards the last row and column
    of grid center values.
    
    Parameters
    ----------
    x : array, shape (N,M)
        Locations of the centers 
    
    Returns
    -------
    xedge : array, shape (N+1,M+1)
    
    """
    xedge = np.zeros((x.shape[0]+1,x.shape[1]+1))
    # interior is a simple average of four adjacent centers
    xedge[1:-1,1:-1] = (x[:-1,:-1] + x[:-1,1:] + x[1:,:-1] + x[1:,1:])/4.0
    
    #         /\
    #        /\/\
    #       / /\ \
    #      /\/  \/\
    #     / /\  /\ \
    #    /\/  \/  \/\
    #   / /\  /\  /\ \
    #  /\/  \/  \/  \/\
    #4 \/\  /\  /\  /\/ 4
    # 3 \ \/  \/  \/ / 3 
    #    \/\  /\  /\/
    #   2 \ \/  \/ / 2  
    #      \/\  /\/
    #     1 \ \/ / 1
    #        \/\/
    #       0 \/ 0 = center ID of 0th dimension
    #
    
    # calculate the deltas along each edge, excluding corners
    xedge[1:-1,0] = xedge[1:-1, 1] - (xedge[1:-1, 2] - xedge[1:-1, 1])/2.0
    xedge[1:-1,-1]= xedge[1:-1,-2] - (xedge[1:-1,-3] - xedge[1:-1,-2])/2.0
    xedge[0,1:-1] = xedge[1,1:-1]  - (xedge[2,1:-1]  - xedge[1,1:-1])/2.0 
    xedge[-1,1:-1]= xedge[-2,1:-1] - (xedge[-3,1:-1] - xedge[-2,1:-1])/2.0
    
    # now do the corners
    xedge[0,0]  = xedge[1, 1] - (xedge[2, 2] - xedge[1, 1])/2.0
    xedge[0,-1] = xedge[1,-2] - (xedge[2,-3] - xedge[1,-2])/2.0
    xedge[-1,0] = xedge[-2,1] - (xedge[-3,2] - xedge[-2,1])/2.0 
    xedge[-1,-1]= xedge[-2,-2]- (xedge[-3,-3]- xedge[-2,-2])/2.0
    
    return xedge


def centers_to_edges(x):
    """ 
    Create a length N+1 vector of edge locations from a
    length N vector of center locations.
    
    In the interior, the edge positions set to the midpoints
    of the values in x. For the outermost edges, the closest 
    dx is assumed to apply. This matters for polar
    meshes, where one edge of the grid becomes a point at the
    polar coordinate origin; dx/2 is a half-hearted way of
    trying to prevent negative ranges.
    
    Useful when plotting with pcolor, which requires
    X, Y of shape (N+1) and grid center values of shape (N).
    Otherwise, pcolor silently discards the last row and column
    of grid center values.
    
    Parameters
    ----------
    x : array, shape (N,)
        Locations of the centers 
    
    Returns
    -------
    xedge : array, shape (N+1,)
    
    """
    xedge=np.zeros(x.shape[0]+1)
    xedge[1:-1] = (x[:-1] + x[1:])/2.0
    dx = np.mean(np.abs(xedge[2:-1] - xedge[1:-2]))
    xedge[0] = xedge[1] - dx/2.0
    xedge[-1] = xedge[-2] + dx/2.0
    return xedge        


def corner_to_point(corner, point):
    """
    Return the x, y distances in meters from a corner to a point.

    Assumes a spherical earth model.

    Parameters
    ----------
    corner : (float, float)
        Latitude and longitude in degrees of the corner.
    point : (float, float)
        Latitude and longitude in degrees of the point.

    Returns
    -------
    x, y : floats
        Distances from the corner to the point in meters.

    """
    Re = 6371.0 * 1000.0
    Rc = ax_radius(point[0], units='degrees')
    #print Rc/Re
    y = ((point[0] - corner[0]) / 360.0) * PI * 2.0 * Re
    x = ((point[1] - corner[1]) / 360.0) * PI * 2.0 * Rc
    return x, y


def ax_radius(lat, units='radians'):
    """
    Return the radius of a constant latitude circle for a given latitude.

    Parameters
    ----------
    lat : float
        Latitude at which to calculate constant latitude circle (parallel)
        radius.
    units : 'radians' or 'degrees'
        Units of lat, either 'radians' or 'degrees'.

    Returns
    -------
    R : float
        Radius in meters of a constant latitude circle (parallel).

    """
    Re = 6371.0 * 1000.0
    if units == 'degrees':
        const = PI / 180.0
    else:
        const = 1.0
    R = Re * np.sin(PI / 2.0 - abs(lat * const))
    return R
