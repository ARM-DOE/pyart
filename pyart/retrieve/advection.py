"""
pyart.retrieve.advection
=========================


.. autosummary::
    :toctree: generated/

    grid_displacement_pc
    grid_shift
    add_grids

"""

import numpy as np
from netCDF4 import datetime
import scipy
import copy
from ..config import get_fillvalue


#Based off work by Christoph Gohlke <http://www.lfd.uci.edu/~gohlke/>
def grid_displacement_pc(grid1, grid2, var, level, return_value = 'pixels'):
    """
    Calculate the grid displacement using Phase correlation.
    See:
    http://en.wikipedia.org/wiki/Phase_correlation
    implimentation inspired by Christoph Gohlke:
    http://www.lfd.uci.edu/~gohlke/code/imreg.py.html

    Notes: Make sure valid min is set to something sensible.
    missing values are set to valid_min. Works best on fields which have
    background areas set to a constant area, eg rain rate, correlation etc..

     Parameters
    ----------
    grid1,grid2 : Grids
        Py-ART Grid objects seperated in time and be square in x/y.
    var : string
        Variable to calculate advection from. var must be in both grid1
        and grid2.fields.keys()
    level : integer
        The level of the 3D grid to use in the calculation
    return_value : string
        'pixels', 'distance' or 'velocity'. Distance in pixels (default)
        or meters or velocity vector in m/s

    Returns
    -------
    displacement : two-tuple
         integers if pixels, otherwise floats. Result of the calculation


    """
    #create copies of the data

    ige1 = copy.deepcopy(grid2.fields[var]['data'][level,:,:])
    ige2 = copy.deepcopy(grid1.fields[var]['data'][level,:,:])

    #set anything below valid min to valid min, this helps nuke _fill_values

    ige1[np.where(ige1 < grid1.fields[var]['valid_min'])] = \
                                                  grid1.fields[var]['valid_min']
    ige2[np.where(ige2 < grid2.fields[var]['valid_min'])] = \
                                                  grid2.fields[var]['valid_min']

    # discrete fast fourier transformation and complex conjugation of image 2

    image1FFT = np.fft.fft2(ige1)
    image2FFT = np.conjugate(np.fft.fft2(ige2))

    # inverse fourier transformation of product -> equal to cross correlation

    imageCCor = np.real( np.fft.ifft2( (image1FFT*image2FFT) ) )

    # Shift the zero-frequency component to the center of the spectrum

    imageCCorShift = np.fft.fftshift(imageCCor)
    # determine the distance of the maximum from the center

    row, col = ige1.shape

    #find the peak in the correlation

    yShift, xShift = np.unravel_index( np.argmax(imageCCorShift), (row,col) )

    yShift -= int(row/2)
    xShift -= int(col/2)

    if return_value == 'pixels':
        tbr = (xShift, yShift)
    elif return_value == 'distance':
        dx = grid1.axes['x_disp']['data'][1] - grid1.axes['x_disp']['data'][0]
        dy = grid1.axes['y_disp']['data'][1] - grid1.axes['y_disp']['data'][0]
        x_movement = xShift * dx
        y_movement = yShift * dy
        tbr = (x_movement, y_movement)
    elif return_value == 'velocity':
        dx = grid1.axes['x_disp']['data'][1] - grid1.axes['x_disp']['data'][0]
        dy = grid1.axes['y_disp']['data'][1] - grid1.axes['y_disp']['data'][0]
        x_movement = xShift * dx
        y_movement = yShift * dy
        t1 = netCDF4.num2date(grid1.axes['time']['data'][0],
                grid1.axes['time']['units'])
        t2 = netCDF4.num2date(grid2.axes['time']['data'][0],
                grid2.axes['time']['units'])
        dt = (t2 - t1).seconds
        u = x_movement/dt
        v = y_movement/dt
        tbr = (u, v)
    else:
        tbr = (xShift, yShift)
    return tbr

def grid_shift(grid, advection, trim_edges = 0., field_list=None):
    """
    Use scipy.ndimage to shift a grid by a certain number of pixels
     Parameters
    ----------
    grid: Grid
        Py-ART Grid object.
    advection : two-tuple of floats
        Pixels to shift the image by.
    trim_edges: integer
        edges to cut off the grid and axes, both x and y. Defaults to zero.

    Returns
    -------
    return_grid : Grid
         grid with fields shifted and, if requested, subset.
    """

    new_grid = copy.deepcopy(grid)

    #grab the x and y axis and trim

    x_g = grid.axes['x_disp']['data'].copy()
    y_g = grid.axes['y_disp']['data'].copy()
    if trim_edges !=0.:
        new_grid.axes['x_disp']['data'] = x_g[trim_edges:-trim_edges]
        new_grid.axes['y_disp']['data'] = y_g[trim_edges:-trim_edges]
    else:
        new_grid.axes['x_disp']['data'] = x_g
        new_grid.axes['y_disp']['data'] = y_g

    #now shift each field. Replacing masking is tricky..
    #either use valid min and max or use  user set
    if field_list == None:
        field_list =  grid.fields.keys()

    for field in field_list:
        image_data = grid.fields[field]['data'].copy()
        try:
            ndv = image_data.fill_value
            image_data = image_data.filled(np.nan)
        except:
            ndv = get_fillvalue()
        new_image = scipy.ndimage.interpolation.shift(image_data,
                [0,advection[0],advection[1]], prefilter=False)
        image_masked = np.ma.fix_invalid(new_image, copy=False)
        image_masked.set_fill_value(ndv)
        if trim_edges != 0.:
            new_grid.fields[field]['data'] = image_masked[:,
                    trim_edges:-trim_edges, trim_edges:-trim_edges]
        else:
            new_grid.fields[field]['data'] = image_masked
    return new_grid

def add_grids(grids, weights=None, fields=None):
    """
    Add a list of grids together. Note: must be same dimensionality


     Parameters
    ----------
    grids : list of Grids
        A python list objects to be aggregated.
    weights : list of floats
        The weights for the grids. Defaults to [1.0]*len(grids).
        Must be same length as grids
    fields : list of strings
        List of fields to be aggregated, defaults to all grids.

    Returns
    -------
    return_grid : Grid
         Aggregated Grid Resulting fields are SUM(grids[n] * weights[n])


    """

    #default is to aggregate all grids

    if fields == None:
        fields = grids[0].fields.keys()

    #default weights are 1.0, so a simple ad

    if weights == None:
        weights = [1.0]*len(grids)

    return_grid = copy.deepcopy(grids[0])
    return_grid.fields =  {k: return_grid.fields.get(k, None) for k in fields}

    for fld in fields:
        datas = [this_grid.fields[fld]['data'] for this_grid in grids]
        mean_grid = np.ma.average(datas, axis = 0, weights = weights)
        return_grid.fields[fld]['data'] = mean_grid

    return return_grid

