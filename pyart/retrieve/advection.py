"""
pyart.retrieve.advection
=========================


.. autosummary::
    :toctree: generated/

"""

import numpy as np
from netCDF4 import datetime
from skimage.transform import warp, AffineTransform

#Based off work by Christoph Gohlke <http://www.lfd.uci.edu/~gohlke/>
def grid_displacememt_pc(grid1, grid2, var, level, return_value = 'pixels'):
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



