"""
Advection calculations.

"""

import copy

import numpy as np
from scipy.ndimage import interpolation
from netCDF4 import num2date

from ..config import get_fillvalue


# Based off work by Christoph Gohlke <http://www.lfd.uci.edu/~gohlke/>


def grid_displacement_pc(grid1, grid2, field, level, return_value='pixels'):
    """
    Calculate the grid displacement using phase correlation.

    See:
    http://en.wikipedia.org/wiki/Phase_correlation

    Implementation inspired by Christoph Gohlke:
    http://www.lfd.uci.edu/~gohlke/code/imreg.py.html

    Note that the grid must have the same dimensions in x and y and assumed to
    have constant spacing in these dimensions.

    Parameters
    ----------
    grid1, grid2 : Grid
        Py-ART Grid objects separated in time and square in x/y.
    field : string
        Field to calculate advection from. Field must be in both grid1
        and grid2.
    level : integer
        The vertical (z) level of the grid to use in the calculation.
    return_value : str, optional
        'pixels', 'distance' or 'velocity'. Distance in pixels (default)
        or meters or velocity vector in m/s.

    Returns
    -------
    displacement : two-tuple
         Calculated displacement in units of y and x. Value returned in
         integers if pixels, otherwise floats.

    """
    # create copies of the data
    field_data1 = grid1.fields[field]['data'][level].copy()
    field_data2 = grid2.fields[field]['data'][level].copy()

    # replace fill values with valid_min or minimum value in array
    if 'valid_min' in grid1.fields[field]:
        min_value1 = grid1.fields[field]['valid_min']
    else:
        min_value1 = field_data1.min()
    field_data1 = np.ma.filled(field_data1, min_value1)

    if 'valid_min' in grid2.fields[field]:
        min_value2 = grid2.fields[field]['valid_min']
    else:
        min_value2 = field_data2.min()
    field_data2 = np.ma.filled(field_data2, min_value2)

    # discrete fast fourier transformation and complex conjugation of field 2
    image1fft = np.fft.fft2(field_data1)
    image2fft = np.conjugate(np.fft.fft2(field_data2))

    # inverse fourier transformation of product -> equal to cross correlation
    imageccor = np.real(np.fft.ifft2((image1fft*image2fft)))

    # shift the zero-frequency component to the center of the spectrum
    imageccorshift = np.fft.fftshift(imageccor)

    # determine the distance of the maximum from the center
    # find the peak in the correlation
    row, col = field_data1.shape
    yshift, xshift = np.unravel_index(np.argmax(imageccorshift), (row, col))
    yshift -= int(row/2)
    xshift -= int(col/2)

    dx = grid1.x['data'][1] - grid1.x['data'][0]
    dy = grid1.y['data'][1] - grid1.y['data'][0]
    x_movement = xshift * dx
    y_movement = yshift * dy

    if return_value == 'pixels':
        displacement = (yshift, xshift)
    elif return_value == 'distance':
        displacement = (y_movement, x_movement)
    elif return_value == 'velocity':
        t1 = num2date(grid1.time['data'][0], grid1.time['units'])
        t2 = num2date(grid2.time['data'][0], grid2.time['units'])
        dt = (t2 - t1).total_seconds()
        u = x_movement/dt
        v = y_movement/dt
        displacement = (v, u)
    else:
        displacement = (yshift, xshift)
    return displacement


def grid_shift(grid, advection, trim_edges=0, field_list=None):
    """
    Shift a grid by a certain number of pixels.

    Parameters
    ----------
    grid: Grid
        Py-ART Grid object.
    advection : two-tuple of floats
        Number of Pixels to shift the image by.
    trim_edges: integer, optional
        Edges to cut off the grid and axes, both x and y. Defaults to zero.
    field_list : list, optional
        List of fields to include in new grid. None, the default, includes all
        fields from the input grid.

    Returns
    -------
    shifted_grid : Grid
         Grid with fields shifted and, if requested, subset.

    """
    if trim_edges == 0:
        trim_slice = slice(None, None)
    else:
        trim_slice = slice(int(trim_edges), -int(trim_edges))

    shifted_grid = copy.deepcopy(grid)

    # grab the x and y axis and trim
    shifted_grid.x['data'] = grid.x['data'][trim_slice].copy()
    shifted_grid.y['data'] = grid.y['data'][trim_slice].copy()

    # shift each field.
    if field_list is None:
        field_list = grid.fields.keys()

    for field in field_list:

        # copy data and fill with nans
        data = grid.fields[field]['data'].copy()
        data = np.ma.filled(data, np.nan)

        # shift the data
        shifted_data = interpolation.shift(
            data, [0, advection[0], advection[1]], prefilter=False)

        # mask invalid, trim and place into grid
        shifted_data = np.ma.fix_invalid(
            shifted_data, copy=False, fill_value=get_fillvalue())
        shifted_data = shifted_data[:, trim_slice, trim_slice]
        shifted_grid.fields[field]['data'] = shifted_data

    return shifted_grid
