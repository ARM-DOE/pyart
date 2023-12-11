"""
Classification of Precipitation Echoes in Radar Data. 

Created on Thu Oct 12 23:12:19 2017
@author: Bhupendra Raut
@modifed: 02/13/2020
@references: 10.1109/TGRS.2020.2965649

.. autosummary::
    getWTClass
    label_classes
    reflectivity_to_rainrate
    calc_scale_break
    sum_conv_wavelets
    atwt2d
"""

import numpy as np
from numpy import log, floor
import sys


def get_reclass(grid,
                refl_field,
                level, 
                zr_a,
                zr_b,
                conv_wt_threshold, 
                tran_wt_threshold, 
                conv_scale_km, 
                min_dbz_threshold,
                conv_dbz_threshold,
                conv_core_threshold):
    """
    Compute ATWT described as Raut et al (2008) and classify radar echoes 
    using scheme of Raut et al (2020).
    First, convert dBZ to rain rates using standard Z-R relationship or user given coefiecient. This is to
    transform the normally distributed dBZ to gamma-like distribution, enhancing the structure of the field.

    Parameters:
    ===========
    dbz_data: ndarray
        3D array containing radar data. Last dimension should be levels.
    res_km: float
        Resolution of the radar data in km
    conv_scale_km: float
        Approximate scale break (in km) between convective and stratiform scales

    Returns:
    ========
    wt_class: ndarray
        Precipitation type classification: 0. N/A 1. stratiform/non-convective, 
        2. convective cores and 3. moderate+transitional (mix) convective 
        regions. 
    """

    # Extract grid data, save mask and get the resolution
    try:
        dbz_data = grid.fields[refl_field]['data'][level, :, :]
    except:
        dbz_data = grid.fields[refl_field]['data'][:, :]

    radar_mask = np.ma.getmask(dbz_data)

    # Warning: dx and dy are considred to be same (res_km).
    res_km = (grid.x["data"][1] - grid.x["data"][0])/1000

    try:
        dbz_data = dbz_data.filled(0)  # In case it's a masked array.
    except Exception:
        pass

    # save the mask for missing data.
    dbz_data[np.isnan(dbz_data)] = 0
    dbz_data_t = reflectivity_to_rainrate(dbz_data, acoeff=zr_a, bcoeff=zr_b)  # transform the dbz data
    
    # get scale break in pixels
    scale_break = calc_scale_break(res_km, conv_scale_km)
    wt_sum = sum_conv_wavelets(dbz_data_t, scale_break)

    wt_class = label_classes(wt_sum, 
                             dbz_data, 
                             conv_wt_threshold, 
                             tran_wt_threshold, 
                             min_dbz_threshold,
                             conv_dbz_threshold,
                             conv_core_threshold)

    wt_class_ma = np.ma.masked_where(radar_mask, wt_class) # add mask back
    wt_class_ma = wt_class_ma.squeeze()
 

    return wt_class_ma




def label_classes(wt_sum, 
                  dbz_data,
                  conv_wt_threshold, 
                  tran_wt_threshold, 
                  min_dbz_threshold,
                  conv_dbz_threshold,
                  conv_core_threshold):
    """ 
    Labels classes using given thresholds:
    - 0. no precipitation, marked by the given min_dbz threshold.
    - 1. stratiform,
    - 2. transitional/intermediate regions,
    - 3. intense/heavy convective.
    
    Following hard coded values are optimised and validated using C-band radars
    over Darwin, Australia (2.5 km grid spacing) and tested for Solapur, India (1km grid spacing) [Raut et al. 2020]. 
    conv_wt_threshold = 5  # WT value more than this is strong convection
    tran_wt_threshold = 2  # WT value for moderate convection
    min_dbz_threshold = 10  # pixels below this value are not classified. 
    conv_dbz_threshold = 30  # pixel below this value are not convective. This works for most cases.
git push -u 
    Parameters:
    ===========
    wt_sum: ndarray    
        Integrated wavelet transform
    vol_data: ndarray
        Array, vector or matrix of data

    Returns:
    ========
    wt_class: ndarray
        Precipitation type classification.
    """
 

    # I first used negative numbers to annotate the categories. Then multiply it by -1.
    wt_class = np.where((wt_sum >= tran_wt_threshold) & (dbz_data >= conv_core_threshold), -3, 0)
    wt_class = np.where((wt_sum >= conv_wt_threshold) & (dbz_data >= conv_dbz_threshold), -3, 0)
    wt_class = np.where((wt_sum < conv_wt_threshold) & (wt_sum >= tran_wt_threshold)
                        & (dbz_data >= conv_dbz_threshold), -2, wt_class)
    wt_class = np.where((wt_class == 0) & (dbz_data >= min_dbz_threshold), -1, wt_class)

    wt_class = -1 * wt_class
    wt_class = np.where((wt_class == 0), np.nan, wt_class)

    return wt_class.astype(np.int32)



def reflectivity_to_rainrate(dbz, acoeff, bcoeff):
    """
    Uses standard values for ZRA=200 and ZRB=1.6.

    Parameters:
    ===========
    dbz: ndarray
        Array, vector or matrix of reflectivity in dBZ.
    acoeff: float
        Z = a*R^b a coefficient.
    bcoeff: float
        Z = a*R^b b coefficient.

    Returns:
    ========
    rr: ndarray
        Rain rate in (mm/h)
    """
    rr = ((10.0 ** (dbz / 10.0)) / acoeff) ** (1.0 / bcoeff)
    return rr



def calc_scale_break(res_km, conv_scale_km):
    """
    Compute scale break for convection and stratiform regions. WT will be
    computed upto this scale and features will be designated as convection.

    Parameters:
    ===========
    res_km: float
        resolution of the image.
    conv_scale_km: float
        expected size of spatial variations due to convection.

    Returns:
    ========
    dyadic: int
        integer scale break in pixels.
    """
    scale_break = log((conv_scale_km / res_km)) / log(2) + 1
    return int(round(scale_break))



def sum_conv_wavelets(vol_data, conv_scale):
    """
    Returns sum of WT upto given scale. Works with both 2d scans and
    volumetric data.

    Parameters:
    ===========
    vol_data: ndarray
        Array, vector or matrix of data.
    conv_scale: float
        Expected size of spatial variations due to convection.

    Returns:
    ========
    wt_sum: ndarray
        Integrated wavelet transform.
    """
    dims = vol_data.shape

    # if data is 2d
    # if len(dims) == 2 or any(dims==1):
    wt, bg = atwt2d(vol_data, max_scale=conv_scale)
    wt_sum = np.sum(wt, axis=(0))
    """ else:  # else for volume data
        num_levels = min(dims) # too many assumptions here for height of the volume.
        wt_sum = np.zeros(dims)

        for lev in range(num_levels):
            if vol_data[:, :, lev].max < 1:
                next()  # this needs reviewing
            wt = atwt2d(vol_data[lev, :, :], max_scale=conv_scale)

            # sum all the WT scales.
            wt_sum[lev, :, :] = np.sum(wt, axis=(0)) """

    # Only positive WT corresponds to convection in radar
    # wt_sum[wt_sum<0] = 0 
    return wt_sum


def atwt2d(data2d, max_scale=-1):
    """
    Computes a trous wavelet transform (ATWT). Computes ATWT of the 2d array
    up to max_scale. If max_scale is outside the boundaries, number of scales
    will be reduced.

    Data is mirrored at the boundaries. 'Negative WT are removed. Not tested
    for non-square data.

    @authors: Bhupendra A. Raut and Dileep M. Puranik
    @references: Press et al. (1992) Numerical Recipes in C.

    Parameters:
    ===========
    data2d: ndarray
        2D image as array or matrix.
    max_scale:
        Computes wavelets up to max_scale. Leave blank for maximum possible
        scales.

    Returns:
    ===========
    tuple of ndarray
        ATWT of input image and the final smoothed image or background image.
    """
    
    if not isinstance(data2d, np.ndarray):
        sys.exit("the input is not a numpy array")

    data2d = data2d.squeeze()
    
    dims = data2d.shape
    min_dims = np.min(dims)
    max_possible_scales = int(floor(log(min_dims) / log(2)))

    if max_scale < 0 or max_possible_scales <= max_scale:
        max_scale = max_possible_scales - 1 

    ny = dims[0]
    nx = dims[1]    

    # For saving wt components
    wt = np.zeros((max_scale, ny, nx))

    temp1 = np.zeros(dims)
    temp2 = np.zeros(dims)

    sf = (0.0625, 0.25, 0.375)  # scaling function

    # start wavelet loop
    for scale in range(1, max_scale + 1):
        # print(scale)
        x1 = 2 ** (scale - 1)
        x2 = 2 * x1

        # Row-wise smoothing
        for i in range(0, nx):
            # find the indices for prev and next points on the line
            prev2 = abs(i - x2)
            prev1 = abs(i - x1)
            next1 = i + x1
            next2 = i + x2

            # If these indices are outside the image, "mirror" them
            # Sometime this causes issues at higher scales.
            if next1 > nx - 1:
                next1 = 2 * (nx - 1) - next1

            if next2 > nx - 1:
                next2 = 2 * (nx - 1) - next2

            if prev1 < 0 or prev2 < 0:
                prev1 = next1
                prev2 = next2

            for j in range(0, ny):
                left2 = data2d[j, prev2]
                left1 = data2d[j, prev1]
                right1 = data2d[j, next1]
                right2 = data2d[j, next2]
                temp1[j, i] = (sf[0] * (left2 + right2)
                               + sf[1] * (left1 + right1)
                               + sf[2] * data2d[j, i])

        # Column-wise smoothing
        for i in range(0, ny):
            prev2 = abs(i - x2)
            prev1 = abs(i - x1)
            next1 = i + x1
            next2 = i + x2

            # If these indices are outside the image use next values
            if next1 > ny - 1:
                next1 = 2 * (ny - 1) - next1
            if next2 > ny - 1:
                next2 = 2 * (ny - 1) - next2
            if prev1 < 0 or prev2 < 0:
                prev1 = next1
                prev2 = next2

            for j in range(0, nx):
                top2 = temp1[prev2, j]
                top1 = temp1[prev1, j]
                bottom1 = temp1[next1, j]
                bottom2 = temp1[next2, j]
                temp2[i, j] = (sf[0] * (top2 + bottom2)
                               + sf[1] * (top1 + bottom1)
                               + sf[2] * temp1[i, j])

        wt[scale - 1, :, :] = data2d - temp2
        data2d[:] = temp2


    return wt, data2d