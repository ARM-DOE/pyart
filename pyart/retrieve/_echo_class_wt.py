"""
Classification of Precipitation Echoes in Radar Data.

Created on Thu Oct 12 23:12:19 2017
@author: Bhupendra Raut
@modifed: 11/19/2023
@references: 10.1109/TGRS.2020.2965649

.. autosummary::
    wavelet_reclass
    label_classes
    calc_scale_break
    atwt2d
"""

import numpy as np


def wavelet_reclass(
    grid,
    refl_field,
    level,
    zr_a,
    zr_b,
    core_wt_threshold,
    conv_wt_threshold,
    scale_break,
    min_reflectivity,
    conv_min_refl,
    conv_core_threshold,
):
    """
    Compute ATWT described as Raut et al (2008) and classify radar echoes using scheme of Raut et al (2020).
    First, convert dBZ to rain rates using standard Z-R relationship or user given coefficients. This is to
    transform the normally distributed dBZ to gamma-like distribution, enhancing the structure of the field.

    Parameters
    ----------
    dbz_data : ndarray
        2D array containing radar data. Last dimension should be levels.
    res_km : float
        Resolution of the radar data in km
    scale_break : int
        Calculated scale break between convective and stratiform scales. Dyadically spaced in grid pixels.

    Returns
    -------
    wt_class : ndarray
        Precipitation type classification: 0. N/A 1. stratiform/non-convective,
        2. convective cores and 3. moderate+transitional (mix) convective
        regions.
    """

    # Extract grid data, save mask and get the resolution in km
    try:
        dbz_data = grid.fields[refl_field]["data"][level, :, :]
    except:
        dbz_data = grid.fields[refl_field]["data"][:, :]

    # save the radar original mask for missing data.
    radar_mask = np.ma.getmask(dbz_data)

    wt_sum = conv_wavelet_sum(dbz_data, zr_a, zr_b, scale_break)

    wt_class = label_classes(
        wt_sum,
        dbz_data,
        core_wt_threshold,
        conv_wt_threshold,
        min_reflectivity,
        conv_min_refl,
        conv_core_threshold,
    )

    wt_class_ma = np.ma.masked_where(radar_mask, wt_class)  # add mask back
    wt_class_ma = wt_class_ma.squeeze()

    return wt_class_ma


def conv_wavelet_sum(dbz_data, zr_a, zr_b, scale_break):
    """
    Computes the sum of wavelet transform components for convective scales from dBZ data.

    Parameters
    ------------
    dbz_data : ndarray
        2D array containing radar dBZ data.
    zr_a, zr_b : float
        Coefficients for the Z-R relationship.
    res_km : float
        Resolution of the radar data in km.
    scale_break : int
        Calculated scale break (in pixels) between convective and stratiform scales

    Returns
    ---------
    wt_sum : ndarray
        Sum of convective scale wavelet transform components.
    """
    try:
        dbz_data = dbz_data.filled(0)
    except Exception:
        pass

    dbz_data[np.isnan(dbz_data)] = 0
    rr_data = ((10.0 ** (dbz_data / 10.0)) / zr_a) ** (1.0 / zr_b)

    wt, _ = atwt2d(rr_data, max_scale=scale_break)
    wt_sum = np.sum(wt, axis=(0))

    return wt_sum


def label_classes(
    wt_sum,
    dbz_data,
    core_wt_threshold,
    conv_wt_threshold,
    min_reflectivity,
    conv_min_refl,
    conv_core_threshold,
):
    """
    Labels classes using given thresholds:
        - 0: No precipitation or unclassified
        - 1: Stratiform/non-convective regions
        - 2: Transitional and mixed convective regions
        - 3: Convective cores

    Following hard coded values are optimized and validated using C-band radars
    over Darwin, Australia (2.5 km grid spacing) and tested for Solapur, India (1km grid spacing) [Raut et al. 2020].
    core_wt_threshold = 5  # WT value more than this is strong convection
    conv_wt_threshold = 2  # WT value for moderate convection
    min_reflectivity = 10  # pixels below this value are not classified.
     conv_min_refl = 30  # pixel below this value are not convective. This works for most cases.

    Parameters
    -----------
    wt_sum : ndarray
        Integrated wavelet transform
    vol_data : ndarray
        Array, vector or matrix of data

    Returns
    ---------
    wt_class : ndarray
        Precipitation type classification.
    """

    # I first used negative numbers to annotate the categories. Then multiply it by -1.
    wt_class = np.where(
        (wt_sum >= conv_wt_threshold) & (dbz_data >= conv_core_threshold), -3, 0
    )
    wt_class = np.where(
        (wt_sum >= core_wt_threshold) & (dbz_data >= conv_min_refl), -3, 0
    )
    wt_class = np.where(
        (wt_sum < core_wt_threshold)
        & (wt_sum >= conv_wt_threshold)
        & (dbz_data >= conv_min_refl),
        -2,
        wt_class,
    )
    wt_class = np.where((wt_class == 0) & (dbz_data >= min_reflectivity), -1, wt_class)

    wt_class = -1 * wt_class
    wt_class = np.where((wt_class == 0), np.nan, wt_class)

    return wt_class.astype(np.int32)


def calc_scale_break(res_meters, conv_scale_km):
    """
    Compute scale break for convection and stratiform regions. WT will be
    computed upto this scale and features will be designated as convection.

    Parameters
    -----------
    res_meters : float
        resolution of the image.
    conv_scale_km : float
        expected size of spatial variations due to convection.

    Returns
    --------
    dyadic scale break : int
        integer scale break in dyadic scale.
    """
    res_km = res_meters / 1000
    scale_break = np.log(conv_scale_km / res_km) / np.log(2) + 1

    return int(round(scale_break))


def atwt2d(data2d, max_scale=-1):
    """
    Computes a trous wavelet transform (ATWT). Computes ATWT of the 2D array
    up to max_scale. If max_scale is outside the boundaries, number of scales
    will be reduced.

    Data is mirrored at the boundaries. 'Negative WT are removed. Not tested
    for non-square data.

    @authors: Bhupendra A. Raut and Dileep M. Puranik
    @references: Press et al. (1992) Numerical Recipes in C.

    Parameters
    -----------
    data2d : ndarray
        2D image as array or matrix.
    max_scale :
        Computes wavelets up to max_scale. Leave blank for maximum possible
        scales.

    Returns
    ---------
    tuple of ndarray
        ATWT of input image and the final smoothed image or background image.
    """

    if not isinstance(data2d, np.ndarray):
        raise TypeError("The input data2d must be a numpy array.")

    data2d = data2d.squeeze()

    dims = data2d.shape
    min_dims = np.min(dims)
    max_possible_scales = int(np.floor(np.log(min_dims) / np.log(2)))

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
                temp1[j, i] = (
                    sf[0] * (left2 + right2)
                    + sf[1] * (left1 + right1)
                    + sf[2] * data2d[j, i]
                )

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
                temp2[i, j] = (
                    sf[0] * (top2 + bottom2)
                    + sf[1] * (top1 + bottom1)
                    + sf[2] * temp1[i, j]
                )

        wt[scale - 1, :, :] = data2d - temp2
        data2d[:] = temp2

    return wt, data2d
