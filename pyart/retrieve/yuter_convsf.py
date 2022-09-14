import numpy as np
import scipy.ndimage

def radialDistanceMask_array(mask_array, minradiuskm, maxradiuskm, x_pixsize,
                             y_pixsize, centerx, centery, circular=True):
    """
    Computes a radial distance mask, everything with distance between minradiuskm and maxradiuskm is assigned 1, everything else is assigned 0. This version can handle rectangular arrays and pixels as well as square ones.


    Parameters
    ----------
    mask_array : array
        Array to mask
    minradiuskm, maxradiuskm : float
        The minimum and maximum radius of the non-maked region in kilometers.
    x_pixsize, y_pixsize : float
        The pixel size in the x- and y-dimension in kilometers, respectively
    centerx, centery : int
        The center pixel in the x- and y-dimension, respectively

    Returns
    -------
    mask_array : array
        Rectangular array masked by a radial distance.
    """

    xsize, ysize = mask_array.shape

    for j in np.arange(0, ysize, 1):
        for i in np.arange(0, xsize, 1):
            # compute range to pixel
            if circular:
                x_range_sq = ((centerx - i) * x_pixsize) ** 2
                y_range_sq = ((centery - j) * y_pixsize) ** 2
                range = np.sqrt(x_range_sq + y_range_sq)
            # if circular is False, use square mask
            else:
                x_range = abs(int(np.floor(centerx - i) * x_pixsize))
                y_range = abs(int(np.floor(centery - j) * y_pixsize))

                if x_range > y_range:
                    range = x_range
                else:
                    range = y_range
            # if range is within min and max, set to True
            if (range <= maxradiuskm) and (range >= minradiuskm):
                mask_array[j, i] = 1
            else:
                mask_array[j, i] = 0

    return mask_array


def backgroundIntensity(refl, index_x, index_y, bkgDiameter, mask_array, dBaveraging):
    """
    For a pixel in the array determine the average intensity of the surrounding pixels in the window. The window is passed in as mask_array

    Parameters
    ----------
    refl : array
        Reflectivity array to compute average
    index_x, index_y : int
        x- and y-dimension index, respectively
    bkgDiameter : int
        diameter to compute average in pixels
    mask_array : array
        Array of radial points to use for average
    dBaveraging : bool
        If True, converts dBZ to linear Z before averaging

    Returns
    -------
    mean : float
        The average value for the given index
    """

    running_total = 0
    mean = np.nan
    pixelradius = 0
    numpixels = 0
    minptsAver = 1
    aval = 0
    centerval = 0
    maski = 0
    maskj = 0

    xsize, ysize = refl.shape
    # only compute background value if there is data at the point
    centerval = refl[index_y, index_x]
    if np.ma.is_masked(centerval):
        return mean

    pixelradius = int(np.floor(bkgDiameter / 2))

    # note the window will only be centered on point if windim is odd
    if pixelradius == (bkgDiameter / 2.0):
        print('Warning: windim = {0} is EVEN, background window will not be centered on point\n'.format(bkgDiameter))

    if (index_x >= 0) and (index_x < xsize) and (index_y >= 0) and (index_y < ysize):

        # loop through point in background window
        for j in np.arange(index_y - pixelradius, index_y + pixelradius + 1, 1):
            maski = 0
            for i in np.arange(index_x - pixelradius, index_x + pixelradius + 1, 1):
                # check that point is within bounds and activated in mask
                if (i >= 0) and (i < xsize) and (j >= 0) and (j < ysize) and mask_array[maskj, maski] == 1:
                    aval = refl[j, i]

                    if not np.ma.is_masked(aval):
                        numpixels += 1
                        # converting to linear Z
                        if dBaveraging:
                            running_total += 10 ** (aval / 10)
                        else:
                            running_total += aval
                maski += 1
            maskj += 1

        if numpixels >= minptsAver:
            mean = running_total / numpixels
            # converting to dBZ
            if dBaveraging:
                mean = 10 * np.log10(mean)

    return mean

def backgroundIntensity_array(refl, mask_array, dBaveraging):
    """
    For a pixel in the array determine the average intensity of the surrounding pixels in the window.
    The window is passed in as mask_array

    Parameters
    ----------
    refl : array
        Reflectivity array to compute average
    mask_array : array
        Array of radial points to use for average
    dBaveraging : bool
        If True, converts dBZ to linear Z before averaging

    Returns
    -------
    refl_bkg : array
        Array of average values
    """

    # if dBaverage is true, convert reflectivity to linear Z
    if dBaveraging:
        refl = 10 ** (refl / 10)

    # calculate background reflectivity with circular footprint
    refl_bkg = scipy.ndimage.generic_filter(refl.filled(np.nan), function=np.nanmean, mode='constant',
                                            footprint=mask_array.astype(bool), cval=np.nan)
    # mask where original reflectivity is invalid
    refl_bkg = np.ma.masked_where(refl.mask, refl_bkg)

    # if dBaveraging is true, convert background reflectivity to dBZ
    if dBaveraging:
        refl_bkg = 10 * np.log10(refl_bkg)
        # mask where original reflectivity is invalid
        refl_bkg = np.ma.masked_where(refl.mask, refl_bkg)

    return refl_bkg

def convcore_cos_scheme(zeVal, ze_bkg, minZeDiff, convThresB, alwaysConvThres, CS_CORE):
    # otherconvthres = absconvthres
    # intense = truncZconvthres
    """
    Cosine scheme for determining is convective core

    Parameters
    ----------
    zeVal : float
        Reflectivity value of point
    ze_bkg : float
        Reflectivity value of background value
    minZediff : float
        Minimum difference between zeVal and ze_bkg needed for convective classification
    convThresB : float
        Convective threshold used in the cosine function
    alwaysConvThres : float
        All values above this threshold considered to be convective

    Returns
    -------
    is_core : bool
        Boolean if point is convective (1) or not (0)
    """

    # initialize to not a convective core
    is_core = 0

    # if zval is greater than or equal to intense, set to core
    if zeVal >= alwaysConvThres:
        is_core = CS_CORE
    elif (not np.ma.is_masked(zeVal)) and (not np.ma.is_masked(ze_bkg)):
        # if background is less than zero, set difference to min difference
        if ze_bkg < 0:
            zeDiff = minZeDiff
        # else, use function from Yuter et al. (1997)
        else:
            zeDiff = minZeDiff * np.cos((np.pi * ze_bkg) / (2 * convThresB))
            # if difference is less than zero, set to zero
            if zeDiff < 0:
                zeDiff = 0
        # if value minus background is greater than or equal to difference, set to core
        if (zeVal - ze_bkg) >= zeDiff:
            is_core = CS_CORE

    return is_core

def convcore_cos_scheme_array(refl, ze_bkg, minZeDiff, convThresB, alwaysConvThres, CS_CORE):
    # otherconvthres = absconvthres
    # intense = truncZconvthres
    """
    Cosine scheme for determining is convective core

    Parameters
    ----------
    zeVal : float
        Reflectivity value of point
    ze_bkg : float
        Reflectivity value of background value
    minZediff : float
        Minimum difference between zeVal and ze_bkg needed for convective classification
    convThresB : float
        Convective threshold used in the cosine function
    alwaysConvThres : float
        All values above this threshold considered to be convective

    Returns
    -------
    is_core : bool
        Boolean if point is convective (1) or not (0)
    """

    # initialize entire array to not a convective core
    conv_core_array = np.zeros_like(refl)

    # calculate zeDiff for entire array
    zeDiff = minZeDiff * np.cos((np.pi * ze_bkg) / (2 * convThresB))
    zeDiff[zeDiff < 0] = 0 # where difference less than zero, set to zero
    zeDiff[ze_bkg < 0] = minZeDiff # where background less than zero, set to min diff

    # set values
    conv_core_array[refl >= alwaysConvThres] = CS_CORE # where Z is greater than alwaysConvThres, set to core
    conv_core_array[(refl - ze_bkg) >= zeDiff] = CS_CORE # where difference exceeeds minimum, set to core

    return conv_core_array

def convcore_scaled_array(refl, ze_bkg, minZeFactor, alwaysConvThres, CS_CORE, addition=False):
    # otherconvthres = absconvthres
    # intense = truncZconvthres
    """
    Cosine scheme for determining is convective core

    Parameters
    ----------
    zeVal : float
        Reflectivity value of point
    ze_bkg : float
        Reflectivity value of background value
    minZediff : float
        Minimum difference between zeVal and ze_bkg needed for convective classification
    convThresB : float
        Convective threshold used in the cosine function
    alwaysConvThres : float
        All values above this threshold considered to be convective

    Returns
    -------
    is_core : bool
        Boolean if point is convective (1) or not (0)
    """

    # initialize entire array to not a convective core
    conv_core_array = np.zeros_like(refl)

    # calculate zeDiff for entire array
    if addition:
        zeDiff = minZeFactor + ze_bkg
    else:
        zeDiff = minZeFactor * ze_bkg
    zeDiff[zeDiff < 0] = 0 # where difference less than zero, set to zero
    zeDiff[ze_bkg < 0] = 0 # where background less than zero, set to zero

    # set values
    conv_core_array[refl >= alwaysConvThres] = CS_CORE # where Z is greater than alwaysConvThres, set to core
    conv_core_array[refl >= zeDiff] = CS_CORE # where difference exceeeds minimum, set to core

    return conv_core_array

def classify_conv_strat_array(refl, conv_strat_array, conv_core_array,
                              NOSFCECHO, CONV, SF, WEAKECHO, CS_CORE, MINDBZUSE, WEAKECHOTHRES):
    """
    Does and initial convective stratiform classification

    Parameters
    ----------
    refl : array
        Array of reflectivity values
    conv_strat_array : array
        Array with convective stratiform classifications
    conv_core_array : array
        Array with convective cores
    NOSFCECHO : int
        Value to assign points classified as no surface echo
    CONV : int
        Value to assign points classified as convective
    SF : int
        Value to assign points classified as stratiform
    WEAKECHO : int
        Value to assign points classfied as weak echo
    CS_CORE : int
        Value assigned to convective cores in conv_core_array
    MINDBZUSE : float
        Minimum dBZ value to consider in classification, all values below this will be set to NOSFCECHO
    WEAKECHOTHRES : float
        dBZ threshold for weak echo classification, all values below this will be set to WEAKECHO

    Returns
    -------
    mean : conv_strat_array
        conv_strat_array with initial classifications
    """

    # assuming order so that each point is only assigned one time, no overlapping assignment
    # initially, assign every point to stratiform
    conv_strat_array[:] = SF
    # where reflectivity is masked, set to no surface echo
    conv_strat_array[refl.mask] = NOSFCECHO
    # assign convective cores to CONV
    conv_strat_array[conv_core_array==CS_CORE] = CONV
    # assign reflectivity less than weakechothres to weak echo
    conv_strat_array[refl < WEAKECHOTHRES] = WEAKECHO
    # assign reflectivity less than minimum to no surface echo
    conv_strat_array[refl < MINDBZUSE] = NOSFCECHO

    return conv_strat_array

def init_conv_radius_mask(maxConvDiameter, radius_km, xspacing, yspacing, centerConvMask_x):
    """
    Does and initial convective stratiform classification

    Parameters
    ----------
    maxConvDiameter : int
        maximum convective diameter in kilometers
    radius_km : int
        convective radius in kilometers
    xpacing, yspacing : float
        x- and y-dimension pixel size in meters, respectively
    centerConvMask_x : int
        index of center point

    Returns
    -------
    mean : conv_mask_array
        array masked based on distance of convective diameter
    """

    conv_mask_array = np.zeros((maxConvDiameter, maxConvDiameter))
    conv_mask_array = radialDistanceMask_array(conv_mask_array, 0, radius_km, xspacing, yspacing,
                                               centerConvMask_x, centerConvMask_x, True)

    return conv_mask_array

def assignConvRadiuskm_array(ze_bkg, dBZformaxconvradius, maxConvRadius=5):
    # alternative version for assigning convective radii
    # returns array the same size as ze_bkg with values for convective radii
    """
    Assigns the convective radius in kilometers based on the background reflectivity

    Parameters
    ----------
    ze_bkg : array
        array of background reflectivity values
    dBZformaxconvradius : float
        reflectivity value for maximum convective radius (5 km)
    maxConvRadius : float, optional
        maximum convective radius in kilometers

    Returns
    -------
    convRadiuskm : array
        array of convective radii based on background values and dBZ for max. conv radius
    """

    convRadiuskm = np.ones_like(ze_bkg)

    convRadiuskm[ze_bkg >= (dBZformaxconvradius - 15)] = maxConvRadius - 3
    convRadiuskm[ze_bkg >= (dBZformaxconvradius - 10)] = maxConvRadius - 2
    convRadiuskm[ze_bkg >= (dBZformaxconvradius - 5)] = maxConvRadius - 1
    convRadiuskm[ze_bkg >= dBZformaxconvradius] = maxConvRadius

    return convRadiuskm


def incorporateConvRadius(conv_strat_array, index_x, index_y, maxConvDiameter, conv_mask_array, NOSFCECHO, CONV):
    """
    Assigns the area within the convective radius around the convective core as convective
    function similar to backgroundIntensity except assigns value to mask array points in conv_strat_array
    only passing in legit convcore points to less error checking
    should only be called if convradius > 1

    Parameters
    ----------
    conv_strat_array : array
        Array with convective stratiform classifications
    index_x, index_y : int
        x- and y-dimension index, respectfully
    maxConvDiameter : int
        diameter to assign convective
    conv_mask_array : array
        Masked array of conv_strat_array

    Returns
    -------
    mean : conv_strat_array
        conv_strat_array with convective radii applied
    """

    xsize, ysize = conv_strat_array.shape

    maski = 0
    maskj = 0

    # note the convradius is always odd as it includes the conv core pixel at the center
    pixelradius = int(np.floor(maxConvDiameter / 2))

    # check limits of requested point are within array bounds
    if (index_x >= 0) and (index_x < xsize) and (index_y >= 0) and (index_y < ysize):

        # loop through the points in the square window in full size map
        for j in np.arange(index_y - pixelradius, index_y + pixelradius + 1, 1):
            maski = 0
            for i in np.arange(index_x - pixelradius, index_x + pixelradius + 1, 1):
                # check that point is within bounds, activated in mask, and has echo
                if (i >= 0) and (i < xsize) and (j >= 0) and (j < ysize) and \
                    conv_mask_array[maskj, maski] == 1 and conv_strat_array[j, i] != NOSFCECHO:
                    # assign pixel as convective
                    conv_strat_array[j,i] = CONV
                maski += 1
            maskj += 1

    return conv_strat_array


# def assignConvRadiuskm(ze_bkg, dBZformaxconvradius, maxconvradius=5):
#     # version for single background points
#     # returns single convective radii for given background value
#     if ze_bkg >= dBZformaxconvradius:
#         convRadiuskm = maxconvradius
#     elif ze_bkg >= (dBZformaxconvradius - 5):
#         convRadiuskm = maxconvradius - 1
#     elif ze_bkg >= (dBZformaxconvradius - 10):
#         convRadiuskm = maxconvradius - 2
#     elif ze_bkg >= (dBZformaxconvradius - 15):
#         convRadiuskm = maxconvradius - 3
#     else:
#         convRadiuskm = maxconvradius - 4
#
#     return convRadiuskm
