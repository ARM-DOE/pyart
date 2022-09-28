import numpy as np
import sys
sys.path.append('C:\\Users\\lmtomkin\\Documents\\GitHub\\pyart_convsf\\pyart\\retrieve\\')

import scipy.ndimage

def _revised_conv_strat(refl, dx, dy, alwaysConvThres=42, bkgRad_km=11,
                        useCosine=True, maxDiff=8, zeroDiffCosVal=55,
                        weakEchoThres=5.0, mindBZused=5.0,
                        scalarDiff=1.5, addition=True,
                        dBaveraging=False, applyLgRadialMask=False,
                        lgRadialMask_minRadkm=0, lgRadialMask_maxRadkm=170,
                        dBZforMaxConvRad=30, maxConvRad_km=5.0,
                        incorpConvRad=True):

    """
    We perform the Steiner et al. (1995) algorithm for echo classification
    using only the reflectivity field in order to classify each grid point
    as either convective, stratiform or undefined. Grid points are
    classified as follows,

    0 = No Surface Echo/ Undefined
    1 = Stratiform
    2 = Convective
    3 = Weak Echo

    ref : array
        array of reflectivity values
    x, y : array
        x and y coordinates of reflectivity array, respectively
    dx, dy : float
        The x- and y-dimension resolutions in meters, respectively.
    alwaysConvThres : float, optional
        Threshold for points that are always convective. All values above the threshold are classifed as convective
    minZeDiff : float, optional
        Minimum difference between background average and reflectivity in order to be classified as convective.
        a value in Yuter et al. (2005)
    convThresB : float, optional
        Convective threshold used in cosine function for classifying convective vs. stratiform
        b value in Yuter et al. (2005)
    bkg_rad : float, optional
        Radius to compute background reflectivity in kilometers. Default is 11 km
    maxConvRad_km : float, optional
        Maximum radius around convective cores to classify as convective. Default is 5 km
    weakEchoThres : float, optional
        Threshold for determining weak echo. All values below this threshold will be considered weak echo
    minDBZused : float, optional
        Minimum dBZ value used for classification. All values below this threshold will be considered no surface echo
    applyRadialMask : bool, optional
        Flag to set a radial mask for algorithm
    dBZforMaxConvRadius : float, optional
        dBZ for maximum convective radius. Convective cores with values above this will have the maximum convective radius
    dBaveraging : bool, optional
        True if using dBZ values that need to be converted to linear Z before averaging. False for other types of values

    """

    if maxConvRad_km > 5:
        print("Max conv radius must be less than 5 km, exiting")
        # quit

    # Constants to fill arrays with
    CS_CORE = 3
    NOSFCECHO = 0
    WEAKECHO = 3
    SF = 1
    CONV = 2

    #%% Set up mask arrays
    # prepare for convective mask arrays
    # calculate maximum convective diameter from max. convective radius (input)
    maxConvDiameter = int(np.floor((maxConvRad_km / (dx / 1000)) * 2))
    # if diameter is even, make odd
    if maxConvDiameter % 2 == 0: maxConvDiameter = maxConvDiameter + 1
    # find center point
    centerConvMask_x = int(np.floor(maxConvDiameter / 2))

    # prepare background mask array for computing background average
    # calculate number of pixels for background array given requested background radius and dx
    bkgDiameter_pix = int(np.floor((bkgRad_km / (dx / 1000)) * 2))
    # set diameter to odd if even
    if bkgDiameter_pix % 2 == 0: bkgDiameter_pix = bkgDiameter_pix + 1
    # find center point
    bkg_center = int(np.floor(bkgDiameter_pix / 2))
    # create background array
    bkg_mask_array = np.ones((bkgDiameter_pix, bkgDiameter_pix), dtype=float)
    # mask outside circular region
    bkg_mask_array = radialDistanceMask(bkg_mask_array, minradiuskm=0, maxradiuskm=bkgRad_km, x_pixsize=dx / 1000,
                                        y_pixsize=dy / 1000, centerx=bkg_center, centery=bkg_center, circular=True)

    # Create large mask array for determining where to calculate convective stratiform
    # initialize array with 1 (calculate convective stratiform over entire array)
    mask_array = np.zeros(refl.shape, dtype=float)
    mask_array[:] = 1
    # if True, create radial mask
    if applyLgRadialMask:
        mask_array = radialDistanceMask(mask_array, lgRadialMask_minRadkm, lgRadialMask_maxRadkm, x_pixsize=dx / 1000,
                                        y_pixsize=dy / 1000, centerx=int(np.floor(refl.shape[0] / 2)),
                                        centery=int(np.floor(refl.shape[1] / 2)), circular=True)

    #%% Convective stratiform detection

    # Compute background radius
    refl_bkg = backgroundIntensity(refl, bkg_mask_array, dBaveraging)
    # mask background average
    refl_bkg = np.ma.masked_where(refl.mask, refl_bkg)

    # Get convective core array from cosine scheme, or scalar scheme
    if useCosine:
        conv_core_array = convcore_cos_scheme(refl, refl_bkg, maxDiff, zeroDiffCosVal, alwaysConvThres, CS_CORE)
    else:
        conv_core_array = convcore_scaled(refl, refl_bkg, scalarDiff, alwaysConvThres, CS_CORE, addition=addition)

    # count convective cores
    corecount = np.count_nonzero(conv_core_array)

    # Assign convective radii based on background reflectivity
    convRadiuskm = assignConvRadiuskm(refl_bkg, dBZformaxconvradius=dBZforMaxConvRad, maxConvRadius=maxConvRad_km)

    # Incorporate convective radius using binary dilation
    # Create empty array for assignment
    temp_assignment = np.zeros_like(conv_core_array)

    # Loop through radii
    for radius in np.arange(1, maxConvRad_km+1):
        # create mask array for radius incorporation
        conv_mask_array = init_conv_radius_mask(maxConvDiameter, radius, dx / 1000, dy / 1000, centerConvMask_x)
        # find location of radius
        temp = convRadiuskm == radius
        # get cores for given radius
        temp_core = np.ma.masked_where(~temp, conv_core_array)
        # dilate cores
        temp_dilated = scipy.ndimage.binary_dilation(temp_core.filled(0), conv_mask_array)
        # add to assignment array
        temp_assignment = temp_assignment + temp_dilated

    # add dilated cores to original array
    conv_core_copy = np.copy(conv_core_array)
    conv_core_copy[temp_assignment >= 1] = CS_CORE

    # Now do convective stratiform classification
    conv_strat_array = np.zeros_like(refl)
    conv_strat_array = classify_conv_strat_array(refl, conv_strat_array, conv_core_copy,
                                                 NOSFCECHO, CONV, SF, WEAKECHO, CS_CORE,
                                                 mindBZused, weakEchoThres)

    return refl_bkg, conv_core_array, conv_strat_array

# functions

def radialDistanceMask(mask_array, minradiuskm, maxradiuskm, x_pixsize,
                       y_pixsize, centerx, centery, circular=True):
    """
    Computes a radial distance mask, everything with distance between minradiuskm
    and maxradiuskm is assigned 1, everything else is assigned 0. This version can
    handle rectangular arrays and pixels as well as square ones.

    Parameters
    ----------
    mask_array : array
        Array to mask
    minradiuskm, maxradiuskm : float
        The minimum and maximum radius of the non-masked region in kilometers.
    x_pixsize, y_pixsize : float
        The pixel size in the x- and y-dimension in kilometers, respectively
    centerx, centery : int
        The center pixel in the x- and y-dimension, respectively
    circular : bool
        True returns circular mask

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

def backgroundIntensity(refl, bkg_mask_array, dBaveraging):
    """
    Calculate the background of the given refl array. The footprint used to
    calculate the average for each pixel is given by bkg_mask_array

    Parameters
    ----------
    refl : array
        Reflectivity array to compute average
    bkg_mask_array : array
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
                                            footprint=bkg_mask_array.astype(bool), cval=np.nan)
    # mask where original reflectivity is invalid
    refl_bkg = np.ma.masked_where(refl.mask, refl_bkg)

    # if dBaveraging is true, convert background reflectivity to dBZ
    if dBaveraging:
        refl_bkg = 10 * np.log10(refl_bkg)
        # mask where original reflectivity is invalid
        refl_bkg = np.ma.masked_where(refl.mask, refl_bkg)

    return refl_bkg


def convcore_cos_scheme(refl, refl_bkg, maxDiff, zeroDiffCosVal, alwaysConvThres, CS_CORE):
    """
    Function for assigning convective cores based on a cosine function

    Parameters
    ----------
    refl : array
        Reflectivity values
    refl_bkg : float
        Background average of reflectivity values
    maxDiff : float
        Maximum difference between refl and refl_bkg needed for convective classification
    zeroDiffCosVal : float
        Value where the cosine function returns a zero difference
    alwaysConvThres : float
        All values above this threshold considered to be convective
    CS_CORE : int
        Value assigned to convective pixels

    Returns
    -------
    conv_core_array : array
        Array of booleans if point is convective (1) or not (0)
    """

    # initialize entire array to not a convective core
    conv_core_array = np.zeros_like(refl)

    # calculate zeDiff for entire array
    zDiff = maxDiff * np.cos((np.pi * refl_bkg) / (2 * zeroDiffCosVal))
    zDiff[zDiff < 0] = 0 # where difference less than zero, set to zero
    zDiff[refl_bkg < 0] = maxDiff # where background less than zero, set to min diff

    # set values
    conv_core_array[refl >= alwaysConvThres] = CS_CORE # where Z is greater than alwaysConvThres, set to core
    conv_core_array[(refl - refl_bkg) >= zDiff] = CS_CORE # where difference exceeeds minimum, set to core

    return conv_core_array


def convcore_scaled(refl, refl_bkg, maxDiff, alwaysConvThres, CS_CORE, addition=False):

    """
    Function for assigning convective cores based on a scalar difference

    Parameters
    ----------
    refl : array
        Reflectivity values
    refl_bkg : float
        Background average of reflectivity values
    maxDiff : float
        Maximum difference between refl and refl_bkg needed for convective classification
    alwaysConvThres : float
        All values above this threshold considered to be convective
    CS_CORE : int
        Value assigned to convective pixels
    addition : bool
        Boolean to determine if scalar should be added (True) or multiplied (False)

    Returns
    -------
    conv_core_array : array
        Array of booleans if point is convective (1) or not (0)
    """

    # initialize entire array to not a convective core
    conv_core_array = np.zeros_like(refl)

    # calculate zDiff for entire array
    # if addition, add difference. Else, multiply difference
    if addition:
        zDiff = maxDiff + refl_bkg
    else:
        zDiff = maxDiff * refl_bkg

    zDiff[zDiff < 0] = 0 # where difference less than zero, set to zero
    zDiff[refl_bkg < 0] = 0 # where background less than zero, set to zero

    # set values
    conv_core_array[refl >= alwaysConvThres] = CS_CORE # where Z is greater than alwaysConvThres, set to core
    conv_core_array[refl >= zDiff] = CS_CORE # where difference exceeeds minimum, set to core

    return conv_core_array

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
   conv_mask_array : array
        array masked based on distance of convective diameter
    """

    conv_mask_array = np.zeros((maxConvDiameter, maxConvDiameter))
    conv_mask_array = radialDistanceMask(conv_mask_array, 0, radius_km, xspacing, yspacing,
                                         centerConvMask_x, centerConvMask_x, True)

    return conv_mask_array

def assignConvRadiuskm(refl_bkg, dBZformaxconvradius, maxConvRadius=5):

    """
    Assigns the convective radius in kilometers based on the background reflectivity

    Parameters
    ----------
    refl_bkg : array
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

    convRadiuskm = np.ones_like(refl_bkg)

    convRadiuskm[refl_bkg >= (dBZformaxconvradius - 15)] = maxConvRadius - 3
    convRadiuskm[refl_bkg >= (dBZformaxconvradius - 10)] = maxConvRadius - 2
    convRadiuskm[refl_bkg >= (dBZformaxconvradius - 5)] = maxConvRadius - 1
    convRadiuskm[refl_bkg >= dBZformaxconvradius] = maxConvRadius

    return convRadiuskm

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