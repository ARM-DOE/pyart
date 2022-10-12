import numpy as np
import scipy.ndimage


def _revised_conv_strat(refl, dx, dy, always_core_thres=42, bkg_rad_km=11,
                        use_cosine=True, max_diff=8, zero_diff_cos_val=55,
                        scalar_diff=1.5, use_addition=True, calc_thres=0.75,
                        weak_echo_thres=5.0, min_dBZ_used=5.0,
                        dB_averaging=False, apply_lg_rad_mask=False,
                        lg_rad_mask_min_rad_km=0, lg_rad_mask_max_rad_km=170,
                        val_for_max_conv_rad=30, max_conv_rad_km=5.0):
    """
    We perform the Yuter and Houze (1997) algorithm for echo classification
    using only the reflectivity field in order to classify each grid point
    as either convective, stratiform or undefined. Grid points are
    classified as follows,

    0 = No Surface Echo/ Undefined
    1 = Stratiform
    2 = Convective
    3 = Weak Echo

    refl : array
        array of reflectivity values
    x, y : array
        x and y coordinates of reflectivity array, respectively
    dx, dy : float
        The x- and y-dimension resolutions in meters, respectively.
    always_core_thres : float, optional
        Threshold for points that are always convective. All values above the threshold are classifed as convective
    bkg_rad_km : float, optional
        Radius to compute background reflectivity in kilometers. Default is 11 km
    use_cosine : bool, optional
        Boolean used to determine if cosine scheme should be used for identifying convective cores (True) or a scalar scheme (False)
    max_diff : float, optional
        Maximum difference between background average and reflectivity in order to be classified as convective.
        "a" value in Eqn. B1 in Yuter and Houze (1997)
    zero_diff_cos_val : float, optional
        Value where difference between background average and reflectivity is zero in the cosine function
        "b" value in Eqn. B1 in Yuter and Houze (1997)
    scalar_diff : float, optional
        If using a scalar difference scheme, this value is the multiplier or addition to the background average
    use_addition : bool, optional
        Determines if a multiplier (False) or addition (True) in the scalar difference scheme should be used
    calc_thres : float, optional
        Minimum percentage of points needed to be considered in background average calculation
    weak_echo_thres : float, optional
        Threshold for determining weak echo. All values below this threshold will be considered weak echo
    min_dBZ_used : float, optional
        Minimum dBZ value used for classification. All values below this threshold will be considered no surface echo
    dB_averaging : bool, optional
        True if using dBZ values that need to be converted to linear Z before averaging. False for other types of values
    apply_lg_rad_mask : bool, optional
        Flag to set a large radial mask for algorithm
    lg_rad_mask_min_rad_km, lg_rad_mask_max_rad_km : float, optional
        Values for setting the large radial mask
    val_for_max_conv_rad : float, optional
        dBZ for maximum convective radius. Convective cores with values above this will have the maximum convective radius
    max_conv_rad_km : float, optional
        Maximum radius around convective cores to classify as convective. Default is 5 km

   Returns
    -------
    refl_bkg : array
        Array of background values
    conv_core_array : array
        Array of initial convective cores (identified convective elements without convective radii applied)
    conv_strat_array : array
        Array of convective stratiform classifcation with convective radii applied
    """

    # Constants to fill arrays with
    CS_CORE = 3
    NOSFCECHO = 0
    WEAKECHO = 3
    SF = 1
    CONV = 2

    # %% Set up mask arrays for background average and
    # prepare for convective mask arrays
    # calculate maximum convective diameter from max. convective radius (input)
    max_conv_diameter = int(np.floor((max_conv_rad_km / (dx / 1000)) * 2))
    # if diameter is even, make odd
    if max_conv_diameter % 2 == 0:
        max_conv_diameter = max_conv_diameter + 1
    # find center point
    center_conv_mask_x = int(np.floor(max_conv_diameter / 2))

    # prepare background mask array for computing background average
    # calculate number of pixels for background array given requested background radius and dx
    bkg_diameter_pix = int(np.floor((bkg_rad_km / (dx / 1000)) * 2))
    # set diameter to odd if even
    if bkg_diameter_pix % 2 == 0:
        bkg_diameter_pix = bkg_diameter_pix + 1
    # find center point
    bkg_center = int(np.floor(bkg_diameter_pix / 2))
    # create background array
    bkg_mask_array = np.ones((bkg_diameter_pix, bkg_diameter_pix), dtype=float)
    # mask outside circular region
    bkg_mask_array = create_radial_mask(bkg_mask_array, min_rad_km=0, max_rad_km=bkg_rad_km, x_pixsize=dx / 1000,
                                        y_pixsize=dy / 1000, center_x=bkg_center, center_y=bkg_center, circular=True)

    # Create large mask array for determining where to calculate convective stratiform
    # initialize array with 1 (calculate convective stratiform over entire array)
    mask_array = np.zeros(refl.shape, dtype=float)
    mask_array[:] = 1
    # if True, create radial mask
    if apply_lg_rad_mask:
        mask_array = create_radial_mask(mask_array, lg_rad_mask_min_rad_km, lg_rad_mask_max_rad_km, x_pixsize=dx / 1000,
                                        y_pixsize=dy / 1000, center_x=int(np.floor(refl.shape[0] / 2)),
                                        center_y=int(np.floor(refl.shape[1] / 2)), circular=True)

    # %% Convective stratiform detection

    # Compute background radius
    refl_bkg = calc_bkg_intensity(refl, bkg_mask_array, dB_averaging, calc_thres)
    # mask background average
    refl_bkg = np.ma.masked_where(refl.mask, refl_bkg)

    # Get convective core array from cosine scheme, or scalar scheme
    if use_cosine:
        conv_core_array = convcore_cos_scheme(refl, refl_bkg, max_diff, zero_diff_cos_val, always_core_thres, CS_CORE)
    else:
        conv_core_array = convcore_scalar_scheme(refl, refl_bkg, scalar_diff, always_core_thres, CS_CORE,
                                                 use_addition=use_addition)

    # Assign convective radii based on background reflectivity
    conv_radius_km = assign_conv_radius_km(refl_bkg, val_for_max_conv_rad=val_for_max_conv_rad,
                                           max_conv_rad=max_conv_rad_km)

    # Incorporate convective radius using binary dilation
    # Create empty array for assignment
    temp_assignment = np.zeros_like(conv_core_array)

    # Loop through radii
    for radius in np.arange(1, max_conv_rad_km + 1):
        # create mask array for radius incorporation
        conv_mask_array = create_conv_radius_mask(max_conv_diameter, radius, dx / 1000, dy / 1000, center_conv_mask_x)
        # find location of radius
        temp = conv_radius_km == radius
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
                                                 min_dBZ_used, weak_echo_thres)

    return refl_bkg, conv_core_array, conv_strat_array


# functions

def create_radial_mask(mask_array, min_rad_km, max_rad_km, x_pixsize,
                       y_pixsize, center_x, center_y, circular=True):
    """
    Computes a radial distance mask, everything with distance between minradiuskm
    and maxradiuskm is assigned 1, everything else is assigned 0. This version can
    handle rectangular arrays and pixels as well as square ones.

    Parameters
    ----------
    mask_array : array
        Array to mask
    min_rad_km, max_rad_km : float
        The minimum and maximum radius of the non-masked region in kilometers.
    x_pixsize, y_pixsize : float
        The pixel size in the x- and y-dimension in kilometers, respectively
    center_x, center_y : int
        The center pixel in the x- and y-dimension, respectively
    circular : bool
        True returns circular mask, False returns a rectangular mask.

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
                x_range_sq = ((center_x - i) * x_pixsize) ** 2
                y_range_sq = ((center_y - j) * y_pixsize) ** 2
                circ_range = np.sqrt(x_range_sq + y_range_sq)
            # if circular is False, use square mask
            else:
                x_range = abs(int(np.floor(center_x - i) * x_pixsize))
                y_range = abs(int(np.floor(center_y - j) * y_pixsize))

                if x_range > y_range:
                    circ_range = x_range
                else:
                    circ_range = y_range
            # if range is within min and max, set to True
            if (circ_range <= max_rad_km) and (circ_range >= min_rad_km):
                mask_array[j, i] = 1
            else:
                mask_array[j, i] = 0

    return mask_array


def calc_bkg_intensity(refl, bkg_mask_array, dB_averaging, calc_thres=None):
    """
    Calculate the background of the given refl array. The footprint used to
    calculate the average for each pixel is given by bkg_mask_array

    Parameters
    ----------
    refl : array
        Reflectivity array to compute average
    bkg_mask_array : array
        Array of radial points to use for average
    dB_averaging : bool
        If True, converts dBZ to linear Z before averaging
    calc_thres : float
        Minimum percentage of points needed to be considered in background average calculation

    Returns
    -------
    refl_bkg : array
        Array of average values
    """

    # if dBaverage is true, convert reflectivity to linear Z
    if dB_averaging:
        refl = 10 ** (refl / 10)

    # calculate background reflectivity with circular footprint
    refl_bkg = scipy.ndimage.generic_filter(refl.filled(np.nan), function=np.nanmean, mode='constant',
                                            footprint=bkg_mask_array.astype(bool), cval=np.nan)

    # if calc_thres is not none, then calculate the number of points used to calculate average
    if calc_thres is not None:
        # count valid points
        refl_count = scipy.ndimage.generic_filter(refl.filled(0), function=np.count_nonzero, mode='constant',
                                                footprint=bkg_mask_array.astype(bool), cval=0)
        # find threshold number of points
        val = calc_thres * np.count_nonzero(bkg_mask_array)
        # mask out values
        refl_bkg = np.ma.masked_where(refl_count < val, refl_bkg)

    # mask where original reflectivity is invalid
    refl_bkg = np.ma.masked_where(refl.mask, refl_bkg)

    # if dBaveraging is true, convert background reflectivity to dBZ
    if dB_averaging:
        refl_bkg = 10 * np.log10(refl_bkg)
        # mask where original reflectivity is invalid
        refl_bkg = np.ma.masked_where(refl.mask, refl_bkg)

    return refl_bkg


def convcore_cos_scheme(refl, refl_bkg, max_diff, zero_diff_cos_val, always_core_thres, CS_CORE):
    """
    Function for assigning convective cores based on a cosine function

    Parameters
    ----------
    refl : array
        Reflectivity values
    refl_bkg : array
        Background average of reflectivity values
    max_diff : float
        Maximum difference between refl and refl_bkg needed for convective classification
    zero_diff_cos_val : float
        Value where the cosine function returns a zero difference
    always_core_thres : float
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
    zDiff = max_diff * np.cos((np.pi * refl_bkg) / (2 * zero_diff_cos_val))
    zDiff[zDiff < 0] = 0  # where difference less than zero, set to zero
    zDiff[refl_bkg < 0] = max_diff  # where background less than zero, set to max. diff

    # set values
    conv_core_array[refl >= always_core_thres] = CS_CORE  # where Z is greater than always_core_thres, set to core
    conv_core_array[(refl - refl_bkg) >= zDiff] = CS_CORE  # where difference exceeeds minimum, set to core

    return conv_core_array


def convcore_scalar_scheme(refl, refl_bkg, max_diff, always_core_thres, CS_CORE, use_addition=False):
    """
    Function for assigning convective cores based on a scalar difference

    Parameters
    ----------
    refl : array
        Reflectivity values
    refl_bkg : array
        Background average of reflectivity values
    max_diff : float
        Maximum difference between refl and refl_bkg needed for convective classification
    always_core_thres : float
        All values above this threshold considered to be convective
    CS_CORE : int
        Value assigned to convective pixels
    use_addition : bool
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
    if use_addition:
        zDiff = max_diff + refl_bkg
    else:
        zDiff = max_diff * refl_bkg

    zDiff[zDiff < 0] = 0  # where difference less than zero, set to zero
    zDiff[refl_bkg < 0] = 0  # where background less than zero, set to zero

    # set values
    conv_core_array[refl >= always_core_thres] = CS_CORE  # where Z is greater than always_core_thres, set to core
    conv_core_array[refl >= zDiff] = CS_CORE  # where difference exceeeds minimum, set to core

    return conv_core_array


def create_conv_radius_mask(max_conv_diameter, radius_km, x_spacing, y_spacing, center_conv_mask_x):
    """
    Does and initial convective stratiform classification

    Parameters
    ----------
    max_conv_diameter : int
        maximum convective diameter in kilometers
    radius_km : int
        convective radius in kilometers
    x_spacing, y_spacing : float
        x- and y-dimension pixel size in meters, respectively
    center_conv_mask_x : int
        index of center point

    Returns
    -------
   conv_mask_array : array
        array masked based on distance of convective diameter
    """

    conv_mask_array = np.zeros((max_conv_diameter, max_conv_diameter))
    conv_mask_array = create_radial_mask(conv_mask_array, 0, radius_km, x_spacing, y_spacing, center_conv_mask_x,
                                         center_conv_mask_x, True)

    return conv_mask_array


def assign_conv_radius_km(refl_bkg, val_for_max_conv_rad, max_conv_rad=5):
    """
    Assigns the convective radius in kilometers based on the background values

    Parameters
    ----------
    refl_bkg : array
        array of background reflectivity values
    val_for_max_conv_rad : float
        reflectivity value for maximum convective radius (5 km)
    max_conv_rad : float, optional
        maximum convective radius in kilometers

    Returns
    -------
    convRadiuskm : array
        array of convective radii based on background values and val for max. conv radius
    """

    convRadiuskm = np.ones_like(refl_bkg)

    convRadiuskm[refl_bkg >= (val_for_max_conv_rad - 15)] = max_conv_rad - 3
    convRadiuskm[refl_bkg >= (val_for_max_conv_rad - 10)] = max_conv_rad - 2
    convRadiuskm[refl_bkg >= (val_for_max_conv_rad - 5)] = max_conv_rad - 1
    convRadiuskm[refl_bkg >= val_for_max_conv_rad] = max_conv_rad

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
    conv_strat_array : array
        array of classifications
    """

    # assuming order so that each point is only assigned one time, no overlapping assignment
    # initially, assign every point to stratiform
    conv_strat_array[:] = SF
    # where reflectivity is masked, set to no surface echo
    conv_strat_array[refl.mask] = NOSFCECHO
    # assign convective cores to CONV
    conv_strat_array[conv_core_array == CS_CORE] = CONV
    # assign reflectivity less than weakechothres to weak echo
    conv_strat_array[refl < WEAKECHOTHRES] = WEAKECHO
    # assign reflectivity less than minimum to no surface echo
    conv_strat_array[refl < MINDBZUSE] = NOSFCECHO

    return conv_strat_array
