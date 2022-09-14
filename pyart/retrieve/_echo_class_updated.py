import numpy as np
import sys
sys.path.append('C:\\Users\\lmtomkin\\Documents\\GitHub\\pyart_convsf\\pyart\\retrieve\\')

import yuter_convsf
import time
import scipy.ndimage

def _revised_conv_strat(refl, dx, dy, alwaysConvThres=42, bkgRad_km=11,
                        useCosine=True, minZeDiff=8, convThresB=55, scalarDiff=1.5,
                        weakEchoThres=5.0, minDBZused=5.0, applyLargeRadialMask=False,
                        largeRadialMask_minRadkm=0, largeRadialMask_maxRadkm=170,
                        dBZforMaxConvRadius = 30, maxConvRad_km = 5.0, dBaveraging=False, incorp_rad=False):
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

    #%% Set up arrays and values for convective stratiform algorithm
    # create Ze arrays for under and overestimate
    # need to run through 3 times - need to figure out best way to do this
    # refl = grid.fields['reflectivity']['data'][0,:,:]
    # refl = np.ma.masked_invalid(refl)
    # refl_under = refl.copy() - 5
    # refl_over = refl.copy() + 5
     # loop through 3 times
    t1 = time.time()
    # create empty arrays
    ze_bkg = np.zeros(refl.shape, dtype=float)
    conv_core_array = np.zeros(refl.shape, dtype=float)
    conv_strat_array = np.zeros(refl.shape, dtype=float)
    mask_array = np.zeros(refl.shape, dtype=float)

    # fill with missing (nan)
    ze_bkg[:] = np.nan
    conv_core_array[:] = np.nan
    conv_strat_array[:] = np.nan
    mask_array[:] = np.nan

    # Constants to fill arrays with
    CS_CORE = 3
    NOSFCECHO = 0
    WEAKECHO = 3
    SF = 1
    CONV = 2
    t2 = time.time()
    print("Time to set up arrays: {0} seconds".format(t2-t1))
    #%% Set up mask arrays
    t1 = time.time()
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
    bkg_mask_array = yuter_convsf.radialDistanceMask_array(bkg_mask_array, minradiuskm=0, maxradiuskm=bkgRad_km,
                                                           x_pixsize=dx / 1000, y_pixsize=dy / 1000,
                                                           centerx=bkg_center, centery=bkg_center, circular=True)

    # Create large mask array for determining where to calculate convective stratiform
    # initialize array with 1 (calculate convective stratiform over entire array)
    mask_array[:] = 1
    # if True, create radial mask
    if applyLargeRadialMask:
        mask_array = yuter_convsf.radialDistanceMask_array(mask_array, largeRadialMask_minRadkm, largeRadialMask_maxRadkm,
                                                           x_pixsize=dx/1000, y_pixsize=dy/1000, centerx=int(np.floor(refl.shape[0] / 2)),
                                                           centery=int(np.floor(refl.shape[1] / 2)), circular=True)

    t2 = time.time()
    print("Time to set up mask arrays: {0} seconds".format(t2 - t1))
    #%% Convective stratiform detection
    t1 = time.time()
    # Compute background radius
    ze_bkg = yuter_convsf.backgroundIntensity_array(refl, bkg_mask_array, dBaveraging)
    # mask background average
    ze_bkg = np.ma.masked_where(refl.mask, ze_bkg)

    t2 = time.time()
    print("Time to calculate background average: {0} seconds".format(t2 - t1))

    t1 = time.time()
    # Get convective core array from cosine scheme, or scalar scheme
    if useCosine:
        conv_core_array = yuter_convsf.convcore_cos_scheme_array(refl, ze_bkg, minZeDiff, convThresB, alwaysConvThres,
                                                                 CS_CORE)
    else:
        conv_core_array = yuter_convsf.convcore_scaled_array(refl, ze_bkg, scalarDiff, alwaysConvThres, CS_CORE)

    t2 = time.time()
    print("Time to get convective cores: {0} seconds".format(t2 - t1))
    # count convective cores
    corecount = np.count_nonzero(conv_core_array)

    t1 = time.time()
    # Do an initial assignment of convsf array
    conv_strat_array = yuter_convsf.classify_conv_strat_array(refl, conv_strat_array, conv_core_array,
                                                              NOSFCECHO, CONV, SF, WEAKECHO, CS_CORE,
                                                              minDBZused, weakEchoThres)
    t2 = time.time()
    print("Time to do initial convective stratiform assignment: {0} seconds".format(t2 - t1))

    # Assign convective radii based on background reflectivity
    convRadiuskm = yuter_convsf.assignConvRadiuskm_array(ze_bkg, dBZformaxconvradius=dBZforMaxConvRadius, maxConvRadius=maxConvRad_km)

    if incorp_rad:
        t1 = time.time()
        # Loop through array for a final time to incorporate the convective radii
        for j in np.arange(0, convRadiuskm.shape[1], 1):

            for i in np.arange(0, convRadiuskm.shape[0], 1):

                # if point is a convective core, find radius, get convective mask radius and incorporate radius
                if conv_core_array[j, i] == CS_CORE:
                    convRadius = np.floor(convRadiuskm[j,i])

                    if convRadius == 1:
                        conv_mask_array = yuter_convsf.init_conv_radius_mask(maxConvDiameter, 1, dx / 1000, dy / 1000, centerConvMask_x)
                    elif convRadius > 1:
                        if convRadius == 2:
                            conv_mask_array = yuter_convsf.init_conv_radius_mask(maxConvDiameter, 2, dx / 1000, dy / 1000, centerConvMask_x)
                        elif convRadius == 3:
                            conv_mask_array = yuter_convsf.init_conv_radius_mask(maxConvDiameter, 3, dx / 1000, dy / 1000, centerConvMask_x)
                        elif convRadius == 4:
                            conv_mask_array = yuter_convsf.init_conv_radius_mask(maxConvDiameter, 4, dx / 1000, dy / 1000, centerConvMask_x)
                        elif convRadius == 5:
                            conv_mask_array = yuter_convsf.init_conv_radius_mask(maxConvDiameter, 5, dx / 1000, dy / 1000, centerConvMask_x)

                    conv_strat_array = yuter_convsf.incorporateConvRadius(conv_strat_array, i, j, maxConvDiameter,
                                                                          conv_mask_array, NOSFCECHO, CONV)

        t2 = time.time()
        print("Time to apply convective radius: {0} seconds".format(t2 - t1))

    return ze_bkg, conv_core_array, conv_strat_array

    # go through map array once to compute the background intensity and find convective cores
    # pixel by pixel method
    # for j in np.arange(0, ysize_refl, 1):
    #     for i in np.arange(0, xsize_refl, 1):
    #
    #         # compute background value
    #         ze_bkg[j, i] = yuter_convsf.backgroundIntensity(refl, i, j, bkgDiameter_pix, bkg_mask_array, dBaveraging)
    #
    #         # identify convective cores
    #         if use_cosine:
    #             conv_core_array[j, i] = yuter_convsf.convcore_cos_scheme(refl[j, i], ze_bkg[j, i], minZeDiff, convThresB, alwaysConvThres, CS_CORE)

    # else:
    #     mask_array = yuter_convsf.radialDistanceMask_array(mask_array, minRad_km, maxRad_km,
    #                                                        x_pixsize=dx/1000, y_pixsize=dy/1000, centerx=int(np.floor(xsize_refl / 2)),
    #                                                        centery=int(np.floor(ysize_refl / 2)), circular=False)

    # initialize mask arrays (1-5 km)
    # convmaskarray1 = yuter_convsf.init_conv_radius_mask(maxConvDiameter, 1, dx / 1000, dy / 1000, centerConvMask_x)
    # convmaskarray2 = yuter_convsf.init_conv_radius_mask(maxConvDiameter, 2, dx / 1000, dy / 1000, centerConvMask_x)
    # convmaskarray3 = yuter_convsf.init_conv_radius_mask(maxConvDiameter, 3, dx / 1000, dy / 1000, centerConvMask_x)
    # convmaskarray4 = yuter_convsf.init_conv_radius_mask(maxConvDiameter, 4, dx / 1000, dy / 1000, centerConvMask_x)
    # convmaskarray5 = yuter_convsf.init_conv_radius_mask(maxConvDiameter, 5, dx / 1000, dy / 1000, centerConvMask_x)
