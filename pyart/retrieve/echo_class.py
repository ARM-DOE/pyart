"""
pyart.retrieve.echo_class
=========================

Functions for echo classification

.. autosummary::
    :toctree: generated/

    steiner_conv_strat
    hydroclass_semisupervised
    _standardize
    _assign_to_class

"""

import numpy as np

from ..config import get_fillvalue, get_field_name, get_metadata
from ..exceptions import MissingOptionalDependency
try:
    from . import _echo_steiner
    _F90_EXTENSIONS_AVAILABLE = True
except ImportError:
    _F90_EXTENSIONS_AVAILABLE = False

from warnings import warn

def steiner_conv_strat(grid, dx=None, dy=None, intense=42.0,
                       work_level=3000.0, peak_relation='default',
                       area_relation='medium', bkg_rad=11000.0,
                       use_intense=True, fill_value=None,
                       refl_field=None):
    """
    Partition reflectivity into convective-stratiform using the Steiner et
    al. (1995) algorithm.

    Parameters
    ----------
    grid : Grid
        Grid containing reflectivity field to partition.

    Other Parameters
    ----------------
    dx, dy : float
        The x- and y-dimension resolutions in meters, respectively.  If None
        the resolution is determined from the first two axes values.
    intense : float
        The intensity value in dBZ. Grid points with a reflectivity
        value greater or equal to the intensity are automatically
        flagged as convective. See reference for more information.
    work_level : float
        The working level (separation altitude) in meters. This is the height
        at which the partitioning will be done, and should minimize bright band
        contamination. See reference for more information.
    peak_relation : 'default' or 'sgp'
        The peakedness relation. See reference for more information.
    area_relation : 'small', 'medium', 'large', or 'sgp'
        The convective area relation. See reference for more information.
    bkg_rad : float
        The background radius in meters. See reference for more information.
    use_intense : bool
        True to use the intensity criteria.
    fill_value : float
         Missing value used to signify bad data points. A value of None
         will use the default fill value as defined in the Py-ART
         configuration file.
    refl_field : str
         Field in grid to use as the reflectivity during partitioning. None
         will use the default reflectivity field name from the Py-ART
         configuration file.

    Returns
    -------
    eclass : dict
        Steiner convective-stratiform classification dictionary.

    References
    ----------
    Steiner, M. R., R. A. Houze Jr., and S. E. Yuter, 1995: Climatological
    Characterization of Three-Dimensional Storm Structure from Operational
    Radar and Rain Gauge Data. J. Appl. Meteor., 34, 1978-2007.
    """
    # check that Fortran extensions is available
    if not _F90_EXTENSIONS_AVAILABLE:
        raise MissingOptionalDependency(
            "Py-ART must be built on a system with a Fortran compiler to "
            "use the steiner_conv_strat function.")

    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')

    # parse dx and dy
    if dx is None:
        dx = grid.x['data'][1] - grid.x['data'][0]
    if dy is None:
        dy = grid.y['data'][1] - grid.y['data'][0]

    # Get coordinates
    x = grid.x['data']
    y = grid.y['data']
    z = grid.z['data']

    # Get reflectivity data
    ze = np.ma.copy(grid.fields[refl_field]['data'])
    ze = np.ma.filled(ze, fill_value).astype(np.float64)

    # Call Fortran routine
    eclass = _echo_steiner.classify(
        ze, x, y, z, dx=dx, dy=dy, bkg_rad=bkg_rad, work_level=work_level,
        intense=intense, peak_relation=peak_relation,
        area_relation=area_relation, use_intense=use_intense,
        fill_value=fill_value)

    return {'data': eclass.astype(np.int32),
            'standard_name': 'echo_classification',
            'long_name': 'Steiner echo classification',
            'valid_min': 0,
            'valid_max': 2,
            'comment_1': ('Convective-stratiform echo '
                          'classification based on '
                          'Steiner et al. (1995)'),
            'comment_2': ('0 = Undefined, 1 = Stratiform, '
                          '2 = Convective')}


def hydroclass_semisupervised(radar, mass_centers=None,
                              weights=np.array([1., 1., 1., 0.75, 0.5]),
                              refl_field=None, zdr_field=None, rhv_field=None,
                              kdp_field=None, temp_field=None,
                              hydro_field=None):
    """
    Classifies precipitation echoes following the approach by
    Besic et al (2016)

    Parameters
    ----------
    radar : radar
        radar object

    Other Parameters
    ----------------
    mass_centers : matrix
        The centroids for each variable and hydrometeor class in (nclasses,
        nvariables)
    weights : vector
        The weight given to each variable.
    refl_field, zdr_field, rhv_field, kdp_field, temp_field : str
        Inputs. Field names within the radar object which represent the
        horizonal reflectivity, the differential reflectivity, the copolar
        correlation coefficient, the specific differential phase and the
        temperature field. A value of None for any of these parameters will
        use the default field name as defined in the Py-ART configuration
        file.
    hydro_field : str
        Output. Field name which represents the hydrometeor class field.
        A value of None will use the default field name as defined in the
        Py-ART configuration file.

    Returns
    -------
    hydro : dict
        hydrometeor classification field

    References
    ----------
    Besic, N., Figueras i Ventura, J., Grazioli, J., Gabella, M., Germann, U.,
    and Berne, A.: Hydrometeor classification through statistical clustering
    of polarimetric radar measurements: a semi-supervised approach,
    Atmos. Meas. Tech., 9, 4425-4445, doi:10.5194/amt-9-4425-2016, 2016

    """
    nclasses = 9
    nvariables = 5
    lapse_rate = -6.5

    # select the centroids as a function of frequency band
    if mass_centers is None:
        mass_centers = np.zeros((nclasses, nvariables))
        # assign coefficients according to radar frequency
        if 'frequency' in radar.instrument_parameters:
            freq = radar.instrument_parameters['frequency']['data'][0]
            # C band
            if freq >= 4e9 and freq < 8e9:
                freq_band = 'C'
                # centroids derived for MeteoSwiss Albis radar
                #                     Zh      ZDR     kdp   RhoHV   delta_Z
                mass_centers[0, :] = [
                    13.5829,  0.4063, 0.0497, 0.9868,  1330.3]  # DS
                mass_centers[1, :] = [
                    02.8453,  0.2457, 0.0000, 0.9798,  0653.8]  # CR
                mass_centers[2, :] = [
                    07.6597,  0.2180, 0.0019, 0.9799, -1426.5]  # LR
                mass_centers[3, :] = [
                    31.6815,  0.3926, 0.0828, 0.9978,  0535.3]  # GR
                mass_centers[4, :] = [
                    39.4703,  1.0734, 0.4919, 0.9876, -1036.3]  # RN
                mass_centers[5, :] = [
                    04.8267, -0.5690, 0.0000, 0.9691,  0869.8]  # VI
                mass_centers[6, :] = [
                    30.8613,  0.9819, 0.1998, 0.9845, -0066.1]  # WS
                mass_centers[7, :] = [
                    52.3969,  2.1094, 2.4675, 0.9730, -1550.2]  # MH
                mass_centers[8, :] = [
                    50.6186, -0.0649, 0.0946, 0.9904,  1179.9]  # IH/HDG
            # X band
            elif freq >= 8e9 and freq <= 12e9:
                freq_band = 'X'
                # centroids derived for MeteoSwiss DX50 radar
                #      Zh      ZDR     kdp   RhoHV   delta_Z
                mass_centers[0, :] = [
                    19.0770,  0.4139, 0.0099, 0.9841,  1061.7]  # DS
                mass_centers[1, :] = [
                    03.9877,  0.5040, 0.0000, 0.9642,  0856.6]  # CR
                mass_centers[2, :] = [
                    20.7982,  0.3177, 0.0004, 0.9858, -1375.1]  # LR
                mass_centers[3, :] = [
                    34.7124, -0.3748, 0.0988, 0.9828,  1224.2]  # GR
                mass_centers[4, :] = [
                    33.0134,  0.6614, 0.0819, 0.9802, -1169.8]  # RN
                mass_centers[5, :] = [
                    08.2610, -0.4681, 0.0000, 0.9722,  1100.7]  # VI
                mass_centers[6, :] = [
                    35.1801,  1.2830, 0.1322, 0.9162, -0159.8]  # WS
                mass_centers[7, :] = [
                    52.4539,  2.3714, 1.1120, 0.9382, -1618.5]  # MH
                mass_centers[8, :] = [
                    44.2216, -0.3419, 0.0687, 0.9683,  1272.7]  # IH/HDG
            else:
                if freq < 4e9:
                    freq_band = 'C'
                    # centroids derived for MeteoSwiss Albis radar
                    #       Zh      ZDR     kdp   RhoHV   delta_Z
                    mass_centers[0, :] = [
                        13.5829,  0.4063, 0.0497, 0.9868,  1330.3]  # DS
                    mass_centers[1, :] = [
                        02.8453,  0.2457, 0.0000, 0.9798,  0653.8]  # CR
                    mass_centers[2, :] = [
                        07.6597,  0.2180, 0.0019, 0.9799, -1426.5]  # LR
                    mass_centers[3, :] = [
                        31.6815,  0.3926, 0.0828, 0.9978,  0535.3]  # GR
                    mass_centers[4, :] = [
                        39.4703,  1.0734, 0.4919, 0.9876, -1036.3]  # RN
                    mass_centers[5, :] = [
                        04.8267, -0.5690, 0.0000, 0.9691,  0869.8]  # VI
                    mass_centers[6, :] = [
                        30.8613,  0.9819, 0.1998, 0.9845, -0066.1]  # WS
                    mass_centers[7, :] = [
                        52.3969,  2.1094, 2.4675, 0.9730, -1550.2]  # MH
                    mass_centers[8, :] = [
                        50.6186, -0.0649, 0.0946, 0.9904,  1179.9]  # IH/HDG
                else:
                    freq_band = 'X'
                    # centroids derived for MeteoSwiss DX50 radar
                    #       Zh      ZDR     kdp   RhoHV   delta_Z
                    mass_centers[0, :] = [
                        19.0770,  0.4139, 0.0099, 0.9841,  1061.7]  # DS
                    mass_centers[1, :] = [
                        03.9877,  0.5040, 0.0000, 0.9642,  0856.6]  # CR
                    mass_centers[2, :] = [
                        20.7982,  0.3177, 0.0004, 0.9858, -1375.1]  # LR
                    mass_centers[3, :] = [
                        34.7124, -0.3748, 0.0988, 0.9828,  1224.2]  # GR
                    mass_centers[4, :] = [
                        33.0134,  0.6614, 0.0819, 0.9802, -1169.8]  # RN
                    mass_centers[5, :] = [
                        08.2610, -0.4681, 0.0000, 0.9722,  1100.7]  # VI
                    mass_centers[6, :] = [
                        35.1801,  1.2830, 0.1322, 0.9162, -0159.8]  # WS
                    mass_centers[7, :] = [
                        52.4539,  2.3714, 1.1120, 0.9382, -1618.5]  # MH
                    mass_centers[8, :] = [
                        44.2216, -0.3419, 0.0687, 0.9683,  1272.7]  # IH/HDG
                warn('Radar frequency out of range. ' +
                     'Centroids only valid for C or X band. ' +
                     freq_band + ' band centroids will be applied')
        else:
            # centroids derived for MeteoSwiss Albis radar
            #       Zh      ZDR     kdp   RhoHV   delta_Z
            mass_centers[0, :] = [
                13.5829,  0.4063, 0.0497, 0.9868,  1330.3]  # DS
            mass_centers[1, :] = [
                02.8453,  0.2457, 0.0000, 0.9798,  0653.8]  # CR
            mass_centers[2, :] = [
                07.6597,  0.2180, 0.0019, 0.9799, -1426.5]  # LR
            mass_centers[3, :] = [
                31.6815,  0.3926, 0.0828, 0.9978,  0535.3]  # GR
            mass_centers[4, :] = [
                39.4703,  1.0734, 0.4919, 0.9876, -1036.3]  # RN
            mass_centers[5, :] = [
                04.8267, -0.5690, 0.0000, 0.9691,  0869.8]  # VI
            mass_centers[6, :] = [
                30.8613,  0.9819, 0.1998, 0.9845, -0066.1]  # WS
            mass_centers[7, :] = [
                52.3969,  2.1094, 2.4675, 0.9730, -1550.2]  # MH
            mass_centers[8, :] = [
                50.6186, -0.0649, 0.0946, 0.9904,  1179.9]  # IH/HDG
            warn('Radar frequency unknown. ' +
                 'Default coefficients for C band will be applied')

    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if zdr_field is None:
        zdr_field = get_field_name('differential_reflectivity')
    if rhv_field is None:
        rhv_field = get_field_name('cross_correlation_ratio')
    if kdp_field is None:
        kdp_field = get_field_name('specific_differential_phase')
    if temp_field is None:
        temp_field = get_field_name('temperature')
    if hydro_field is None:
        hydro_field = get_field_name('radar_echo_classification')

    # extract fields and parameters from radar
    if refl_field in radar.fields:
        refl = radar.fields[refl_field]['data']
    else:
        raise KeyError('Field not available: ' + refl_field)
    if zdr_field in radar.fields:
        zdr = radar.fields[zdr_field]['data']
    else:
        raise KeyError('Field not available: ' + zdr_field)
    if rhv_field in radar.fields:
        rhohv = radar.fields[rhv_field]['data']
    else:
        raise KeyError('Field not available: ' + rhv_field)
    if kdp_field in radar.fields:
        kdp = radar.fields[kdp_field]['data']
    else:
        raise KeyError('Field not available: ' + kdp_field)
    if temp_field in radar.fields:
        temp = radar.fields[temp_field]['data']
    else:
        raise KeyError('Field not available: ' + temp_field)

    # convert temp in relative height respect to iso0
    relh = temp*(1000./lapse_rate)

    # data mask
    # mask_zh = refl.mask

    # standardize data
    refl_std = _standardize(refl, 'Zh')
    zdr_std = _standardize(zdr, 'ZDR')
    kdp_std = _standardize(kdp, 'KDP')
    rhohv_std = _standardize(rhohv, 'RhoHV')
    relh_std = _standardize(relh, 'relH')

    # standardize centroids
    mc_std = np.zeros((nclasses, nvariables))
    mc_std[:, 0] = _standardize(mass_centers[:, 0], 'Zh')
    mc_std[:, 1] = _standardize(mass_centers[:, 1], 'ZDR')
    mc_std[:, 2] = _standardize(mass_centers[:, 2], 'KDP')
    mc_std[:, 3] = _standardize(mass_centers[:, 3], 'RhoHV')
    mc_std[:, 4] = _standardize(mass_centers[:, 4], 'relH')

    # assign to class
    hydroclass_data, min_dist = _assign_to_class(
        refl_std, zdr_std, kdp_std, rhohv_std, relh_std, mc_std,
        weights=weights)

    # prepare output fields
    hydro = get_metadata(hydro_field)
    hydro['data'] = hydroclass_data

    return hydro


def _standardize(data, field_name):
    """
    Streches the radar data to -1 to 1 interval

    Parameters
    ----------
    data : array
        radar field

    field_name : str
        type of field (relH, Zh, ZDR, KDP or RhoHV)

    Returns
    -------
    field_std : dict
        standardized radar data
    """
    if field_name == 'relH':
        field_std = 2./(1.+np.ma.exp(-0.005*data))-1.
        return field_std

    if field_name == 'Zh':
        mx = 60.
        mn = -10.
    elif field_name == 'ZDR':
        mx = 5.
        mn = -5.
    elif field_name == 'KDP':
        mx = 7.
        mn = -10.

        data[data < -0.5] = -0.5
        data = 10.*np.ma.log10(data+0.6)
    elif field_name == 'RhoHV':
        mx = -5.23
        mn = -50.
        data = 10.*np.ma.log10(1.-data)

    mask = np.ma.getmaskarray(data)
    field_std = 2.*(data-mn)/(mx-mn)-1.
    field_std[data < mn] = -1.
    field_std[data > mx] = 1.
    field_std[mask] = np.ma.masked

    return field_std


def _assign_to_class(zh, zdr, kdp, rhohv, relh, mass_centers,
                    weights=np.array([1., 1., 1., 0.75, 0.5])):
    """
    assigns an hydrometeor class to a radar range bin computing
    the distance between the radar variables an a centroid

    Parameters
    ----------
    zh,zdr,kdp,rhohv,relh : radar field
        variables used for assigment normalized to [-1, 1] values

    mass_centers : matrix
        centroids normalized to [-1, 1] values

    weights : array
        optional. The weight given to each variable

    Returns
    -------
    hydroclass : int array
        the index corresponding to the assigned class
    mind_dist : float array
        the minimum distance to the centroids
    """
    # prepare data
    nrays = zh.shape[0]
    nbins = zdr.shape[1]
    nclasses = mass_centers.shape[0]
    nvariables = mass_centers.shape[1]

    data = np.ma.array([zh, zdr, kdp, rhohv, relh])
    weights_mat = np.broadcast_to(
        weights.reshape(nvariables, 1, 1),
        (nvariables, nrays, nbins))
    dist = np.ma.zeros((nclasses, nrays, nbins), dtype='float64')

    # compute distance: masked entries will not contribute to the distance
    for i in range(nclasses):
        centroids_class = mass_centers[i, :]
        centroids_class = np.broadcast_to(
            centroids_class.reshape(nvariables, 1, 1),
            (nvariables, nrays, nbins))
        dist[i, :, :] = np.ma.sqrt(np.ma.sum(
            ((centroids_class-data)**2.)*weights_mat, axis=0))

    # use very large fill_value so that masked entries will be sorted at the
    # end. There should not be any masked entry anyway
    class_vec = dist.argsort(axis=0, fill_value=10e40)

    # get minimum distance. Acts as a confidence value
    dist_sorted = dist.sort(axis=0, fill_value=10e40)
    min_dist = dist[0, :, :]

    # Entries with non-valid reflectivity values are set to 0 (No class)
    mask = np.ma.getmaskarray(zh)
    hydroclass = class_vec[0, :, :]+1
    hydroclass[mask] = 0

    return hydroclass, min_dist
