"""
pyart.io.cdm
============

Functions for accessing Common Data Model (CDM) NEXRAD Level 2 files.

.. autosummary::
    :toctree: generated/

    read_nexrad_cdm
    _gen_vnames
    _radar_from_cdm

"""

import netCDF4
import numpy as np

from .radar import Radar
from .common import get_metadata

# XXX move this into a nexrad_common module
# METADATA specific to NEXRAD fields
NEXRAD_METADATA = {

    'reflectivity': {
        'units': 'dBZ',
        'standard_name': 'equivalent_reflectivity_factor',
        'long_name': 'equivalent_reflectivity_factor',
        'valid_max': 94.5,
        'valid_min': -32.0,
        'coordinates': 'elevation azimuth range'},

    'velocity': {
        'units': 'meters_per_second',
        'standard_name': (
            'radial_velocity_of_scatterers_away_from_instrument'),
        'long_name': (
            'radial_velocity_of_scatterers_away_from_instrument'),
        'valid_max': 95.0,          # XXX
        'valid_min': -95.0,         # XXX
        'coordinates': 'elevation azimuth range'},

    'spectrum_width': {
        'units': 'meters_per_second',
        'standard_name': 'spectrum_width',
        'long_name': 'spectrum_width',
        'valid_max': 63.0,
        'valid_min': -63.5,
        'coordinates': 'elevation azimuth range'},

    'differential_reflectivity': {
        'units': 'dB',
        'standard_name': 'log_differential_reflectivity_hv',
        'long_name': 'log_differential_reflectivity_hv',
        'valid_max': 7.9375,
        'valid_min': -7.8750,
        'coordinates': 'elevation azimuth range'},

    'differential_phase': {
        'units': 'degrees',
        'standard_name': 'differential_phase_hv',
        'long_name': 'differential_phase_hv',
        'valid_max': 360.0,
        'valid_min': 0.0,
        'coordinates': 'elevation azimuth range'},

    'correlation_coefficient': {
        'units': 'ratio',
        'standard_name': 'cross_correlation_ratio_hv',
        'long_name': 'cross_correlation_ratio_hv',
        'valid_max': 1.0,
        'valid_min': 0.0,
        'coordinates': 'elevation azimuth range'},
}


def read_nexrad_cdm(filename):
    """
    Read a Common Data Model (CDM) NEXRAD Level 2 file.

    Parameters
    ----------
    filename : str
        File name or URL of a Common Data Model (CDM) NEXRAD Level 2 file.
        This file can be created using NetCDF Java Library tools [1]_,
        or a URL of a OPeNDAP file on the UCAR THREDDS Data Server [2]_.

    Returns
    -------
    refl_hi : Radar
        Radar object containing the super resolution reflectivity field.
    dopl_hi : Radar
        Radar object containing the super resolution doppler fields present.
        Possible fields in the radar are: velocity, spectrum_width,
        differential_reflectivity, differential_phase, and
        correlation_coefficient.
    refl_sd : Radar
        Radar object containing the standard resolution reflectivity field.
    dopl_sd : Radar
        Radar object containing the standard resolution doppler fields.

    References
    ----------
    .. [1] http://www.unidata.ucar.edu/software/netcdf-java/documentation.htm
    .. [2] http://thredds.ucar.edu/thredds/catalog.html

    """
    dataset = netCDF4.Dataset(filename)
    dattrs = dataset.ncattrs()
    if 'cdm_data_type' not in dattrs or dataset.cdm_data_type != 'RADIAL':
        raise IOError('%s is not a valid CDM NetCDF file' % (filename))

    # Might need to add a check to see if all fields/resolution are present.
    refl_hi = _radar_from_cdm(dataset, _gen_vnames('refl', True))
    dopl_hi = _radar_from_cdm(dataset, _gen_vnames('doppler', True))
    refl_sd = _radar_from_cdm(dataset, _gen_vnames('refl', False))
    dopl_sd = _radar_from_cdm(dataset, _gen_vnames('doppler', False))
    return refl_hi, dopl_hi, refl_sd, dopl_sd


# mappings from CDM dataset variables to Radar object field names
CDM_FIELD_MAPPING = {
    'Reflectivity_HI': 'reflectivity',
    'RadialVelocity_HI': 'velocity',
    'SpectrumWidth_HI': 'spectrum_width',
    'DifferentialReflectivity_HI': 'differential_reflectivity',
    'DifferentialPhase_HI': 'differential_phase',
    'CorrelationCoefficient_HI': 'correlation_coefficient',
    'Reflectivity': 'reflectivity',
    'RadialVelocity': 'velocity',
    'SpectrumWidth': 'spectrum_width',
    'DifferentialReflectivity': 'differential_reflectivity',
    'DifferentialPhase': 'differential_phase',
    'CorrelationCoefficient': 'correlation_coefficient'
}


def _gen_vnames(field, high_res):
    """
    Return a dictionary of variables/dimension names in a thredds dataset.

    Parameters
    ----------
    field : 'refl' or 'doppler'
        Type of fields to generate dictionary for.
    high_res: bool
        True if the fields are high resolution.

    Returns
    -------
    vnames : dic
        Dictionary of variable and dimension names in the thredds server
        dataset for a given field type.  For use as the vnames parameter in
        :py:func:`_radar_from_cdm`.

    """
    if field == 'refl':
        suffix = 'R'
        fields = ['Reflectivity']
    else:
        suffix = 'V'
        fields = ['RadialVelocity', 'SpectrumWidth',
                  'DifferentialReflectivity', 'CorrelationCoefficient',
                  'DifferentialPhase']
        # The RadialVelocity and SpectrumWidth fields have
        # variables/coordinate suffixed by a V, where as the
        # DifferentialReflectivity, CorrelationCoefficient, and
        # DifferentialPhase fields have their own suffixes (D, C, P)
        # These all appear to have the values, but  this has not
        # been verified and may by incorrect.  Going to this assumption
        # these fields are refered to in this function as doppler fields.

    if high_res:
        suffix += '_HI'
        fields = [f + '_HI' for f in fields]

    # build the variables as foo{R, V, D, C, P}{_HI}.
    keys = ['scan', 'radial', 'gate', 'time', 'elevation', 'azimuth',
            'distance', 'numRadials', 'numGates']
    vnames = dict([(k, k + suffix) for k in keys])
    vnames['fields'] = fields
    return vnames


def _radar_from_cdm(dataset, vnames):
    """
    Return a Radar from a CDM dataset using the variable and dimension names
    in the vnames dictionary.
    """
    dvars = dataset.variables

    # dimensionality
    nsweeps = len(dataset.dimensions[vnames['scan']])
    rays_per_sweep = len(dataset.dimensions[vnames['radial']])
    ngates = len(dataset.dimensions[vnames['gate']])
    nrays = nsweeps * rays_per_sweep

    # time
    time = get_metadata('time')
    # time in msecs, convert to seconds
    time['data'] = dvars[vnames['time']][:].reshape(nrays) / 1000.

    # change unit string from 'msecs since ...' to 'seconds since ...'
    units = dvars[vnames['time']].units
    time['units'] = 'seconds ' + units[6:]

    # _range
    _range = get_metadata('range')
    _range['data'] = dvars[vnames['distance']][:]

    # fields
    fields = {}
    for field in vnames['fields']:

        if field not in dvars:  # do nothing if field is not present
            continue

        field_name = CDM_FIELD_MAPPING[field]
        field_dic = NEXRAD_METADATA[field_name].copy()

        # extract the field
        fvar = dvars[field]
        fvar.set_auto_maskandscale(False)
        raw_fdata = fvar[:]

        # mask, reshape, scale and offset
        raw_fdata = np.ma.masked_less_equal(raw_fdata, 1)
        raw_fdata = raw_fdata.reshape(-1, raw_fdata.shape[-1])

        if 'scale_factor' in fvar.ncattrs():
            scale = fvar.scale_factor
        else:
            scale = 1.0

        if 'add_offset' in fvar.ncattrs():
            add_offset = fvar.add_offset
        else:
            add_offset = 0.0

        field_dic['data'] = raw_fdata * scale + add_offset

        # add the field dictionary to the fields dictionary
        fields[field_name] = field_dic

    # metadata
    metadata = {'original_container': 'NEXRAD Level 2 Archive - Thredds'}
    for key in dataset.ncattrs():
        metadata[key] = dataset.getncattr(key)

    # additional required CF/Radial metadata, blank if not present
    if 'Title' in dataset.ncattrs():
        metadata['title'] = dataset.Title
    else:
        metadata['title'] = ''
    metadata['institution'] = ''
    metadata['references'] = ''
    metadata['source'] = ''
    metadata['comment'] = ''
    metadata['instrument_name'] = ''

    # scan_type
    scan_type = 'ppi'

    # latitude, longitude, altitude
    latitude = get_metadata('latitude')
    longitude = get_metadata('longitude')
    altitude = get_metadata('altitude')
    latitude['data'] = np.array([dataset.StationLatitude], dtype='float64')
    longitude['data'] = np.array([dataset.StationLongitude], dtype='float64')
    altitude['data'] = np.array([dataset.StationElevationInMeters],
                                dtype='float64')

    # sweep_number, sweep_mode, sweep_start_ray_index,sweep_end_ray_index
    sweep_number = get_metadata('sweep_number')
    sweep_mode = get_metadata('sweep_mode')
    sweep_start_ray_index = get_metadata('sweep_start_ray_index')
    sweep_end_ray_index = get_metadata('sweep_end_ray_index')

    nsweeps = len(dataset.dimensions['scanR_HI'])
    rays_per_sweep = len(dataset.dimensions[vnames['radial']])

    sweep_number['data'] = np.arange(nsweeps, dtype='int32')
    sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])
    ssri = np.arange(nsweeps, dtype='int32') * rays_per_sweep
    seri = np.arange(1, nsweeps + 1, dtype='int32') * rays_per_sweep - 1
    sweep_start_ray_index['data'] = ssri
    sweep_end_ray_index['data'] = seri

    # azimuth, elevation, fixed_angle
    azimuth = get_metadata('azimuth')
    elevation = get_metadata('elevation')
    fixed_angle = get_metadata('fixed_angle')

    azimuth['data'] = dvars[vnames['azimuth']][:].reshape(-1)
    elevation['data'] = dvars[vnames['elevation']][:].reshape(-1)
    fixed_angle['data'] = dvars[vnames['elevation']][:].reshape(-1)

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=None)
