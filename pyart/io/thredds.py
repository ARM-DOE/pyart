"""

"""

import netCDF4
import numpy as np

from .radar import Radar
from .common import get_metadata, make_time_unit_str


def read_nexrad_thredds_url(url):
    """
    Read a NEXRAD Level II dataset from the Thredds server returning a Radar
    object.

    """
    dataset = netCDF4.Dataset(url)
    return _radar_from_thredds_dataset(dataset)


NEXRADFIELDS = {
    'Reflectivity_HI': 'reflectivity',
    'RadialVelocity_HI': 'velocity',
    'SpectrumWidth_HI': 'spectrum_width',
    'DifferentialReflectivity_HI': 'differential_reflectivity',
    'DifferentialPhase_HI': 'differential_phase',
    'CorrelationCoefficient_HI': 'correlation_coefficient'
}

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


def _radar_from_thredds_dataset(dataset, high_res=True):
    """
    """
    dvars = dataset.variables

    # time
    time = get_metadata('time')
    # time in msecs since, convert to seconds since
    time['data'] = dvars['timeR_HI'][:].reshape(-1) / 1000.
    units = dvars['timeR_HI'].units
    time['units'] = 'seconds ' + units[6:]

    # _range
    _range = get_metadata('range')
    _range['data'] = dvars['distanceR_HI'][:]

    # fields
    high_res_fields = [
        'Reflectivity_HI', 'RadialVelocity_HI',
        'SpectrumWidth_HI', 'DifferentialReflectivity_HI',
        'CorrelationCoefficient_HI', 'DifferentialPhase_HI']
    fields = {}
    for moment in high_res_fields:
        if moment not in dvars:
            continue    # moment is not present, continue on
        field_name = NEXRADFIELDS[moment]
        field_dic = NEXRAD_METADATA[field_name].copy()

        # extract the field
        fvar = dvars[moment]
        fvar.set_auto_maskandscale(False)
        raw_fdata = fvar[:]

        # mask, reshape, scale and offset
        raw_fdata = np.ma.masked_less_equal(raw_fdata, 1)
        raw_fdata = raw_fdata.reshape(-1, raw_fdata.shape[-1])

        if 'scale' in fvar.ncattrs():
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
    metadata = {'original_container': 'NEXRAD Level II Archive - Thredds'}
    # additional required CF/Radial metadata set to blank strings
    metadata['title'] = ''
    metadata['institution'] = ''
    metadata['references'] = ''
    metadata['source'] = ''
    metadata['comment'] = ''
    metadata['instrument_name'] = 'NEXRAD'

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

    # sweep_number, sweep_mode, sweep_start_ray_index,
    # sweep_end_ray_index
    sweep_number = get_metadata('sweep_number')
    sweep_mode = get_metadata('sweep_mode')
    sweep_start_ray_index = get_metadata('sweep_start_ray_index')
    sweep_end_ray_index = get_metadata('sweep_end_ray_index')

    nsweeps = len(dataset.dimensions['scanR_HI'])
    rays_per_sweep = len(dataset.dimensions['radialR_HI'])

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

    azimuth['data'] = dvars['azimuthR_HI'][:].reshape(-1)
    elevation['data'] = dvars['elevationR_HI'][:].reshape(-1)
    fixed_angle['data'] = dvars['elevationR_HI'][:].reshape(-1)

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=None)
