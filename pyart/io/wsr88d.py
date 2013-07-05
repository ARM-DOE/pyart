"""

"""

import numpy as np

from .radar import Radar
from .common import get_metadata, make_time_unit_str
from .wsr88d_level2 import Archive2File


def read_wsr88d(filename):
    """
    """
    afile = Archive2File(filename)
    return [_radar_from_archive2(afile, i) for i in range(afile.nscans)]


NEXRADFIELDS = {
    'REF': 'reflectivity',
    'VEL': 'velocity',
    'SW': 'spectrum_width',
    'ZDR': 'differential_reflectivity',
    'PHI': 'differential_phase',
    'RHO': 'correlation_coefficient'
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


def _radar_from_archive2(afile, scan_num):
    """
    """
    # time
    time = get_metadata('time')
    time_start, _time = afile.get_scan_time(scan_num)
    time['data'] = _time
    time['units'] = make_time_unit_str(time_start)

    # _range
    _range = get_metadata('range')
    _range['data'] = afile.get_scan_range(scan_num, 'REF')

    # fields
    fields = {}
    for moment in afile.get_scan_moment_names(scan_num):
        field_name = NEXRADFIELDS[moment]
        fields[field_name] = NEXRAD_METADATA[field_name].copy()
        #fields[field_name] = get_metadata(field_name)
        fields[field_name]['data'] = afile.get_scan_moment(scan_num, moment)
    nrays = fields['reflectivity']['data'].shape[0]

    # metadata
    metadata = {'original_container': 'NEXRAD Level II Archive'}
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

    lat, lon, alt = afile.get_location()
    latitude['data'] = np.array([lat], dtype='float64')
    longitude['data'] = np.array([lon], dtype='float64')
    altitude['data'] = np.array([alt], dtype='float64')

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
    # sweep_end_ray_index
    sweep_number = get_metadata('sweep_number')
    sweep_mode = get_metadata('sweep_mode')
    fixed_angle = get_metadata('fixed_angle')
    sweep_start_ray_index = get_metadata('sweep_start_ray_index')
    sweep_end_ray_index = get_metadata('sweep_end_ray_index')

    sweep_number['data'] = np.arange(1, dtype='int32')
    sweep_mode['data'] = np.array(1 * ['azimuth_surveillance'])
    sweep_start_ray_index['data'] = np.array([0], dtype='int32')
    sweep_end_ray_index['data'] = np.array([nrays - 1], dtype='int32')

    # azimuth, elevation
    azimuth = get_metadata('azimuth')
    elevation = get_metadata('elevation')
    azimuth['data'] = np.array(afile.get_scan_azimuth_angles(scan_num))
    elev = np.array(afile.get_scan_elevation_angles(scan_num))
    elevation['data'] = np.array([elev.mean()])
    fixed_angle['data'] = np.array([elev.mean()])

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=None)
