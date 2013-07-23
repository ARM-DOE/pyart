"""
pyart.io.nexrad

"""

import numpy as np

from .radar import Radar
from .common import get_metadata, make_time_unit_str
from .nexrad_level2_common import Archive2File


def read_nexrad_level2_archive(filename):
    """
    Read a NEXRAD Level 2 Archive file.

    Parameters
    ----------
    filename : str
        Filename of NEXRAD Level 2 Archive file.

    Supports NCDC and Motherload files.

    Returns
    -------
    radars : list
        List of radar objects.

    """
    afile = Archive2File(filename)
    scan_info = afile.get_scan_info()

    # super resolution reflectivity radar
    if len(scan_info['REF'][1]['scans']) != 0:
        moments = ['REF']
        scans = scan_info['REF'][1]['scans']
        max_gates = max(scan_info['REF'][1]['ngates'])
        refl_hi = _radar_from_archive2(afile, moments, scans, max_gates)
    else:
        refl_hi = None

    # standard resolution reflectivity radar
    if len(scan_info['REF'][2]['scans']) != 0:
        moments = ['REF']
        scans = scan_info['REF'][2]['scans']
        max_gates = max(scan_info['REF'][2]['ngates'])
        refl_std = _radar_from_archive2(afile, moments, scans, max_gates)
    else:
        refl_std = None

    # super resolution doppler radar
    if len(scan_info['VEL'][1]['scans']) != 0:
        moments = ['VEL', 'SW']
        if len(scan_info['ZDR'][1]['scans']) != 0:
            moments = ['VEL', 'SW', 'ZDR', 'PHI', 'RHO']
            dualpol_offset = (scan_info['ZDR'][1]['scans'][0] -
                              scan_info['VEL'][1]['scans'][0])
        else:
            dualpol_offset = 0
        scans = scan_info['VEL'][1]['scans']
        max_gates = max(scan_info['VEL'][1]['ngates'])
        doppler_hi = _radar_from_archive2(afile, moments, scans, max_gates,
                                          dualpol_offset)
    else:
        doppler_hi = None

    # standard resolution doppler radar
    if len(scan_info['VEL'][2]['scans']) != 0:
        moments = ['VEL', 'SW']
        if len(scan_info['ZDR'][2]['scans']) != 0:
            moments = ['VEL', 'SW', 'ZDR', 'PHI', 'RHO']
        scans = scan_info['VEL'][2]['scans']
        max_gates = max(scan_info['VEL'][2]['ngates'])
        doppler_std = _radar_from_archive2(afile, moments, scans, max_gates)
    else:
        doppler_std = None

    return refl_hi, doppler_hi, refl_std, doppler_std


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
        'valid_max': 126.0,          # or 63.0
        'valid_min': -127.0,         # or -63.6
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


def _radar_from_archive2(afile, moments, scans, max_gates, dualpol_offset=0):
    """
    Create a radar object from a Archive2File object.
    """

    # time
    time = get_metadata('time')
    time_start, _time = afile.get_time_scans(scans)
    time['data'] = _time
    time['units'] = make_time_unit_str(time_start)

    # _range
    _range = get_metadata('range')
    _range['data'] = afile.get_scan_range(scans[0], moments[0])

    # fields
    fields = {}
    for moment in moments:
        field_name = NEXRADFIELDS[moment]
        fields[field_name] = NEXRAD_METADATA[field_name].copy()
        fields[field_name] = get_metadata(field_name)
        if moment in ['ZDR', 'RHO', 'PHI']:
            offset_scans = [s + dualpol_offset for s in scans]
            fdata = afile.get_data(offset_scans, moment, max_gates)
        else:
            fdata = afile.get_data(scans, moment, max_gates)
        fields[field_name]['data'] = fdata

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

    nsweeps = len(scans)
    sweep_number['data'] = np.arange(len(scans), dtype='int32')
    sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])
    rays_per_scan = afile.get_nrays(scans[0])
    ssri = np.arange(nsweeps, dtype='int32') * rays_per_scan
    sweep_start_ray_index['data'] = ssri
    sweep_end_ray_index['data'] = ssri + (rays_per_scan - 1)

    # azimuth, elevation
    azimuth = get_metadata('azimuth')
    elevation = get_metadata('elevation')
    azimuth['data'] = afile.get_azimuth_angles_scans(scans)
    elev = afile.get_elevation_angles_scans(scans)
    elevation['data'] = elev.astype('float32')
    fixed_angle['data'] = elev.astype('float32')

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=None)
