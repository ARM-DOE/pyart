#! /usr/bin/env python
"""
Make a small netCDF CF/Radial file containing a single PPI scan.

Single field and scan is converted from sigmet file XSW110520105408.RAW7HHF
"""

import pyart

radar = pyart.io.read_rsl('XSW110520105408.RAW7HHF')

time_slice = slice(None, 400, 10)
range_slice = slice(None, None, 16)
sweep_slice = slice(None, 1)

# remove all but the reflectivity_horizontal fields
rf_field = radar.fields['reflectivity']
rf_data = rf_field['data']
rf_field['data'] = rf_data[time_slice, range_slice]
radar.fields = {'reflectivity_horizontal': rf_field}

radar.nsweeps = 1
radar.nray = 40
radar.ngates = 42

# truncate the range based variables
radar.range['data'] = radar.range['data'][range_slice]

# truncate the time based variables
radar.time['data'] = radar.time['data'][time_slice]
radar.azimuth['data'] = radar.azimuth['data'][time_slice]
radar.elevation['data'] = radar.elevation['data'][time_slice]
radar.instrument_parameters['prt']['data'] = \
    radar.instrument_parameters['prt']['data'][time_slice]

radar.instrument_parameters['unambiguous_range']['data'] = \
    radar.instrument_parameters['unambiguous_range']['data'][time_slice]

radar.instrument_parameters['nyquist_velocity']['data'] = \
    radar.instrument_parameters['nyquist_velocity']['data'][time_slice]

# truncate the sweep based variables
radar.sweep_number['data'] = radar.sweep_number['data'][sweep_slice]
radar.fixed_angle['data'] = radar.fixed_angle['data'][sweep_slice]
radar.sweep_start_ray_index['data'] = \
    radar.sweep_start_ray_index['data'][sweep_slice]
radar.sweep_end_ray_index['data'] = \
    radar.sweep_end_ray_index['data'][sweep_slice]
radar.sweep_end_ray_index['data'][0] = 39
radar.sweep_mode['data'] = radar.sweep_mode['data'][sweep_slice]

radar.sweep_number['data'] = radar.sweep_number['data'][sweep_slice]

radar.instrument_parameters['prt_mode']['data'] = \
    radar.instrument_parameters['prt_mode']['data'][sweep_slice]

# adjust metadata
radar.metadata = {
    'Conventions': 'CF/Radial instrument_parameters',
    'version': '1.2',
    'title': 'Py-ART Example PPI CF/Radial file',
    'institution': ('United States Department of Energy - Atmospheric '
                    'Radiation Measurement (ARM) program'),
    'references': 'none',
    'source': 'ARM SGP XSAPR Radar',
    'history': 'created by jhelmus on evs348532 at 2013-05-22T12:34:56',
    'comment': 'none',
    'instrument_name': 'xsapr-sgp'}

pyart.io.write_cfradial('example_cfradial_ppi.nc', radar)
