"""
Configuration file for the Python ARM Radar Toolkit (Py-ART)

The values for a number of Py-ART parameters and the default metadata created
when reading files, correcting fields, etc. is controlled by this single
Python configuration file.

Py-ART's configuration can be modified by setting the environment variable
PYART_CONFIG to point to a configuration file with formatting similar to this
file.

The recommended method for changing these defaults is for users to copy this
file into their home directory, rename it to .pyart_config.py, make any
desired changes, and adjust their login scripts to set the PYART_CONFIG
environment variable to point to .pyart_config.py in their home directory.

Py-ART's configuration can also be modified within a script or shell session
using the load_config functions, such modification will last until the end of
the script/session or until a new configuration is loaded.

"""

##############################################################################
##############################################################################
# Simple configuration
#
# Adjust the values of the variable (right hand side of the equal sign) in
# this section for an easy method of customizing Py-ART.  Do not change the
# variable names (the left hand side of the equal sign). More advanced
# settings are based upon these variables.  Most users will find that
# adjusting this section is all that is needed.
##############################################################################
##############################################################################

# The default fill value for masked arrays and _FillValue keys
fill_value = -9999.0

# Field names used when reading in radar files and in the various correction
# and retrieval algorithms. The comments in this section provide additional
# information about the fields in that section.

# Radar reflectivity fields, DZ
reflectivity = 'reflectivity'
corrected_reflectivity = 'corrected_reflectivity'
total_power = 'total_power'

# Mean Doppler velocity fields, VEL
velocity = 'velocity'
corrected_velocity = 'corrected_velocity'
simulated_velocity = 'simulated_velocity'

# Spectral width fields, SW
spectrum_width = 'spectrum_width'

# Differential reflectivity fields, ZDR
differential_reflectivity = 'differential_reflectivity'
corrected_differential_reflectivity = 'corrected_differential_reflectivity'

# Cross correlation ratio, correlation coefficient, RhoHV
cross_correlation_ratio = 'cross_correlation_ratio'

# Normalized coherent power, signal quality index, SQI, NCP
normalized_coherent_power = 'normalized_coherent_power'

# Differential phase shift, PhiDP
differential_phase = 'differential_phase'
unfolded_differential_phase = 'unfolded_differential_phase'
corrected_differential_phase = 'corrected_differential_phase'

# Specific differential phase shift, KDP
specific_differential_phase = 'specific_differential_phase'
corrected_specific_differential_phase = 'corrected_specific_differential_phase'

# Linear depolarization ration (h - horizontal, v - vertical), LDR
linear_depolarization_ratio = 'linear_polarization_ratio'
linear_depolarization_ratio_h = 'linear_depolarization_ratio_h'
linear_depolarization_ratio_v = 'linear_depolarization_ratio_v'

# Misc fields
signal_to_noise_ratio = 'signal_to_noise_ratio'
rain_rate = 'rain_rate'
radar_estimated_rain_rate = 'radar_estimated_rain_rate'
radar_echo_classification = 'radar_echo_classification'
specific_attenuation = 'specific_attenuation'
specific_differential_attenuation = 'specific_differential_attenuation'
clutter_filter_power_removed = 'clutter_filter_power_removed'

# Textures
differential_phase_texture = 'differential_phase_texture'

# Wind retrieval fields
eastward_wind_component = 'eastward_wind_component'
northward_wind_component = 'northward_wind_component'
vertical_wind_component = 'vertical_wind_component'

# profile variables
height = 'height'
interpolated_profile = 'interpolated_profile'
height_over_iso0 = 'height_over_iso0'

# path integrated attenuation
path_integrated_attenuation = 'path_integrated_attenuation'
path_integrated_differential_attenuation = ('path_integrated' +
                                            'differential_attenuation')

# End of Simple Configuration section

##############################################################################
##############################################################################
# Advanced Configuration
#
# Most users will not want to make any changes in this section.  For users
# who want a more fine-grained control over Py-ART's configuration this
# section provides access to these controls.  The layout of this section can
# be changed, the only requirement for a valid configuration file is that
# the ALL CAPITALIZED variable must must be present with the formatting
# present in this file.  These required variables are:
#
# FILL_VALUE, DEFAULT_METADATA, FILE_SPECIFIC_METADATA, FIELD_MAPPINGS,
# DEFAULT_FIELD_NAMES
#
# This section makes generous use of the variables in the Simple Configuration
# section, this is not required, but simplifies and enforces uniformity on
# the configuration.
##############################################################################
##############################################################################


##############################################################################
# Parameters
#
# Various parameters used in Py-ART.
##############################################################################

# the default fill value for masked arrays and the _FillValue key
FILL_VALUE = fill_value

# The DEFAULT_FIELD_NAMES controls the field names which are used in the
# correction and retrieval algorithms in Py-ART. The keys of the dictionary
# are "internal" names which cannot change, the values are the field names
# which will be used in the algorithms by default. For best results use the
# names defined by the variables in simple configuration section which are
# also used in the DEFAULT_METADATA and FIELD_MAPPINGS variable. If you
# choose to change a field name the names should also be changed in the
# DEFAULT_METADATA and FIELD_MAPPINGS variable. This is not required but
# highly suggested.

DEFAULT_FIELD_NAMES = {
    # Internal field name (do not change): field name used (can change)
    'reflectivity': reflectivity,
    'corrected_reflectivity': corrected_reflectivity,
    'total_power': total_power,
    'velocity': velocity,
    'corrected_velocity': corrected_velocity,
    'simulated_velocity': simulated_velocity,
    'spectrum_width': spectrum_width,
    'differential_reflectivity': differential_reflectivity,
    'corrected_differential_reflectivity':
    corrected_differential_reflectivity,
    'cross_correlation_ratio': cross_correlation_ratio,
    'normalized_coherent_power': normalized_coherent_power,
    'differential_phase': differential_phase,
    'unfolded_differential_phase': unfolded_differential_phase,
    'corrected_differential_phase': corrected_differential_phase,
    'specific_differential_phase': specific_differential_phase,
    'corrected_specific_differential_phase':
    corrected_specific_differential_phase,
    'linear_depolarization_ratio': linear_depolarization_ratio,
    'linear_depolarization_ratio_h': linear_depolarization_ratio_h,
    'linear_depolarization_ratio_v': linear_depolarization_ratio_v,
    'signal_to_noise_ratio': signal_to_noise_ratio,
    'rain_rate': rain_rate,
    'radar_estimated_rain_rate': radar_estimated_rain_rate,
    'radar_echo_classification': radar_echo_classification,
    'specific_attenuation': specific_attenuation,
    'differential_phase_texture': differential_phase_texture,
    'eastward_wind_component': eastward_wind_component,
    'northward_wind_component': northward_wind_component,
    'vertical_wind_component': vertical_wind_component,
    'height': height,
    'height_over_iso0': height_over_iso0,
    'interpolated_profile': interpolated_profile,
    'path_integrated_attenuation': path_integrated_attenuation,
    'specific_differential_attenuation': specific_differential_attenuation,
    'path_integrated_differential_attenuation':
    path_integrated_differential_attenuation,
    'clutter_filter_power_removed': clutter_filter_power_removed,
}


##############################################################################
# Default metadata
#
# The DEFAULT_METADATA dictionary contains dictionaries which provide the
# default radar attribute and field metadata. When reading in a file with
# Py-ART the FILE_SPECIFIC_METADATA variable is first queued for a metadata
# dictionary, if it is not found then the metadata in DEFAULT_METADATA is
# utilized.
##############################################################################

DEFAULT_METADATA = {

    # Metadata for radar attributes. These closely follow the CF/Radial
    # standard
    'azimuth': {
        'units': 'degrees',
        'standard_name': 'beam_azimuth_angle',
        'long_name': 'azimuth_angle_from_true_north',
        'axis': 'radial_azimuth_coordinate',
        'comment': 'Azimuth of antenna relative to true north'},

    'elevation': {
        'units': 'degrees',
        'standard_name': 'beam_elevation_angle',
        'long_name': 'elevation_angle_from_horizontal_plane',
        'axis': 'radial_elevation_coordinate',
        'comment': 'Elevation of antenna relative to the horizontal plane'},

    'scan_rate': {
        'units': 'degrees_per_second',
        'long_name': 'Antenna angle scan rate'},

    'range': {
        'units': 'meters',
        'standard_name': 'projection_range_coordinate',
        'long_name': 'range_to_measurement_volume',
        'axis': 'radial_range_coordinate',
        'spacing_is_constant': 'true',
        'comment': (
            'Coordinate variable for range. Range to center of each bin.')},

    'time': {
        'units': 'seconds',
        'standard_name': 'time',
        'long_name': 'time_in_seconds_since_volume_start',
        'calendar': 'gregorian',
        'comment': ('Coordinate variable for time. '
                    'Time at the center of each ray, in fractional seconds '
                    'since the global variable time_coverage_start')},

    'metadata': {
        'Conventions': 'CF/Radial instrument_parameters',
        'version': '1.3',
        'title': '',
        'institution': '',
        'references': '',
        'source': '',
        'history': '',
        'comment': '',
        'instrument_name': ''},


    # Metadata for radar sweep information dictionaries
    'sweep_number': {
        'units': 'count',
        'standard_name': 'sweep_number',
        'long_name': 'Sweep number'},

    'sweep_mode': {
        'units': 'unitless',
        'standard_name': 'sweep_mode',
        'long_name': 'Sweep mode',
        'comment': ('Options are: "sector", "coplane", "rhi", '
                    '"vertical_pointing", "idle", "azimuth_surveillance", '
                    '"elevation_surveillance", "sunscan", "pointing", '
                    '"manual_ppi", "manual_rhi"')},

    'fixed_angle': {
        'long_name': 'Target angle for sweep',
        'units': 'degrees',
        'standard_name': 'target_fixed_angle'},

    'sweep_start_ray_index': {
        'long_name': 'Index of first ray in sweep, 0-based',
        'units': 'count'},

    'sweep_end_ray_index': {
        'long_name': 'Index of last ray in sweep, 0-based',
        'units': 'count'},

    'rays_per_sweep': {
        'long_name': 'Number of rays in each sweep',
        'units': 'count'},

    'target_scan_rate': {
        'long_name': 'Target scan rate for sweep',
        'units': 'degrees_per_second',
        },

    'rays_are_indexed': {
        'long_name': 'Flag for indexed rays',
        'units': 'unitless',
        'options': ('true: rays are indexed to a regular grid, ' +
                    'false: rays are not indexed to a regular grid'),
        },

    'ray_angle_res': {
        'long_name': 'Angular resolution between rays',
        'units': 'degrees',
        'comment': 'Only applicable when rays_are_indexed variable is true',
        },

    # Metadata for radar location attributes
    'latitude': {
        'long_name': 'Latitude',
        'standard_name': 'Latitude',
        'units': 'degrees_north'},

    'longitude': {
        'long_name': 'Longitude',
        'standard_name': 'Longitude',
        'units': 'degrees_east'},

    'altitude': {
        'long_name': 'Altitude',
        'standard_name': 'Altitude',
        'units': 'meters',
        'positive': 'up'},

    'gate_x': {
        'long_name': 'Cartesian x location of gate with origin at the radar',
        'units': 'meters'},

    'gate_y': {
        'long_name': 'Cartesian y location of gate with origin at the radar',
        'units': 'meters'},

    'gate_z': {
        'long_name': 'Cartesian z location of gate with origin at the radar',
        'units': 'meters'},

    'gate_edge_x': {
        'long_name': 'Cartesian x location of the edges of each gate',
        'units': 'meters'},

    'gate_edge_y': {
        'long_name': 'Cartesian y location of the edges of each gate',
        'units': 'meters'},

    'gate_edge_z': {
        'long_name': 'Cartesian z location of the edges of each gate',
        'units': 'meters'},

    'gate_longitude': {
        'long_name': 'Longitude of radar gate.',
        'units': 'degrees_north'},

    'gate_latitude': {
        'long_name': 'Latitude of radar gate',
        'units': 'degrees_east'},

    'gate_altitude': {
        'long_name': 'Altitude of radar gate',
        'units': 'meters'},

    # Metadata for instrument_parameter dictionary
    'prt_mode': {
        'comments': ('Pulsing mode Options are: "fixed", "staggered", '
                     '"dual". Assumed "fixed" if missing.'),
        'meta_group': 'instrument_parameters',
        'long_name': 'Pulsing mode',
        'units': 'unitless'},

    'nyquist_velocity': {
        'units': 'meters_per_second',
        'comments': "Unambiguous velocity",
        'meta_group': 'instrument_parameters',
        'long_name': 'Nyquist velocity'},

    'prt': {
        'units': 'seconds',
        'comments': ("Pulse repetition time. For staggered prt, "
                     "also see prt_ratio."),
        'meta_group': 'instrument_parameters',
        'long_name': 'Pulse repetition time'},

    'unambiguous_range': {
        'units': 'meters',
        'comments': 'Unambiguous range',
        'meta_group': 'instrument_parameters',
        'long_name': 'Unambiguous range'},

    'pulse_width': {
        'units': 'seconds',
        'comments': 'Pulse width',
        'meta_group': 'instrument_parameters',
        'long_name': 'Pulse width'},

    'prt_ratio': {
        'units': 'unitless',
        'meta_group': 'instrument_parameters',
        'long_name': 'Pulse repetition frequency ratio'},

    'frequency': {
        'units': 's-1',
        'meta_group': 'instrument_parameters',
        'long_name': 'Radiation frequency'},

    'n_samples': {
        'units': 'unitless',
        'meta_group': 'instrument_parameters',
        'long_name': 'Number of samples used to compute moments'},

    # non-standard parameter for specifying the PRF high/low for each ray
    'prf_flag': {
        'units': 'unitless',
        'comments': "PRF used to collect ray. 0 for high PRF, 1 for low PRF.",
        'meta_group': 'instrument_parameters',
        'long_name': 'PRF flag'},

    # Metadata for radar_parameter sub-convention
    'radar_beam_width_h': {
        'units': 'degrees',
        'meta_group': 'radar_parameters',
        'long_name': 'Antenna beam width H polarization'},

    'radar_beam_width_v': {
        'units': 'degrees',
        'meta_group': 'radar_parameters',
        'long_name': 'Antenna beam width V polarization'},

    # Metadata for airborne radar parameters
    'rotation': {
        'units': 'degrees',
        'standard_name': 'ray_rotation_angle_relative_to_platform',
        'long_name': 'Ray rotation angle relative to platform'},

    'tilt': {
        'units': 'degrees',
        'standard_name': 'ray_tilt_angle_relative_to_platform',
        'long_name': 'Ray tilt angle relative to platform'},

    'roll': {
        'units': 'degrees',
        'standard_name': 'platform_roll_angle',
        'long_name': 'Platform roll angle'},

    'drift': {
        'units': 'degrees',
        'standard_name': 'platform_drift_angle',
        'long_name': 'Platform drift angle'},

    'heading': {
        'units': 'degrees',
        'standard_name': 'platform_heading_angle',
        'long_name': 'Platform heading angle'},

    'pitch': {
        'units': 'degrees',
        'standard_name': 'platform_pitch_angle',
        'long_name': 'Platform pitch angle'},

    'georefs_applied': {
        'units': 'unitless',
        'standard_name': 'georefs_have_been_applied_to_ray',
        'long_name': 'Geoferences have been applied to ray',
        'comment': '1 if georefs have been applied, 0 otherwise'},

    # Reflectivity fields
    reflectivity: {
        'units': 'dBZ',
        'standard_name': 'equivalent_reflectivity_factor',
        'long_name': 'Reflectivity',
        'coordinates': 'elevation azimuth range'},

    corrected_reflectivity: {
        'units': 'dBZ',
        'standard_name': 'corrected_equivalent_reflectivity_factor',
        'long_name': 'Corrected reflectivity',
        'coordinates': 'elevation azimuth range'},

    total_power: {
        'units': 'dBZ',
        'standard_name': 'equivalent_reflectivity_factor',
        'long_name': 'Total power',
        'coordinates': 'elevation azimuth range'},

    # Velocity fields
    velocity: {
        'units': 'meters_per_second',
        'standard_name': 'radial_velocity_of_scatterers_away_from_instrument',
        'long_name': 'Mean dopper velocity',
        'coordinates': 'elevation azimuth range'},

    corrected_velocity: {
        'units': 'meters_per_second',
        'standard_name': (
            'corrected_radial_velocity_of_scatterers_away_from_instrument'),
        'long_name': 'Corrected mean doppler velocity',
        'coordinates': 'elevation azimuth range'},

    simulated_velocity: {
        'units': 'meters_per_second',
        'standard_name': (
            'radial_velocity_of_scatterers_away_from_instrument'),
        'long_name': 'Simulated mean doppler velocity',
        'coordinates': 'elevation azimuth range'},

    # Spectrum width fields
    spectrum_width: {
        'units': 'meters_per_second',
        'standard_name': 'doppler_spectrum_width',
        'long_name': 'Doppler spectrum width',
        'coordinates': 'elevation azimuth range'},

    # Dual-polarization fields
    differential_reflectivity: {
        'units': 'dB',
        'standard_name': 'log_differential_reflectivity_hv',
        'long_name': 'Differential reflectivity',
        'coordinates': 'elevation azimuth range'},

    corrected_differential_reflectivity: {
        'units': 'dB',
        'standard_name': 'corrected_log_differential_reflectivity_hv',
        'long_name': 'Corrected differential reflectivity',
        'coordinates': 'elevation azimuth range'},

    cross_correlation_ratio: {
        'units': 'ratio',
        'standard_name': 'cross_correlation_ratio_hv',
        'long_name': 'Cross correlation ratio (RHOHV)',
        'valid_max': 1.0,
        'valid_min': 0.0,
        'coordinates': 'elevation azimuth range'},

    normalized_coherent_power: {
        'units': 'ratio',
        'standard_name': 'normalized_coherent_power',
        'long_name': 'Normalized coherent power',
        'valid_max': 1.0,
        'valid_min': 0.0,
        'comment': 'Also know as signal quality index (SQI)',
        'coordinates': 'elevation azimuth range'},

    differential_phase: {
        'units': 'degrees',
        'standard_name': 'differential_phase_hv',
        'long_name': 'Differential phase (PhiDP)',
        'valid_max': 180.0,
        'valid_min': -180.0,
        'coordinates': 'elevation azimuth range'},

    unfolded_differential_phase: {
        'units': 'degrees',
        'standard_name': 'differential_phase_hv',
        'long_name': 'Unfolded differential phase',
        'coordinates': 'elevation azimuth range'},

    corrected_differential_phase: {
        'units': 'degrees',
        'standard_name': 'differential_phase_hv',
        'long_name': 'Corrected differential phase',
        'coordinates': 'elevation azimuth range'},

    specific_differential_phase: {
        'units': 'degrees/km',
        'standard_name': 'specific_differential_phase_hv',
        'long_name': 'Specific differential phase (KDP)',
        'coordinates': 'elevation azimuth range'},

    corrected_specific_differential_phase: {
        'units': 'degrees/km',
        'standard_name': 'specific_differential_phase_hv',
        'long_name': 'Corrected specific differential phase (KDP)',
        'coordinates': 'elevation azimuth range'},


    # Depolarization ratio fields
    linear_depolarization_ratio: {
        'units': 'dB',
        'standard_name': 'log_linear_depolarization_ratio_hv',
        'long_name': 'Linear depolarization ratio',
        'coordinates': 'elevation azimuth range'},

    linear_depolarization_ratio_h: {
        'units': 'dB',
        'standard_name': 'log_linear_depolarization_ratio_h',
        'long_name': 'Linear depolarization ratio horizontal',
        'coordinates': 'elevation azimuth range'},

    linear_depolarization_ratio_v: {
        'units': 'dB',
        'standard_name': 'log_linear_depolarization_ratio_v',
        'long_name': 'Linear depolarization ratio vertical',
        'coordinates': 'elevation azimuth range'},

    # Misc fields
    signal_to_noise_ratio: {
        'units': 'dB',
        'standard_name': 'signal_to_noise_ratio',
        'long_name': 'Signal to noise ratio',
        'coordinates': 'elevation azimuth range'},

    rain_rate: {
        'units': 'kg/m2/s',
        'standard_name': 'rain_rate',
        'long_name': 'Rain rate',
        'coordinates': 'elevation azimuth range'},

    radar_estimated_rain_rate: {
        'units': 'mm/hr',
        'standard_name': 'radar_estimated_rain_rate',
        'long_name': 'Radar estimated rain rate',
        'coordinates': 'elevation azimuth range'},

    radar_echo_classification: {
        'units': 'legend',
        'standard_name': 'radar_echo_classification',
        'long_name': 'Radar Echo classification',
        'coordinates': 'elevation azimuth range'},

    specific_attenuation: {
        'units': 'dB/km',
        'standard_name': 'specific_attenuation',
        'long_name': 'Specific attenuation',
        'valid_min': 0.0,
        'valid_max': 1.0,
        'coordinates': 'elevation azimuth range'},

    specific_differential_attenuation: {
        'units': 'dB/km',
        'long_name': 'Specific Differential Attenuation',
        'coordinates': 'elevation azimuth range'},

    path_integrated_attenuation: {
        'units': 'dB',
        'long_name': 'Path Integrated Attenuation',
        'coordinates': 'elevation azimuth range'},

    path_integrated_differential_attenuation: {
        'units': 'dB',
        'long_name': 'Path Integrated Differential Attenuation',
        'coordinates': 'elevation azimuth range'},

    # Textures
    differential_phase_texture: {
        'units': 'degrees',
        'standard_name': 'differential_phase_hv_texture',
        'long_name': 'Texture of differential phase (PhiDP)',
        'coordinates': 'elevation azimuth range'},

    # Wind retrieval fields
    eastward_wind_component: {
        'units': 'meters_per_second',
        'standard_name': 'eastward_wind_component',
        'long_name': 'Eastward wind component'},

    northward_wind_component: {
        'units': 'meters_per_second',
        'standard_name': 'northward_wind_component',
        'long_name': 'Northward wind component'},

    vertical_wind_component: {
        'units': 'meters_per_second',
        'standard_name': 'vertical_wind_component',
        'long_name': 'Vertical wind component'},

    # profile variables
    height: {
        'units': 'm',
        'long_name': 'Height of radar beam',
        'standard_name': 'height'},

    height_over_iso0: {
        'units': 'm',
        'long_name': 'Height of radar beam'},

    interpolated_profile: {
        'long_name': 'Interpolated profile',
        'standard_name':  'interpolated_profile',
        'units': 'unknown'},

    clutter_filter_power_removed: {
        'units': 'dB',
        'long_name': 'Clutter filter power removed',
        'standard_name': 'clutter_filter_power_removed',
        'coordinates': 'elevation azimuth range'},

    # Grid metadata

    'grid_time': {
        'units': 'seconds',
        'standard_name': 'time',
        'long_name': 'Time of grid',
        'calendar': 'gregorian'},

    'origin_longitude': {
        'long_name': 'Longitude at grid origin',
        'units': 'degrees_east',
        'standard_name': 'longitude',
        'valid_min': -180.,
        'valid_max': 180.},

    'origin_latitude': {
        'long_name': 'Latitude at grid origin',
        'units': 'degrees_north',
        'standard_name': 'latitude',
        'valid_min': -90.,
        'valid_max': 90.},

    'origin_altitude': {
        'long_name': 'Altitude at grid origin',
        'units': 'm',
        'standard_name': 'altitude'},

    'x': {
        'standard_name': 'projection_x_coordinate',
        'long_name': 'X distance on the projection plane from the origin',
        'axis': 'X',
        'units': 'm'},

    'y': {
        'standard_name': 'projection_y_coordinate',
        'long_name': 'Y distance on the projection plane from the origin',
        'axis': 'Y',
        'units': 'm'},

    'z': {
        'standard_name': 'projection_z_coordinate',
        'long_name': 'Z distance on the projection plane from the origin',
        'axis': 'Z',
        'units': 'm',
        'positive': 'up'},

    'point_x': {
        'long_name': 'Cartesian x distance of each grid point from the origin',
        'units': 'meters'},

    'point_y': {
        'long_name': 'Cartesian y distance of each grid point from the origin',
        'units': 'meters'},

    'point_z': {
        'long_name': 'Cartesian z distance of each grid point from the origin',
        'positive': 'up',
        'units': 'meters'},

    'point_longitude': {
        'long_name': 'Longitude of each grid point',
        'units': 'degrees_north'},

    'point_latitude': {
        'long_name': 'Latitude of each grid point',
        'units': 'degrees_east'},

    'point_altitude': {
        'long_name': 'Altitude of each grid point',
        'units': 'meters'},

    'radar_latitude': {
        'long_name': 'Latitude of radars used to make the grid.',
        'units': 'degrees_north', },

    'radar_longitude': {
        'long_name': 'Longitude of radars used to make the grid.',
        'units': 'degrees_east', },

    'radar_altitude': {
        'long_name': 'Altitude of radars used to make the grid.',
        'units': 'm', },

    'radar_time': {
        'calendar': 'gregorian',
        'long_name': 'Time in seconds of the volume start for each radar'},

    'radar_name': {
        'long_name': 'Name of radar used to make the grid', },

}


##############################################################################
# File specific metadata
#
# These dictionaries define metadata that is to be used only when reading in
# a given type of file.  This metadata is used in place of the
# DEFAULT_METADATA when it is avialable.  The main use of these variable
# is to define field specific data, it is safe to leave some/all of these
# empty if the default metadata is acceptable.
##############################################################################

# Metadata for Sigmet/IRIS files
sigmet_metadata = {}

# Metadata for NEXRAD Level II files (Archive and CDM files)
nexrad_metadata = {
    reflectivity: {
        'units': 'dBZ',
        'standard_name': 'equivalent_reflectivity_factor',
        'long_name': 'Reflectivity',
        'valid_max': 94.5,
        'valid_min': -32.0,
        'coordinates': 'elevation azimuth range'},

    velocity: {
        'units': 'meters_per_second',
        'standard_name': 'radial_velocity_of_scatterers_away_from_instrument',
        'long_name': 'Mean doppler Velocity',
        'valid_max': 95.0,
        'valid_min': -95.0,
        'coordinates': 'elevation azimuth range'},

    spectrum_width: {
        'units': 'meters_per_second',
        'standard_name': 'doppler_spectrum_width',
        'long_name': 'Spectrum Width',
        'valid_max': 63.0,
        'valid_min': -63.5,
        'coordinates': 'elevation azimuth range'},

    differential_reflectivity: {
        'units': 'dB',
        'standard_name': 'log_differential_reflectivity_hv',
        'long_name': 'log_differential_reflectivity_hv',
        'valid_max': 7.9375,
        'valid_min': -7.8750,
        'coordinates': 'elevation azimuth range'},

    differential_phase: {
        'units': 'degrees',
        'standard_name': 'differential_phase_hv',
        'long_name': 'differential_phase_hv',
        'valid_max': 360.0,
        'valid_min': 0.0,
        'coordinates': 'elevation azimuth range'},

    cross_correlation_ratio: {
        'units': 'ratio',
        'standard_name': 'cross_correlation_ratio_hv',
        'long_name': 'Cross correlation_ratio (RHOHV)',
        'valid_max': 1.0,
        'valid_min': 0.0,
        'coordinates': 'elevation azimuth range'},

    clutter_filter_power_removed: {
        'units': 'dB',
        'long_name': 'Clutter filter power removed',
        'standard_name': 'clutter_filter_power_removed',
        'valid_max': 73.0,
        'valid_min': 0.0,
        'coordinates': 'elevation azimuth range'},
}

# Metadata for NEXRAD Level 3 Products
nexrad_level3_metadata = {

    radar_estimated_rain_rate: {
        'units': 'inches',
        'standard_name': 'radar_estimated_rain_rate',
        'long_name': 'Radar estimated rain rate',
        'coordinates': 'elevation azimuth range'},

    radar_echo_classification: {
        'units': 'legend',
        'standard_name': 'radar_echo_classification',
        'long_name': 'Radar echo classification',
        'options': ('0: Below Threshold (ND), '
                    '10: Biological (BI), '
                    '20: Anomalous Propagation/Group Clutter (GC), '
                    '30: Ice Crystals (IC), '
                    '40: Dry Snow (DS), '
                    '50: Wet Snow (WS), '
                    '60: Light and/or Moderate Rain (RA), '
                    '70: Heavy Rain (HR), '
                    '80: Big Drops (rain) (BD), '
                    '90: Graupel (GR), '
                    '100: Hail, possibly with rain (HA), '
                    '110: Large Hail (LH), '
                    '120: Giant Hail (GH), '
                    '140: Unknown Classification (UK), '
                    '150: Range Folded (RH)'),
        'coordinates': 'elevation azimuth range'},
}

# Metadata for CF/Radial files
cfradial_metadata = {}

# Metadata for MDV files
mdv_metadata = {}

# Metadata for RSL files
rsl_metadata = {}

# Metadata for CSU-CHILL, CHL files
chl_metadata = {}

FILE_SPECIFIC_METADATA = {      # Required
    'sigmet': sigmet_metadata,
    'nexrad_archive': nexrad_metadata,
    'nexrad_cdm': nexrad_metadata,
    'nexrad_level3': nexrad_level3_metadata,
    'cfradial': cfradial_metadata,
    'mdv': mdv_metadata,
    'rsl': rsl_metadata,
    'chl': chl_metadata,
}

##############################################################################
# Field name mapping
#
# These dictionaries map file field names or data types to a radar field
# name.  These are used to populate the radar.fields dictionary during a read
# in Py-ART.  A value of None will not include that field in the radar object.
# These can be over-ridden on a per-read basis using the field_mapping
# parameter, or using setting the file_field_names parameter to True.
##############################################################################

# Sigmet/IRIS file field mapping
# Note that multiple sigmet fields map to the same radar field, if
# more than one of these fields are present the radar field will be
# overwritten with the last sigmet field.
sigmet_field_mapping = {
    # Sigmet data type :field name              # (Data_type) Description
    'XHDR': None,                               # (0) Extended Header
    'DBT': total_power,                         # (1) Total Power
    'DBZ': reflectivity,                        # (2) Reflectivity
    'VEL': velocity,                            # (3) Velocity
    'WIDTH': spectrum_width,                    # (4) Width
    'ZDR': differential_reflectivity,           # (5) Diff. reflectivity
    'DBZC': corrected_reflectivity,             # (7) Corrected reflectivity
    'DBT2': total_power,                        # (8) Total Power
    'DBZ2': reflectivity,                       # (9) Reflectivity
    'VEL2': velocity,                           # (10) Velocity
    'WIDTH2': spectrum_width,                   # (11) Width
    'ZDR2': differential_reflectivity,          # (12) Diff. reflectivity
    'RAINRATE2': radar_estimated_rain_rate,     # (13) Rainfall rate
    'KDP': specific_differential_phase,         # (14) KDP (diff. phase)
    'KDP2': specific_differential_phase,        # (15) KDP (diff. phase)
    'PHIDP': differential_phase,                # (16) PhiDP (diff. phase)
    'VELC': corrected_velocity,                 # (17) Corrected velocity
    'SQI': normalized_coherent_power,           # (18) SQI
    'RHOHV': cross_correlation_ratio,           # (19) RhoHV
    'RHOHV2': cross_correlation_ratio,          # (20) RhoHV
    'DBZC2': corrected_reflectivity,            # (21) Corrected Reflec.
    'VELC2': corrected_velocity,                # (21) Corrected Velocity
    'SQI2': normalized_coherent_power,          # (23) SQI
    'PHIDP2': differential_phase,               # (24) PhiDP (diff. phase)
    'LDRH': linear_depolarization_ratio_h,      # (25) LDR xmt H, rcv V
    'LDRH2': linear_depolarization_ratio_h,     # (26) LDR xmt H, rcv V
    'LDRV': linear_depolarization_ratio_v,      # (27) LDR xmt V, rcv H
    'LDRV2': linear_depolarization_ratio_v,     # (28) LDR xmt V, rcv H
    'HEIGHT': None,                             # (32) Height (1/10 km)
    'VIL2': None,                               # (33) Linear Liquid
    'RAW': None,                                # (34) Raw Data
    'SHEAR': None,                              # (35) Wind Shear
    'DIVERGE2': None,                           # (36) Divergence
    'FLIQUID2': None,                           # (37) Floated liquid
    'USER': None,                               # (38) User type
    'OTHER': None,                              # (39) Unspecified
    'DEFORM2': None,                            # (40) Deformation
    'VVEL2': None,                              # (41) Vertical velocity
    'HVEL2': None,                              # (42) Horizontal velocity
    'HDIR2': None,                              # (43) Horiz. wind direction
    'AXDIL2': None,                             # (44) Axis of dilation
    'TIME2': None,                              # (45) Time in seconds
    'RHOH': None,                               # (46) Rho, xmt H, rcv V
    'RHOH2': None,                              # (47) Rho, xmt H, rcv V
    'RHOV': None,                               # (48) Rho, xmt V, rcv H
    'RHOV2': None,                              # (49) Rho, xmt V, rcv H
    'PHIH': None,                               # (50) Phi, xmt H, rcv V
    'PHIH2': None,                              # (51) Phi, xmt H, rcv V
    'PHIV': None,                               # (52) Phi, xmt V, rcv H
    'PHIV2': None,                              # (53) Phi, xmt V, rcv H
    'USER2': None,                              # (54) User type
    'HCLASS': radar_echo_classification,        # (55) Hydrometeor class
    'HCLASS2': radar_echo_classification,       # (56) Hydrometeor class
    'ZDRC': corrected_differential_reflectivity,
                                                # (57) Corrected diff. refl.
    'ZDRC2': corrected_differential_reflectivity,
                                                # (58) Corrected diff. refl.
    'UNKNOWN_59': None,                         # Unknown field
    'UNKNOWN_60': None,                         # Unknown field
    'UNKNOWN_61': None,                         # Unknown field
    'UNKNOWN_62': None,                         # Unknown field
    'UNKNOWN_63': None,                         # Unknown field
    'UNKNOWN_64': None,                         # Unknown field
    'UNKNOWN_65': None,                         # Unknown field
    'UNKNOWN_66': None,                         # Unknown field
    # there may be more field, add as needed
}


# NEXRAD Level II Archive files
nexrad_archive_field_mapping = {
    # NEXRAD field: radar field name
    'REF': reflectivity,
    'VEL': velocity,
    'SW': spectrum_width,
    'ZDR': differential_reflectivity,
    'PHI': differential_phase,
    'RHO': cross_correlation_ratio,
    'CFP': clutter_filter_power_removed
}

# NEXRAD Level II CDM files
nexrad_cdm_field_mapping = {
    # CDM variable name (without _HI): radar field name
    'Reflectivity': reflectivity,
    'RadialVelocity': velocity,
    'SpectrumWidth': spectrum_width,
    'DifferentialReflectivity': differential_reflectivity,
    'DifferentialPhase': differential_phase,
    'CorrelationCoefficient': cross_correlation_ratio
}

# NEXRAD Level 3 Product files.
nexrad_level3_mapping = {
    # Message code : field name         # Product name
    19: reflectivity,                   # Base Reflectivity
    20: reflectivity,                   # Base Reflectivity
    25: velocity,                       # Base Velocity
    27: velocity,                       # Base Velocity
    28: spectrum_width,                 # Base Spectrum Width
    30: spectrum_width,                 # Base Spectrum Width
    32: reflectivity,                   # Digital Hybrid Scan Reflectivity
    34: None,                           # Clutter Filter Control
    56: velocity,                       # Storm Relative Mean Radial Velocity
    78: radar_estimated_rain_rate,      # Surface Rainfall Accum. (1 hr)
    79: radar_estimated_rain_rate,      # Surface Rainfall Accum. (3 hr)
    80: radar_estimated_rain_rate,      # Storm Total Rainfall Accumulation
    94: reflectivity,                   # Base Reflectivity Data Array
    99: velocity,                       # Base Velocity Data Array
    134: None,                          # High Resolution VIL
    135: None,                          # Enhanced Echo Tops
    138: radar_estimated_rain_rate,     # Digital Storm Total Precipitation
    159: differential_reflectivity,     # Digital Differential Reflectivity
    161: cross_correlation_ratio,       # Digital Correlation Coefficient
    163: specific_differential_phase,   # Digital Specific Differential Phase
    165: radar_echo_classification,     # Digital Hydrometeor Classification
    169: radar_estimated_rain_rate,     # One Hour Accumulation
    170: radar_estimated_rain_rate,     # Digital Accumulation Array
    171: radar_estimated_rain_rate,     # Storm Total Accumulation
    172: radar_estimated_rain_rate,     # Digital Storm Total Accumulation
    173: radar_estimated_rain_rate,     # Digital User-Selectable Accum.
    174: radar_estimated_rain_rate,     # Digital 1 hr Diff. Accum.
    175: radar_estimated_rain_rate,     # Digital Storm Total Diff. Accum.
    176: radar_estimated_rain_rate,     # Digital Inst. Precipitation Rate
    177: radar_echo_classification,     # Hybrid Hydrometeor Classification
    181: reflectivity,                  # Base Reflectivity
    182: velocity,                      # Base Velocity
    186: reflectivity,                  # Base Reflectivity
}

# MDV files
mdv_field_mapping = {
    # MDV moment: radar field name
    'DBZ_F': reflectivity,
    'VEL_F': velocity,
    'WIDTH_F': spectrum_width,
    'ZDR_F': differential_reflectivity,
    'RHOHV_F': cross_correlation_ratio,
    'NCP_F': normalized_coherent_power,
    'KDP_F': specific_differential_phase,
    'PHIDP_F': differential_phase,
    'VEL_COR': corrected_velocity,
    'PHIDP_UNF': unfolded_differential_phase,
    'KDP_SOB': corrected_specific_differential_phase,
    'DBZ_AC': corrected_reflectivity,
    # repeated integer moments
    'DBZ': reflectivity,
    'VEL': velocity,
    'WIDTH': spectrum_width,
    'ZDR': differential_reflectivity,
    'RHOHV': cross_correlation_ratio,
    'NCP': normalized_coherent_power,
    'KDP': specific_differential_phase,
    'PHIDP': differential_phase,
}

# CF/Radial files
cfradial_field_mapping = {}

# RSL files
# Note that multiple RSL field map to the same radar field, if
# more than one of these fields are present in the RSL data structure
# the radar field will be overwritten with the last field.
rsl_field_mapping = {
    # RSL 2 letter field: radar field           # RSL description
    'DZ': reflectivity,                         # reflectivity
    'VR': velocity,                             # velocity
    'SW': spectrum_width,                       # spectrum width
    'CZ': corrected_reflectivity,               # corrected reflectivity
    'ZT': reflectivity,                         # uncorrected reflectivity
    'DR': differential_reflectivity,            # differential reflectivity
    'LR': differential_reflectivity,            # another diff. reflectivity
    'ZD': differential_reflectivity,            # another diff. reflectivity
    'DM': None,                                 # received power
    'RH': cross_correlation_ratio,              # RhoHV
    'PH': differential_phase,                   # PhiDP
    'XZ': None,                                 # X-band reflectivity
    'CD': corrected_differential_reflectivity,  # Corrected DR.
    'MZ': None,                                 # DZ mask
    'MD': None,                                 # DR Mask
    'ZE': corrected_reflectivity,               # edited reflectivity
    'VE': corrected_velocity,                   # edited velocity
    'KD': specific_differential_phase,          # specific diff. phase
    'TI': None,                                 # TIME (unknown)
    'DX': None,                                 # ???
    'CH': None,                                 # ???
    'AH': None,                                 # ???
    'CV': None,                                 # ???
    'AV': None,                                 # ???
    'SQ': normalized_coherent_power,            # Signal Quality Index (sigmet)
    'VS': None,                                 # Radial Vel. combined
    'VL': None,                                 # Radial Vel. combined
    'VG': None,                                 # Radial Vel. combined
    'VT': None,                                 # Radial Vel. combined
    'NP': normalized_coherent_power,            # Normalized Coherent Power
    'HC': radar_echo_classification,            # Hydroclass
    'VC': None,                                 # Radial Vel. Corrected.
    'V2': None,                                 # Radial Vel cut 2
    'S2': None,                                 # Spectrum width cut 2
    'V3': None,                                 # Radial Vel cut 3
    'S3': None,                                 # Spectrum width cut 3
}

chl_field_mapping = {
    # Chill field name : radar field name
    'Z': reflectivity,
    'V': velocity,
    'W': spectrum_width,
    'ZDR': differential_reflectivity,
    'LDRH': linear_depolarization_ratio_h,
    'LDRV': linear_depolarization_ratio_v,
    b'\xce\xa8 DP'.decode('utf-8'): differential_phase,
    'KDP': specific_differential_phase,
    b'\xcf\x81 HV'.decode('utf-8'): cross_correlation_ratio,
    'NCP': normalized_coherent_power,
    # These fields are not mapped by default
    'H Re(lag 1)': None,    # Real part of lag-1 correlation, H Channel
    'V Re(lag 2)': None,    # Real part of lag-2 correlation, V Channel
    'VAvgQ': None,          # Average Q, V Channel
    'V Im(lag 1)': None,    # Imaginary part of lag-1 correlation, V Channel
    'HAvgQ': None,          # Average Q, H Channel
    'H Im(lag 2)': None,    # Imaginary part of lag-2 correlation, H Channel
    'V lag 0': None,        # Absolute value of lag-0 correlation, V Channel
    'H lag 0': None,        # Absolute value of lag-0 correlation, H Channel
    'H lag 0 cx': None,     # Absolute value of lag-0 cross correlation,
                            # H Channel
    'H Im(lag 1)': None,    # Imaginary part of lag-1 correlation, H Channel
    'H Re(lag 2)': None,    # Real part of lag-2 correlation, H Channel
    'V lag 0 cx': None,     # Absolute value of lag-0 cross correlation,
                            # V Channel
    'V Re(lag 1)': None,    # Real part of lag-1 correlation, V Channel
    'V Im(lag 2)': None,    # Imaginary part of lag-2 correlation, V Channel
    'HV lag 0 I': None,     # Real part of cross channel correlation at lag-0
    'HV lag 0 Q': None,     # Imaginary part of cross channel correlation at
                            # lag-0
    'VAvgI': None,          # Average I, V Channel
    'HAvgI': None,          # Average I, H Channel
    b'\xcf\x81 HCX'.decode('utf-8'): None,   # H Co to Cross Correlation
    b'\xcf\x81 VCX'.decode('utf-8'): None,   # V Co to Cross Correlation
}

# GAMIC HDF5 files
gamic_field_mapping = {
    # Description of GAMIC fields taken from Radx source code file:
    # GamicHdf5RadxFile.cc

    # General notes on field names
    # ----------------------------
    # * An 'U' prefix indicated the moment was corrected from a timeseries
    #   that was not clutter corrected.  Fields without a 'U' prefix are
    #   calculated from a timeseries that has been clutter corrected.
    # * An 'A' prefix indicates that the moment has been corrected
    #   for rainfall attenuation.
    # * 'h' and 'v' endings indicate a the horizontal and vertical channels

    # GAMIC field name: radar field name
    'I': None,          # In phase signal
    'Q': None,          # Quadrature signal

    'Z': corrected_reflectivity,
                        # Reflectivity, corrected for 2nd trip & clutter
    'UZ': reflectivity,
                        # Uncorrected reflectivity
    'Zh': corrected_reflectivity,
    'ZH': corrected_reflectivity,
                        # Corrected reflectivity, horizontal channel
    'Zv': None,
    'ZV': None,         # Corrected reflectivity, vertical channel
    'UZh': reflectivity,
    'UH': reflectivity,
                        # Uncorrected reflectivity, horizontal channel
    'UZv': None,
    'UV': None,         # Uncorrected reflectivity, vertical channel
    'AZh': None,        # Refl., rainfall atten. & clutter corrected, h chan.

    # 'F' in velocity fields indicate folded velocities (no de-aliasing),
    # velocities fields without an F have been unfolded.
    'V': corrected_velocity,
                        # Unfolded velocity from corrected timeseries
    'VF': velocity,     # Folded velocity from corr. t.s.
    'UV': None,         # Unfolded velocity from uncorrected timeseries
    'UVF': None,        # Folded velcity from uncorr. t.s.
    'Vh': corrected_velocity,
    'VH': corrected_velocity,
                        # Velocity from corr. t.s., horizontal channel
    'Vv': None,
    'VV': None,         # Velocity from corr. t.s., vertical channel
    'UVh': None,        # Velocity from uncorr. t.s., horizontal channel
    'UVv': None,        # Velocity from uncorr. t.s., vertical channel
    'VFh': velocity,
                        # Folded velocity, corr. t.s., horizontal channel
    'VFv': None,        # Folded velocity, corr. t.s., vertical channel
    'UnV': None,        # Folded velocity from uncorrected timeseries
    'UnVFh': None,      # Folded velocity, uncorr. t.s., horizontal channel
    'UnVFv': None,      # Folded velocity, uncorr. t.s., vertical channel

    # 'C' in spectral width fields indicate that the field has been
    # corrected for decorrelation causes by antenna rotation
    'W': spectrum_width,
                        # Spectral width from clutter corrected time series.
    'UW': None,         # Spectral width from uncorrected timeseries.
    'CW': None,         # Spec. width, antenna rotation corrected, corr. t.s.
    'UCW': None,        # Spec. width, antenna rotation corrected, uncorr t.s.
    'Wh': spectrum_width,
    'WH': spectrum_width,
                        # Spectral width, corr t.s., horizontal channel
    'Wv': None,
    'WV': None,         # Spectral width, corr t.s., vertical channel
    'UWh': None,        # Spectral width, uncorr t.s., horizontal channel
    'UWv': None,        # Spectral width, uncorr t.s., vertical channel
    'CWh': None,        # Spec. width, antenna rot. corr., horizontal channel
    'CWv': None,        # Spec. width, antenna rot. corr., vertical channel

    'SQI': normalized_coherent_power,
                        # Signal quality index
    'SQIh': normalized_coherent_power,
                        # Signal quality index, horizontal channel
    'SQIv': None,       # Signal quality index, vertical channel

    'CCOR': None,       # Clutter power correction
    'CCORh': None,      # Clutter power correction, horizontal channel
    'CCORv': None,      # Clutter power correction, vertical channel

    'SIGPOW': None,     # Singal Power
    'SNR': None,        # Raw signal to noise ratio
    'SNRh': None,       # Raw signal to noise ratio, horizontal channel
    'SNRv': None,       # Raw signal to noise ration, vertical channel

    'DFT': None,        # Signal spectrum amplitude
    'DFTh': None,       # Signal spectrum amplitude, horizontal channel
    'DFTv': None,       # Signal spectrum amplitude, vertical channel

    'LOG': None,        # Logarithmic amplitude 10*log Isq + Qsq
    'LOGh': None,       # Logarithmic amplitude, horizontal channel
    'LOGv': None,       # Logarithmic amplitude, vertical channel

    'CMAP': None,       # Censor map

    # A '1' ending on differential reflectivity fields indicates that the
    # moment has calculated using a 1st LAG algorithm.
    'ZDR': corrected_differential_reflectivity,
                        # Differential reflectivity from corrected timeseries
    'UZDR': differential_reflectivity,
                        # Diff. refl. from uncorrected timeseries
    'AZDR': None,       # Diff. refl., rainfall atten. corr., corr t.s.
    'ZDR1': None,       # Diff. refl., corr. t.s., 1st LAG algo.
    'UZDR1': None,      # Diff. refl., uncorr. t.s., 1st LAG algo.
    'AZDR1': None,      # Diff. refl., rain. atten. corr., corr. t.s., 1st LAG

    'PHI': corrected_differential_phase,
                        # Differential phase from corrected timeseries
    'PHIDP': corrected_differential_phase,
                        # Differential phase from corrected timeseries
    'UPHIDP': differential_phase,
                        # Diff. phase from uncorrected timeseries.
    'PHIH': None,       # Diff. phase, corr. t.s., horizontal channel
    'UPHIH': None,      # Diff. phase, uncorr t.s., hortizontal channel

    'KDP': specific_differential_phase,
                        # Specific differential phase

    'RHO': cross_correlation_ratio,
                        # Cross correlation coefficient from corrected t.s.
    'RHOHV': cross_correlation_ratio,
                        # Cross correlation coefficient from corrected t.s.
    'URHOHV': None,     # Cross correlation coefficient from uncorr. t.s.
    'RHOH': None,       # Cross corr., corr. t.s., horizontal transmit only
    'URHOH': None,      # Cross corr., uncorr t.s., horizontal transmit only

    'LDR': linear_depolarization_ratio,
                        # Linear depolarization ratio from corr. t.s.
    'ULDR': None,       # Linear depolarization ratio  from uncorr. t.s.
}

# UF field mappings
uf_field_mapping = {
    # UF 2 letter field: radar field
    'DZ': reflectivity,
    'CZ': corrected_reflectivity,
    'ZE': corrected_reflectivity,
    'ZT': total_power,
    'UZ': total_power,
    'VR': velocity,
    'VE': corrected_velocity,
    'SW': spectrum_width,
    'ZD': differential_reflectivity,
    'DR': corrected_differential_reflectivity,
    'CD': corrected_differential_reflectivity,
    'LR': linear_depolarization_ratio,
    'PH': differential_phase,
    'KD': specific_differential_phase,
    'RH': cross_correlation_ratio,
    'SQ': normalized_coherent_power,
    'NP': normalized_coherent_power,
    'HC': radar_echo_classification,
}

# Mapping used when writing UF files
write_uf_mapping = {
    # UF 2 letter field: radar field
    reflectivity: 'DZ',
    corrected_reflectivity: 'CZ',
    total_power: 'ZT',
    velocity: 'VR',
    corrected_velocity: 'VE',
    spectrum_width: 'SW',
    differential_reflectivity: 'ZD',
    corrected_differential_reflectivity: 'DR',
    linear_depolarization_ratio: 'LR',
    differential_phase: 'PH',
    specific_differential_phase: 'KD',
    cross_correlation_ratio: 'RH',
    normalized_coherent_power: 'SQ',
    radar_echo_classification: 'HC',
}

FIELD_MAPPINGS = {                  # Required variable
    'sigmet': sigmet_field_mapping,
    'nexrad_archive': nexrad_archive_field_mapping,
    'nexrad_cdm': nexrad_cdm_field_mapping,
    'nexrad_level3': nexrad_level3_mapping,
    'cfradial': cfradial_field_mapping,
    'mdv': mdv_field_mapping,
    'rsl': rsl_field_mapping,
    'chl': chl_field_mapping,
    'gamic': gamic_field_mapping,
    'uf': uf_field_mapping,
    'write_uf': write_uf_mapping,
}


def velocity_limit(container=None, selection=0):
    import pyart
    if isinstance(container, pyart.core.Radar):
        try:
            if selection >= 0 and selection < container.nsweeps:
                vel = container.get_nyquist_vel(selection,
                                                check_uniform=False)
            else:
                vel = container.get_nyquist_vel(0, check_uniform=False)
            return (-vel, vel)
        except LookupError:
            return (-30., 30.)
    else:
        return (-30., 30.)


def spectrum_width_limit(container=None, selection=0):
    import pyart
    if isinstance(container, pyart.core.Radar):
        try:
            if selection >= 0 and selection < container.nsweeps:
                vel = container.get_nyquist_vel(selection,
                                                check_uniform=False)
            else:
                vel = container.get_nyquist_vel(0, check_uniform=False)
            return (0, vel)
        except LookupError:
            return (0, 30.)
    else:
        return (0, 30.)

# map each field to a colormap

DEFAULT_FIELD_COLORMAP = {
    # field name : colormap
    reflectivity: 'pyart_HomeyerRainbow',
    corrected_reflectivity: 'pyart_HomeyerRainbow',
    total_power: 'pyart_HomeyerRainbow',
    signal_to_noise_ratio: 'pyart_Carbone17',

    velocity: 'pyart_BuDRd18',
    corrected_velocity: 'pyart_BuDRd18',
    simulated_velocity: 'pyart_BuDRd18',
    eastward_wind_component: 'pyart_BuDRd18',
    northward_wind_component: 'pyart_BuDRd18',
    vertical_wind_component: 'pyart_BuDRd18',

    spectrum_width: 'pyart_NWS_SPW',

    normalized_coherent_power: 'pyart_Carbone17',

    differential_reflectivity: 'pyart_RefDiff',
    corrected_differential_reflectivity: 'pyart_RefDiff',

    clutter_filter_power_removed: 'pyart_RefDiff',

    cross_correlation_ratio: 'pyart_RefDiff',

    differential_phase: 'pyart_Wild25',
    unfolded_differential_phase: 'pyart_Wild25',
    corrected_differential_phase: 'pyart_Wild25',

    specific_differential_phase: 'pyart_Theodore16',
    corrected_specific_differential_phase: 'pyart_Theodore16',

    linear_depolarization_ratio: 'pyart_SCook18',
    linear_depolarization_ratio_h: 'pyart_SCook18',
    linear_depolarization_ratio_v: 'pyart_SCook18',

    rain_rate: 'pyart_RRate11',
    radar_estimated_rain_rate: 'pyart_RRate11',

    radar_echo_classification: 'pyart_LangRainbow12',

    specific_attenuation: 'pyart_Carbone17',

    differential_phase_texture: 'pyart_BlueBrown11',
    height: 'pyart_SCook18',
    interpolated_profile: 'pyart_SCook18',

    # Additional reflectivity like fields
    'CZ': 'pyart_HomeyerRainbow',
    'DZ': 'pyart_HomeyerRainbow',
    'AZ': 'pyart_HomeyerRainbow',
    'Z': 'pyart_HomeyerRainbow',
    'dbz': 'pyart_HomeyerRainbow',
    'DBZ': 'pyart_HomeyerRainbow',
    'dBZ': 'pyart_HomeyerRainbow',
    'DBZH': 'pyart_HomeyerRainbow',
    'DBZ_S': 'pyart_HomeyerRainbow',
    'DBZ_K': 'pyart_HomeyerRainbow',
    'reflectivity_horizontal': 'pyart_HomeyerRainbow',
    'corr_reflectivity': 'pyart_HomeyerRainbow',
}

# map each field to a limit or a limit function

DEFAULT_FIELD_LIMITS = {
    # field name : limits
    reflectivity: (-30., 75.),
    corrected_reflectivity: (-30., 75.),
    total_power: (-30., 75.),
    signal_to_noise_ratio: (-20, 30.),

    velocity: velocity_limit,
    corrected_velocity: velocity_limit,
    eastward_wind_component: velocity_limit,
    northward_wind_component: velocity_limit,
    vertical_wind_component: velocity_limit,

    spectrum_width: spectrum_width_limit,

    normalized_coherent_power: (0., 1.),

    differential_reflectivity: (-1., 8.),
    corrected_differential_reflectivity: (-1., 8.),

    cross_correlation_ratio: (0.5, 1.05),

    differential_phase: (-180, 180.),
    unfolded_differential_phase: (-360, 360.),
    corrected_differential_phase: (-360, 360.),

    specific_differential_phase: (-2., 5.),
    corrected_specific_differential_phase: (-2., 5.),

    linear_depolarization_ratio: (-40., 0.),
    linear_depolarization_ratio_h: (-40., 0.),
    linear_depolarization_ratio_v: (-40., 0.),

    rain_rate: (0., 50.),
    radar_estimated_rain_rate: (0., 50.),

    radar_echo_classification: (0, 11),

    specific_attenuation: (0., 10.),

    differential_phase_texture: (0, 180.),

    height: (0, 20000),
    interpolated_profile: (0, 10000),

    # Additional reflectivity like fields
    'CZ': (-10., 65.),
    'DZ': (-10., 65.),
    'AZ': (-10., 65.),
    'Z': (-10., 65.),
    'dbz': (-10., 65.),
    'DBZ': (-10., 65.),
    'dBZ': (-10., 65.),
    'DBZH': (-10., 65.),
    'DBZ_S': (-10., 65.),
    'DBZ_K': (-10., 65.),
    'reflectivity_horizontal': (-10., 65.),
    'corr_reflectivity': (-10., 65.),
}
