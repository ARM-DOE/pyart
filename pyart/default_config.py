"""
Configuration file for the Python ARM Radar Toolkit (Py-ART)

The values for a number of Py-ART parameters and the default metadata created
when reading files, correcting fields, etc. is controlled by this single
Python configuration file.

Py-ART's configuration can be modified by setting the environmental
variable PYART_CONFIG to point to a configuration file with formatting similar
to this file.

The recommended method for changing these defaults is for users to copy this
file into their home directory, rename it to .pyart_config.py, make any
desired changes, and adjust their login scripts to set the PYART_CONFIG
environmental variable to point to .pyart_config.py in their home directory.

Py-ART's configuration can also be modified within a script or shell session
using the load_config functions, such modification will last until a the end
of the script/session or until a new configuration is loaded.

"""

##############################################################################
##############################################################################
# Simple configuration
#
# Adjust the values of the variable (right hand side of the equal sign) in
# this section for an easy method of customizing Py-ART.  Do not change the
# variables names (the left hand side of the equal sign).  More advanced
# settings are based upon these variables.  Most users will find that
# adjusting this section is all that is needed.
##############################################################################
##############################################################################

# The default fill value for masked arrays and _FillValue keys
fill_value = -9999.0

# Field names used when reading in radar files and in the various correction
# and retrieval algorithms. The comments in this section provide additional
# information about the fields in that section.

# radar reflectivity fields,
reflectivity = 'reflectivity'
corrected_reflectivity = 'corrected_reflectivity'
total_power = 'total_power'

# mean doppler velocity fields, VEL
velocity = 'velocity'
corrected_velocity = 'corrected_velocity'

# spectral width field, SW
spectrum_width = 'spectrum_width'

# differential reflectivity fields, ZDR
differential_reflectivity = 'differential_reflectivity'
corrected_differential_reflectivity = 'corrected_differential_reflectivity'

# cross correlation ratio, correlation coefficient, RhoHV
cross_correlation_ratio = 'cross_correlation_ratio'

# normalized coherent power, signal quality index, SQI, NCP
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
linear_depolarization_ratio_h = 'linear_polarization_ratio_h'
linear_depolarization_ratio_v = 'linear_polarization_ratio_v'

# Misc fields
signal_to_noise_ratio = 'signal_to_noise_ratio'
rain_rate = 'rain_rate'
radar_estimated_rain_rate = 'radar_estimated_rain_rate'
radar_echo_classification = 'radar_echo_classification'
specific_attenuation = 'specific_attenuation'

# End of Simple Configuration section

##############################################################################
##############################################################################
# Advanced Configuration
#
# Most users will not want to make any changes in this section.  For users
# who want a more fined grained control over Py-ART's configuration this
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

FILL_VALUE = fill_value     # the default fill value for masked arrays and
                            # the _FillValue key.

# The DEFAULT_FIELD_NAMES controls the field names which are used in
# by the correction and retrieval algorithms in Py-ART.  The keys of the
# dictionary are "internal" names which cannot change, the values are the field
# names which will be used in the algorithms by default.
# which will be used.  For best results use the names defined by the variables
# in simple configuration section which are also used in the DEFAULT_METADATA
# and FIELD_MAPPINGS variable.  If you choose to change a field name
# the names should also be changed in the DEFAULT_METADATA and FIELD_MAPPINGS
# variable. This is not required but highly suggested.

DEFAULT_FIELD_NAMES = {
    # internal field name (do not change): field name used (can change)
    'reflectivity': reflectivity,
    'corrected_reflectivity': corrected_reflectivity,
    'total_power': total_power,
    'velocity': velocity,
    'corrected_velocity': corrected_velocity,
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
}


##############################################################################
# Default metadata
#
# The DEFAULT_METADATA dictionary contains dictionaries which provide
# the default radar attribute and field metadata.  When reading in a file
# with Py-ART the FILE_SPECIFIC_METADATA variable is first queued for a
# metadata dictionary, if it is not found then the metadata in
# DEFAULT_METADATA is utilized.
##############################################################################

DEFAULT_METADATA = {

    # metadata for radar attributes
    # these closely follow the CF/Radial standard
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

    'sweep_mode': {
        'units': 'unitless',
        'standard_name': 'sweep_mode',
        'long_name': 'Sweep mode',
        'comment': ('Options are: "sector", "coplane", "rhi", '
                    '"vertical_pointing", "idle", "azimuth_surveillance", '
                    '"elevation_surveillance", "sunscan", "pointing", '
                    '"manual_ppi", "manual_rhi"')},

    'sweep_number': {
        'units': 'count',
        'standard_name': 'sweep_number',
        'long_name': 'Sweep number'},

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

    # metadata for radar sweep information dictionaries
    'sweep_start_ray_index': {
        'long_name': 'Index of first ray in sweep, 0-based',
        'units': 'count'},

    'sweep_end_ray_index': {
        'long_name': 'Index of last ray in sweep, 0-based',
        'units': 'count'},

    'fixed_angle': {
        'long_name': 'Target angle for sweep',
        'units': 'degrees',
        'standard_name': 'target_fixed_angle'},

    # metadata for radar location attributes
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

    # metadata for instrument_parameter dictionary
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

    # reflectivity fields
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

    # velocity fields
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

    # spectrum width fields
    spectrum_width: {
        'units': 'meters_per_second',
        'standard_name': 'doppler_spectrum_width',
        'long_name': 'Doppler spectrum width',
        'coordinates': 'elevation azimuth range'},

    # Dual Pol fields
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


    # depolarization ratio fields
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

    # misc. fields
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
}

# Metadata for CF/Radial files
cfradial_metadata = {}

# Metadata for MDV files
mdv_metadata = {}

# Metadata for RSL files
rsl_metadata = {}

FILE_SPECIFIC_METADATA = {      # Required
    'sigmet': sigmet_metadata,
    'nexrad_archive': nexrad_metadata,
    'nexrad_cdm': nexrad_metadata,
    'cfradial': cfradial_metadata,
    'mdv': mdv_metadata,
    'rsl': rsl_metadata
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
    'RHO': cross_correlation_ratio
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
    'DBZ_AC': corrected_reflectivity, }

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

FIELD_MAPPINGS = {                  # Required variable
    'sigmet': sigmet_field_mapping,
    'nexrad_archive': nexrad_archive_field_mapping,
    'nexrad_cdm': nexrad_cdm_field_mapping,
    'cfradial': cfradial_field_mapping,
    'mdv': mdv_field_mapping,
    'rsl': rsl_field_mapping,
}

# Locations of known WSR-88D radar stations. If location information is missing in archive, 
# we use these values if we know the code/name of the radar stations
# FIXME: The tuple is representing [Lat, Lon, Projected_X, Projected_Y, Altitude], the projection is WRF Lambert. However, projected coordination is redundant 
KNOWN_STATION_LOCATIONS = {
    "KABR": [45.456000, -98.413000, -1102460.950000, 1207729.970000, 1302.000000],
    "KABX": [35.150000, -106.824000, -1991098.570000, 247315.160000, 5870.000000],
    "KAKQ": [36.984000, -77.007000, 601088.300000, 233965.890000, 112.000000],
    "KAMA": [35.233000, -101.709000, -1548798.600000, 164283.490000, 3587.000000],
    "KAMX": [25.611000, -80.413000, 352954.600000, -1004082.010000, 14.000000],
    "KAPX": [44.907000, -84.720000, -55705.440000, 1067472.310000, 1464.000000],
    "KARX": [43.823000, -91.191000, -564446.880000, 970022.660000, 1276.000000],
    "KATX": [48.194000, -122.496000, -2766592.480000, 1973742.750000, 494.000000],
    "KBBX": [39.496000, -121.632000, -3059184.860000, 1072684.430000, 173.000000],
    "KBGM": [42.200000, -75.985000, 643244.650000, 799854.330000, 1606.000000],
    "KBHX": [40.498000, -124.292000, -3220909.080000, 1256838.540000, 2402.000000],
    "KBIS": [46.771000, -100.761000, -1256157.180000, 1377777.860000, 1658.000000],
    "KBIX": [30.524000, -88.985000, -463640.470000, -468425.440000, 136.000000],
    "KBLX": [45.854000, -108.607000, -1858689.230000, 1403430.480000, 3598.000000],
    "KBMX": [33.172000, -86.770000, -249730.540000, -192424.650000, 645.000000],
    "KBOX": [41.956000, -71.137000, 1033959.260000, 814833.490000, 118.000000],
    "KBRO": [25.916000, -97.419000, -1312147.090000, -888069.810000, 23.000000],
    "KBUF": [42.949000, -78.737000, 418329.530000, 865884.600000, 693.000000],
    "KBYX": [24.598000, -81.703000, 228612.550000, -1117952.210000, 8.000000],
    "KCAE": [33.949000, -81.118000, 257393.460000, -108878.620000, 231.000000],
    "KCBW": [46.039000, -67.806000, 1227137.650000, 1291641.680000, 746.000000],
    "KCBX": [43.491000, -116.236000, -2499326.440000, 1325243.690000, 3061.000000],
    "KCCX": [40.923000, -78.004000, 489808.140000, 650810.890000, 2405.000000],
    "KCLE": [41.413000, -81.860000, 173763.320000, 690642.680000, 763.000000],
    "KCLX": [32.656000, -81.042000, 268332.820000, -247239.970000, 97.000000],
    "KCRI": [35.238000, -97.460000, -1179790.200000, 106073.380000, 1295.000000],
    "KCRP": [27.784000, -97.511000, -1293503.530000, -686593.810000, 45.000000],
    "KCXX": [44.511000, -73.167000, 841179.840000, 1070390.780000, 317.000000],
    "KCYS": [41.152000, -104.806000, -1682939.010000, 838648.150000, 6128.000000],
    "KDAX": [38.501000, -121.678000, -3102565.090000, 974902.720000, 30.000000],
    "KDDC": [37.761000, -99.969000, -1354365.150000, 405807.390000, 2590.000000],
    "KDFX": [29.273000, -100.281000, -1530130.150000, -488228.490000, 1131.000000],
    "KDIX": [39.947000, -74.411000, 792760.010000, 569380.390000, 149.000000],
    "KDLH": [46.837000, -92.210000, -616953.770000, 1303773.960000, 1428.000000],
    "KDMX": [41.731000, -93.723000, -784897.210000, 761762.940000, 981.000000],
    "KDOX": [38.826000, -75.440000, 718352.270000, 441507.990000, 50.000000],
    "KDTX": [42.700000, -83.472000, 42131.920000, 827895.360000, 1072.000000],
    "KDVN": [41.612000, -90.581000, -532578.350000, 727990.540000, 754.000000],
    "KDYX": [32.538000, -99.254000, -1380408.330000, -157037.530000, 1517.000000],
    "KEAX": [38.810000, -94.264000, -861056.950000, 453430.770000, 995.000000],
    "KEMX": [31.894000, -110.630000, -2408666.590000, -6096.980000, 5202.000000],
    "KENX": [42.586000, -74.064000, 792750.690000, 855386.310000, 1826.000000],
    "KEOX": [31.461000, -85.459000, -134262.430000, -378493.860000, 434.000000],
    "KEPZ": [31.873000, -106.698000, -2060385.650000, -97308.530000, 4104.000000],
    "KESX": [35.701000, -114.891000, -2656655.430000, 494276.330000, 4867.000000],
    "KEVX": [30.564000, -85.921000, -178647.700000, -474098.430000, 140.000000],
    "KEWX": [29.704000, -98.028000, -1313359.530000, -474831.960000, 633.000000],
    "KEYX": [35.098000, -117.561000, -2899284.660000, 508211.840000, 2757.000000],
    "KFCX": [37.024000, -80.274000, 320301.070000, 222999.600000, 2868.000000],
    "KFDR": [34.362000, -98.976000, -1325890.270000, 32407.310000, 1267.000000],
    "KFDX": [34.635000, -103.630000, -1727356.900000, 133219.790000, 4650.000000],
    "KFFC": [33.364000, -84.566000, -50916.350000, -175198.570000, 858.000000],
    "KFSD": [43.588000, -96.729000, -1000537.070000, 988566.080000, 1430.000000],
    "KFSX": [34.574000, -111.198000, -2380565.510000, 284423.790000, 7417.000000],
    "KFTG": [39.787000, -104.546000, -1692590.030000, 690433.430000, 5497.000000],
    "KFWS": [32.573000, -97.303000, -1204498.060000, -178870.370000, 683.000000],
    "KGGW": [48.206000, -106.625000, -1653088.200000, 1618908.760000, 2276.000000],
    "KGJX": [39.062000, -108.214000, -2008076.140000, 683247.010000, 9992.000000],
    "KGLD": [39.367000, -101.700000, -1468914.640000, 600531.620000, 3651.000000],
    "KGRB": [44.498000, -88.111000, -319826.320000, 1029431.220000, 682.000000],
    "KGRK": [30.722000, -97.383000, -1238555.390000, -374762.490000, 538.000000],
    "KGRR": [42.894000, -85.545000, -122950.700000, 849737.840000, 778.000000],
    "KGSP": [34.883000, -82.220000, 157181.980000, -11121.080000, 940.000000],
    "KGWX": [33.897000, -88.329000, -386801.410000, -109726.490000, 476.000000],
    "KGYX": [43.891000, -70.256000, 1075259.660000, 1031815.380000, 409.000000],
    "KHDX": [33.076000, -106.123000, -1980398.210000, 16650.970000, 4222.000000],
    "KHGX": [29.472000, -95.079000, -1041374.780000, -534781.070000, 18.000000],
    "KHNX": [36.314000, -119.632000, -3024562.790000, 692771.750000, 243.000000],
    "KHPX": [36.737000, -87.285000, -283430.000000, 190895.830000, 576.000000],
    "KHTX": [34.931000, -86.083000, -183826.000000, -5454.680000, 1760.000000],
    "KICT": [37.655000, -97.443000, -1143126.020000, 362600.110000, 1335.000000],
    "KICX": [37.591000, -112.862000, -2428930.570000, 635444.340000, 10600.000000],
    "KILN": [39.420000, -83.822000, 14841.320000, 474181.590000, 1056.000000],
    "KILX": [40.151000, -89.337000, -440523.800000, 564705.140000, 582.000000],
    "KIND": [39.708000, -86.280000, -189370.310000, 507321.570000, 790.000000],
    "KINX": [36.175000, -95.565000, -1002774.310000, 184744.650000, 668.000000],
    "KIWA": [33.289000, -111.670000, -2458907.400000, 163548.350000, 1353.000000],
    "KIWX": [41.359000, -85.700000, -138140.590000, 684123.960000, 960.000000],
    "KDGX": [32.280000, -89.984000, -545026.590000, -275066.040000, 609.000000],
    "KJAX": [30.485000, -81.702000, 213899.510000, -481837.630000, 33.000000],
    "KJGX": [32.675000, -83.351000, 58868.550000, -249029.760000, 521.000000],
    "KJKL": [37.591000, -83.313000, 58647.390000, 277972.950000, 1364.000000],
    "KLBB": [33.654000, -101.814000, -1588258.570000, -542.050000, 3259.000000],
    "KLCH": [30.125000, -93.216000, -860271.990000, -482840.130000, 13.000000],
    "KLIX": [30.337000, -89.826000, -542960.940000, -484208.050000, 24.000000],
    "KLNX": [41.958000, -100.576000, -1329890.150000, 859513.690000, 2970.000000],
    "KLOT": [41.605000, -88.085000, -330767.290000, 716319.520000, 663.000000],
    "KLRX": [40.740000, -116.803000, -2638644.980000, 1058754.290000, 6744.000000],
    "KLSX": [38.699000, -90.683000, -562032.640000, 415732.630000, 608.000000],
    "KLTX": [33.989000, -78.429000, 497114.640000, -94316.650000, 64.000000],
    "KLVX": [37.975000, -85.944000, -165133.330000, 320597.990000, 719.000000],
    "KLWX": [38.975000, -77.478000, 546555.120000, 444409.380000, 272.000000],
    "KLZK": [34.836000, -92.262000, -729180.580000, 12960.740000, 568.000000],
    "KMAF": [31.943000, -102.189000, -1654940.510000, -174832.260000, 2868.000000],
    "KMAX": [42.081000, -122.717000, -3036410.590000, 1364146.540000, 7513.000000],
    "KMBX": [48.393000, -100.865000, -1233545.570000, 1554883.770000, 1493.000000],
    "KMHX": [34.776000, -76.876000, 629396.060000, -1281.750000, 31.000000],
    "KMKX": [42.968000, -88.551000, -361684.660000, 865125.580000, 958.000000],
    "KMLB": [28.113000, -80.654000, 320073.590000, -734323.130000, 35.000000],
    "KMOB": [30.679000, -88.240000, -393681.210000, -455021.510000, 208.000000],
    "KMPX": [44.849000, -93.566000, -739568.120000, 1096821.460000, 946.000000],
    "KMQT": [46.531000, -87.548000, -268075.340000, 1249391.840000, 1411.000000],
    "KMRX": [36.169000, -83.402000, 51973.600000, 125444.550000, 1337.000000],
    "KMSX": [47.041000, -113.986000, -2214712.080000, 1639499.640000, 7855.000000],
    "KMTX": [41.263000, -112.448000, -2282844.520000, 1003854.550000, 0.000000],
    "KMUX": [37.155000, -121.898000, -3173888.400000, 847919.870000, 3469.000000],
    "KMVX": [47.528000, -97.326000, -989305.510000, 1420847.440000, 986.000000],
    "KMXX": [32.537000, -85.790000, -162625.410000, -262549.660000, 400.000000],
    "KNKX": [32.919000, -117.042000, -2932884.600000, 272486.320000, 955.000000],
    "KNQA": [35.345000, -89.873000, -515390.110000, 52310.370000, 282.000000],
    "KOAX": [41.320000, -96.367000, -1002880.040000, 741638.040000, 1148.000000],
    "KOHX": [36.247000, -86.563000, -222515.670000, 136537.260000, 579.000000],
    "KOKX": [40.866000, -72.864000, 909011.770000, 681119.800000, 85.000000],
    "KOTX": [47.680000, -117.627000, -2450281.750000, 1792153.490000, 2384.000000],
    "KPAH": [37.068000, -88.772000, -409928.770000, 231583.110000, 392.000000],
    "KPBZ": [40.532000, -80.218000, 310676.370000, 599760.920000, 1185.000000],
    "KPDT": [45.691000, -118.853000, -2611309.530000, 1618624.930000, 1515.000000],
    "KPOE": [31.156000, -92.976000, -827869.980000, -374567.590000, 408.000000],
    "KPUX": [38.459000, -104.181000, -1691919.880000, 544496.810000, 5249.000000],
    "KRAX": [35.666000, -78.490000, 481649.480000, 84821.840000, 348.000000],
    "KRGX": [39.754000, -119.462000, -2880932.280000, 1033249.250000, 8299.000000],
    "KRIW": [43.066000, -108.477000, -1923478.170000, 1107123.350000, 5568.000000],
    "KRLX": [38.311000, -81.723000, 192583.500000, 357257.880000, 1080.000000],
    "KRMX": [43.468000, -75.458000, 673575.620000, 940114.970000, 1516.000000],
    "KRTX": [45.715000, -122.965000, -2903124.010000, 1735947.050000, 1572.000000],
    "KSFX": [43.106000, -112.686000, -2244338.280000, 1200177.130000, 4474.000000],
    "KSGF": [37.235000, -93.401000, -804963.920000, 277945.280000, 1278.000000],
    "KSHV": [32.451000, -93.841000, -893562.900000, -228666.630000, 273.000000],
    "KSJT": [31.371000, -100.493000, -1512354.290000, -262499.720000, 1890.000000],
    "KSOX": [33.818000, -117.636000, -2951359.740000, 381125.760000, 3027.000000],
    "KSRX": [35.291000, -94.362000, -908798.990000, 78932.750000, 640.000000],
    "KTBW": [27.706000, -82.402000, 153596.310000, -782442.670000, 41.000000],
    "KTFX": [47.460000, -111.385000, -2015203.260000, 1627850.940000, 3714.000000],
    "KTLH": [30.398000, -84.329000, -30657.440000, -493628.270000, 63.000000],
    "KTLX": [35.333000, -97.278000, -1162569.170000, 114001.550000, 1213.000000],
    "KTWX": [38.997000, -96.233000, -1022976.850000, 492182.690000, 1367.000000],
    "KTYX": [43.756000, -75.680000, 653482.400000, 969741.070000, 1846.000000],
    "KUDX": [44.125000, -102.830000, -1464167.690000, 1122324.120000, 3016.000000],
    "KUEX": [40.321000, -98.442000, -1185726.250000, 657991.210000, 1976.000000],
    "KVAX": [30.890000, -83.002000, 92461.350000, -440340.600000, 178.000000],
    "KVBX": [34.838000, -120.397000, -3143581.500000, 568896.770000, 1233.000000],
    "KVNX": [36.741000, -98.128000, -1214969.770000, 273709.390000, 1210.000000],
    "KVTX": [34.412000, -119.179000, -3058888.010000, 487940.960000, 2726.000000],
    "KVWX": [38.266600, -87.716600, -314474.640000, 356191.370000, 392.000000],
    "KYUX": [32.495000, -114.657000, -2742215.430000, 160789.640000, 174.000000],
    "LPLA": [38.730000, -27.322000, 4510528.670000, 1731855.550000, 3334.000000],
    "PABC": [60.792000, -161.876000, -4137893.050000, 4595614.390000, 162.000000],
    "PACG": [56.853000, -135.529000, -3140584.930000, 3244840.620000, 270.000000],
    "PAEC": [64.511000, -165.295000, -3940826.980000, 5053425.530000, 54.000000],
    "PAHG": [60.726000, -151.351000, -3682172.920000, 4172361.660000, 242.000000],
    "PAIH": [59.461000, -146.303000, -3535284.020000, 3867501.110000, 67.000000],
    "PAKC": [58.679000, -156.629000, -4083350.620000, 4194256.520000, 63.000000],
    "PAPD": [65.035000, -147.502000, -3182624.750000, 4459157.330000, 2593.000000],
    "PGUA": [13.454000, 144.808000, -10743955.540000, 6055945.750000, 264.000000],
    "PHKI": [21.894000, -159.552000, -7011835.270000, 1406494.130000, 179.000000],
    "PHKM": [20.126000, -155.778000, -6857333.310000, 998183.040000, 3812.000000],
    "PHMO": [21.133000, -157.180000, -6891003.660000, 1178361.220000, 1363.000000],
    "PHWA": [19.095000, -155.569000, -6916826.160000, 898110.510000, 1370.000000],
    "RKJK": [35.924000, 126.622000, -8590696.710000, 8210291.310000, 78.000000],
    "RKSG": [36.956000, 127.021000, -8478205.120000, 8182361.090000, 52.000000],
    "RODN": [26.302000, 127.910000, -9615542.820000, 8025415.240000, 218.000000],
    "TJUA": [18.116000, -66.078000, 1903018.010000, -1661338.570000, 2794.000000],
    "KJAN": [32.318000, -90.080000, -553508.400000, -270462.850000, 297.000000],
    "KOUN": [35.236058, -97.462350, -1180023.200000, 105895.290000, 1214.000000],
    "KLGX": [47.115800, -124.106900, -2922938.190000, 1910840.560000, 366.000000]
}
