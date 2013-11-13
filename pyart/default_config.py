# dictionary of standard metadata for various parameters
DEFAULT_METADATA = {
    # metadata for radar fields, assuming a stationary platform
    'reflectivity_horizontal': {
        'units': 'dBZ',
        'standard_name': 'equivalent_reflectivity_factor',
        'long_name': 'equivalent_reflectivity_factor',
        'valid_max': 80.0,
        'valid_min': -45.0,
        'coordinates': 'elevation azimuth range'},

    'reflectivity_horizontal_filtered': {
        'units': 'dBZ',
        'standard_name': 'equivalent_reflectivity_factor_filtered',
        'long_name': 'equivalent_reflectivity_factor_filtered',
        'valid_max': 80.0,
        'valid_min': -45.0,
        'coordinates': 'elevation azimuth range'},

    'mean_doppler_velocity': {
        'units': 'meters_per_second',
        'standard_name': (
            'radial_velocity_of_scatterers_away_from_instrument'),
        'long_name': (
            'radial_velocity_of_scatterers_away_from_instrument'),
        'valid_max': 95.0,
        'valid_min': -95.0,
        'coordinates': 'elevation azimuth range'},

    'diff_phase': {
        'units': 'degrees/km',
        'standard_name': 'specific_differential_phase_hv',
        'long_name': 'specific_differential_phase_hv',
        'valid_max': 20.0,
        'valid_min': -10.0,
        'coordinates': 'elevation azimuth range'},

    'diff_reflectivity': {
        'units': 'dB',
        'standard_name': 'log_differential_reflectivity_hv',
        'long_name': 'log_differential_reflectivity_hv',
        'valid_max': 8.0,
        'valid_min': -6.0,
        'coordinates': 'elevation azimuth range'},

    'copol_coeff': {
        'units': 'ratio',
        'standard_name': 'cross_correlation_ratio_hv',
        'long_name': 'cross_correlation_ratio_hv',
        'valid_max': 1.0,
        'valid_min': 0.0,
        'coordinates': 'elevation azimuth range'},

    'norm_coherent_power': {
        'units': 'ratio',
        'standard_name': 'signal_quality',
        'long_name': 'signal_quality',
        'valid_max': 1.0,
        'valid_min': 0.0,
        'comment': 'Also know as Normalized Coherent Power',
        'coordinates': 'elevation azimuth range'},

    'doppler_spectral_width': {
        'units': 'meters_per_second',
        'standard_name': 'spectrum_width',
        'long_name': 'spectrum_width',
        'valid_max': 45.0,
        'valid_min': 0.0,
        'coordinates': 'elevation azimuth range'},

    'dp_phase_shift': {
        'units': 'degrees',
        'standard_name': 'differential_phase_hv',
        'long_name': 'differential_phase_hv',
        'valid_max': 180.0,
        'valid_min': -180.0,
        'coordinates': 'elevation azimuth range'},

    'corrected_mean_doppler_velocity': {
        'units': 'meters_per_second',
        'standard_name': (
            'radial_velocity_of_scatterers_away_from_instrument'),
        'long_name': (
            'radial_velocity_of_scatterers_away_from_instrument'),
        'valid_max': 45.0,
        'valid_min': -45.0,
        'coordinates': 'elevation azimuth range'},

    'unfolded_dp_phase_shift': {
        'units': 'degrees',
        'standard_name': 'differential_phase_hv',
        'long_name': 'differential_phase_hv',
        'valid_max': 480.0,
        'valid_min': 0.0,
        'coordinates': 'elevation azimuth range'},

    'attenuation_corrected_reflectivity_horizontal': {
        'units': 'dBZ',
        'standard_name': 'equivalent_reflectivity_factor',
        'long_name': 'equivalent_reflectivity_factor',
        'valid_max': 80.0,
        'valid_min': -45.0,
        'coordinates': 'elevation azimuth range'},

    'recalculated_diff_phase': {
        'units': 'degrees/km',
        'standard_name': 'specific_differential_phase_hv',
        'long_name': 'specific_differential_phase_hv',
        'valid_max': 20.0,
        'valid_min': -1.0,
        'least_significant_digit': 3},

    # metadata for radar attributes
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

    # metadata for radar location dictionaries
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

    # metadata for instrument parameter dictionary
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

}


# File specific metadata
SIGMET_METADATA = {}

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

FILE_SPECIFIC_METADATA = {
    'sigmet': SIGMET_METADATA,
    'nexrad_archive': NEXRAD_METADATA,
    'nexrad_cdm': NEXRAD_METADATA,
}


# FIELD_MAPPINGS

# This dictionary maps sigmet data types -> radar field names
# Users can pass their own version of this dictionary to the read_sigmet
# function as the field_names parameter.
SIGMET_FIELD_MAPPING = {
    'XHDR': None,                       # (0) Extended Header, cannot be mapped
    'DBT': None,                        # (1) Total Power
    'DBZ': None,                        # (2) Reflectivity
    'VEL': None,                        # (3) Velocity
    'WIDTH': None,                      # (4) Width
    'ZDR': None,                        # (5) Differential reflectivity
    'DBZC': None,                       # (7) Corrected reflectivity
    'DBT2': None,                       # (8) Total Power
    'DBZ2': 'reflectivity_horizontal',  # (9) Reflectivity
    'VEL2': 'mean_doppler_velocity',    # (10) Velocity
    'WIDTH2': None,                     # (11) Width
    'ZDR2': 'diff_reflectivity',        # (12) Differential reflectivity
    'RAINRATE2': None,                  # (13) Rainfall rate
    'KDP': None,                        # (14) KDP (differential phase)
    'KDP2': 'diff_phase',               # (15) KDP (differential phase)
    'PHIDP': None,                      # (16) PhiDP (differential phase)
    'VELC': None,                       # (17) Corrected velocity
    'SQI': None,                        # (18) SQI
    'RHOHV': None,                      # (19) RhoHV
    'RHOHV2': 'copol_coeff',            # (20) RhoHV
    'DBZC2': 'reflectivity_horizontal_filtered',    # (21) Corrected Reflec.
    'VELC2': 'VELC2',                   # (21) Corrected Velocity
    'SQI2': 'norm_coherent_power',      # (23) SQI
    'PHIDP2': 'dp_phase_shift',         # (24) PhiDP (differential phase)
    'LDRH': None,                       # (25) LDR xmt H, rcv V
    'LDRH2': None,                      # (26) LDR xmt H, rcv V
    'LDRV': None,                       # (27) LDR xmt V, rcv H
    'LDRV2': None,                      # (28) LDR xmt V, rcv H
    'HEIGHT': None,                     # (32) Height (1/10 km)
    'VIL2': None,                       # (33) Linear Liquid
    'RAW': None,                        # (34) Raw Data
    'SHEAR': None,                      # (35) Wind Shear
    'DIVERGE2': None,                   # (36) Divergence
    'FLIQUID2': None,                   # (37) Floated liquid
    'USER': None,                       # (38) User type
    'OTHER': None,                      # (39) Unspecified
    'DEFORM2': None,                    # (40) Deformation
    'VVEL2': None,                      # (41) Vertical velocity
    'HVEL2': None,                      # (42) Horizontal velocity
    'HDIR2': None,                      # (43) Horizontal wind direction
    'AXDIL2': None,                     # (44) Axis of dilation
    'TIME2': None,                      # (45) Time in seconds
    'RHOH': None,                       # (46) Rho, xmt H, rcv V
    'RHOH2': None,                      # (47) Rho, xmt H, rcv V
    'RHOV': None,                       # (48) Rho, xmt V, rcv H
    'RHOV2': None,                      # (49) Rho, xmt V, rcv H
    'PHIH': None,                       # (50) Phi, xmt H, rcv V
    'PHIH2': None,                      # (51) Phi, xmt H, rcv V
    'PHIV': None,                       # (52) Phi, xmt V, rcv H
    'PHIV2': None,                      # (53) Phi, xmt V, rcv H
    'USER2': None,                      # (54) User type
    'HCLASS': None,                     # (55) Hydrometeor class
    'HCLASS2': None,                    # (56) Hydrometeor class
    'ZDRC': None,                       # (57) Corrected diff. refl.
    'ZDRC2': None                       # (58) Corrected diff. refl.
}


NEXRAD_ARCHIVE_FIELD_MAPPING = {
    'REF': 'reflectivity',
    'VEL': 'velocity',
    'SW': 'spectrum_width',
    'ZDR': 'differential_reflectivity',
    'PHI': 'differential_phase',
    'RHO': 'correlation_coefficient'
}

NEXRAD_CDM_FIELD_MAPPING = {
    'Reflectivity': 'reflectivity',
    'RadialVelocity': 'velocity',
    'SpectrumWidth': 'spectrum_width',
    'DifferentialReflectivity': 'differential_reflectivity',
    'DifferentialPhase': 'differential_phase',
    'CorrelationCoefficient': 'correlation_coefficient'
}

MDV_FIELD_MAPPING = {
    'DBZ_F': 'reflectivity_horizontal',
    'VEL_F': 'mean_doppler_velocity',
    'WIDTH_F': 'doppler_spectral_width',
    'ZDR_F': 'diff_reflectivity',
    'RHOHV_F': 'copol_coeff',
    'NCP_F': 'norm_coherent_power',
    'KDP_F': 'diff_phase',
    'PHIDP_F': 'dp_phase_shift',
    'VEL_COR': 'corrected_mean_doppler_velocity',
    'PHIDP_UNF': 'unfolded_dp_phase_shift',
    'KDP_SOB': 'recalculated_diff_phase',
    'DBZ_AC': 'attenuation_corrected_reflectivity_horizontal', }

CFRADIAL_FIELD_MAPPING = {}
RSL_FIELD_MAPPING = {
    'DZ': 'reflectivity_horizontal_filtered',   # reflectivity
    'VR': 'mean_doppler_velocity',              # velocity
    'SW': None,                                 # spectrum width
    'CZ': None,                                 # corrected reflectivity
    'ZT': 'reflectivity_horizontal',            # uncorrected reflectivity
    'DR': None,                                 # differential reflectivity
    'LR': None,                                 # another diff. reflectivity
    'ZD': 'differential_reflectivity',          # another diff. reflectivity
    'DM': None,                                 # received power
    'RH': 'copolor_coefficient',                # RhoHV
    'PH': 'differential_phase_shift',           # PhiDP
    'XZ': None,                                 # X-band reflectivity
    'CD': None,                                 # Corrected DR.
    'MZ': None,                                 # DZ mask
    'MD': None,                                 # DR Mask
    'ZE': None,                                 # edited reflectivity
    'VE': 'corrected_mean_doppler_velocity',    # edited velocity
    'KD': 'specific_differential_phase',        # specific diff. phase
    'TI': None,                                 # TIME (unknown)
    'DX': None,                                 # ???
    'CH': None,                                 # ???
    'AH': None,                                 # ???
    'CV': None,                                 # ???
    'AV': None,                                 # ???
    'SQ': 'normalized_coherent_power',          # Signal Quality Index (sigmet)
    'VS': None,                                 # Radial Vel. combined
    'VL': None,                                 # Radial Vel. combined
    'VG': None,                                 # Radial Vel. combined
    'VT': None,                                 # Radial Vel. combined
    'NP': None,                                 # Normalized Coherent Power
    'HC': None,                                 # Hydroclass
    'VC': None,                                 # Radial Vel. Corrected.
    'V2': None,                                 # Radial Vel cut 2
    'S2': None,                                 # Spectrum width cut 2
    'V3': None,                                 # Radial Vel cut 3
    'S3': None,                                 # Spectrum width cut 3
}

#
FIELD_MAPPINGS = {
    'sigmet': SIGMET_FIELD_MAPPING,
    'nexrad_archive': NEXRAD_ARCHIVE_FIELD_MAPPING,
    'nexrad_cdm': NEXRAD_CDM_FIELD_MAPPING,
    'cfradial': CFRADIAL_FIELD_MAPPING,
    'mdv': MDV_FIELD_MAPPING,
    'rsl': RSL_FIELD_MAPPING,
}
