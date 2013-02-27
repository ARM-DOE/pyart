

def dms_to_d(dms):
    """ Degrees, minutes, seconds to degrees """
    return dms[0] + (dms[1] + dms[2] / 60.0) / 60.0


def csapr_standard_names():
    prop_names = {'DBZ_F': 'reflectivity_horizontal',
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
    return prop_names


def get_mdv_meta(radarobj, field):
    debug = True
    print "go"
    csapr_names = csapr_standard_names()
    moment_fixes = {
        'DBZ_F': {
            'units': 'dBZ',
            'standard_name': 'equivalent_reflectivity_factor',
            'long_name': 'equivalent_reflectivity_factor',
            'valid_max': 80.0,
            'valid_min': -45.0,
            'least_significant_digit': 2},

        'VEL_F': {
            'units': 'm/s',
            'standard_name': (
                'radial_velocity_of_scatterers_away_from_instrument'),
            'long_name': (
                'radial_velocity_of_scatterers_away_from_instrument'),
            'valid_max': 95.0,
            'valid_min': -95.0,
            'least_significant_digit': 2},

        'KDP_F': {
            'units': 'degrees/km',
            'standard_name': 'specific_differential_phase_hv',
            'long_name': 'specific_differential_phase_hv',
            'valid_max': 20.0,
            'valid_min': -10.0,
            'least_significant_digit': 2},

        'ZDR_F': {
            'units': 'dB',
            'standard_name': 'log_differential_reflectivity_hv',
            'long_name': 'log_differential_reflectivity_hv',
            'valid_max': 8.0,
            'valid_min': -6.0,
            'least_significant_digit': 3},

        'RHOHV_F': {
            'units': 'ratio',
            'standard_name': 'cross_correlation_ratio_hv',
            'long_name': 'cross_correlation_ratio_hv',
            'valid_max': 1.0,
            'valid_min': 0.0,
            'least_significant_digit': 5},

        'NCP_F': {
            'units': 'ratio',
            'standard_name': 'signal_quality',
            'long_name': 'signal_quality',
            'valid_max': 1.0,
            'valid_min': 0.0,
            'comment': 'Also know as Normalized Coherent Power',
            'least_significant_digit': 5},

        'WIDTH_F': {
            'units': 'm/s',
            'standard_name': 'spectrum_width',
            'long_name': 'spectrum_width',
            'valid_max': 45.0,
            'valid_min': 0.0,
            'least_significant_digit': 2},

        'PHIDP_F': {
            'units': 'degrees',
            'standard_name': 'differential_phase_hv',
            'long_name': 'differential_phase_hv',
            'valid_max': 180.0,
            'valid_min': -180.0,
            'least_significant_digit': 2},

        'VEL_COR': {
            'units': 'm/s',
            'standard_name': (
                'radial_velocity_of_scatterers_away_from_instrument'),
            'long_name': (
                'radial_velocity_of_scatterers_away_from_instrument'),
            'valid_max': 45.0,
            'valid_min': -45.0,
            'least_significant_digit': 2},

        'PHIDP_UNF': {
            'units': 'degrees',
            'standard_name': 'differential_phase_hv',
            'long_name': 'differential_phase_hv',
            'valid_max': 480.0,
            'valid_min': 0.0,
            'least_significant_digit': 2},

        'DBZ_AC': {
            'units': 'dBZ',
            'standard_name': 'equivalent_reflectivity_factor',
            'long_name': 'equivalent_reflectivity_factor',
            'valid_max': 80.0,
            'valid_min': -45.0,
            'least_significant_digit': 2},

        'KDP_SOB': {
            'units': 'degrees/km',
            'standard_name': 'specific_differential_phase_hv',
            'long_name': 'specific_differential_phase_hv',
            'valid_max': 20.0,
            'valid_min': -1.0,
            'least_significant_digit': 3}}
    return moment_fixes[field]

