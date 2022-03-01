"""
Correct radar fields.

"""

from .dealias import dealias_fourdd
from .attenuation import calculate_attenuation
from .phase_proc import phase_proc_lp, phase_proc_lp_gf
from .attenuation import calculate_attenuation_zphi
from .attenuation import calculate_attenuation_philinear
# for backwards compatibility GateFilter available in the correct namespace
from ..filters.gatefilter import GateFilter, moment_based_gate_filter
from .unwrap import dealias_unwrap_phase
from .region_dealias import dealias_region_based
from .despeckle import find_objects, despeckle_field
from .bias_and_noise import correct_noise_rhohv, correct_bias

__all__ = [s for s in dir() if not s.startswith('_')]
