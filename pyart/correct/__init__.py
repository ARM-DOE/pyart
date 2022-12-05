"""
Correct radar fields.

"""

# for backwards compatibility GateFilter available in the correct namespace
from ..filters.gatefilter import GateFilter, moment_based_gate_filter
from .attenuation import (
    calculate_attenuation,
    calculate_attenuation_philinear,
    calculate_attenuation_zphi,
)
from .bias_and_noise import correct_bias, correct_noise_rhohv
from .dealias import dealias_fourdd
from .despeckle import despeckle_field, find_objects
from .phase_proc import phase_proc_lp, phase_proc_lp_gf
from .region_dealias import dealias_region_based
from .unwrap import dealias_unwrap_phase

__all__ = [s for s in dir() if not s.startswith("_")]
