"""
Correct radar fields.

"""

# for backwards compatibility GateFilter available in the correct namespace
from ..filters.gatefilter import GateFilter, moment_based_gate_filter  # noqa
from .attenuation import calculate_attenuation  # noqa
from .attenuation import calculate_attenuation_philinear  # noqa
from .attenuation import calculate_attenuation_zphi  # noqa
from .bias_and_noise import calc_zdr_offset  # noqa
from .bias_and_noise import calc_cloud_mask, calc_noise_floor, correct_bias  # noqa
from .bias_and_noise import (
    correct_noise_rhohv,  # noqa
    cloud_threshold,  # noqa
    range_correction,  # noqa
)  # noqa
from .dealias import dealias_fourdd  # noqa
from .despeckle import despeckle_field, find_objects  # noqa
from .phase_proc import phase_proc_lp, phase_proc_lp_gf  # noqa
from .region_dealias import dealias_region_based  # noqa
from .unwrap import dealias_unwrap_phase  # noqa

__all__ = [s for s in dir() if not s.startswith("_")]
