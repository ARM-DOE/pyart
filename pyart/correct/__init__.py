"""
========================================
Radar Corrections (:mod:`pyart.correct`)
========================================

.. currentmodule:: pyart.correct

Py-ART has the ability to perform a number of common corrections on radar
moments and data.

.. autosummary::
    :toctree: generated/

    dealias_fourdd
    dealias_unwrap_phase
    dealias_region_based
    calculate_attenuation
    phase_proc_lp
    GateFilter
    moment_based_gate_filter
    find_time_in_interp_sonde

"""

try:
    from .dealias import dealias_fourdd, find_time_in_interp_sonde
except ImportError:
    pass
from .attenuation import calculate_attenuation
from .phase_proc import phase_proc_lp
from .filters import GateFilter, moment_based_gate_filter
from .unwrap import dealias_unwrap_phase
from .region_dealias import dealias_region_based

__all__ = [s for s in dir() if not s.startswith('_')]
