"""
========================================
Radar Retrievals (:mod:`pyart.retrieve`)
========================================

.. currentmodule:: pyart.retrieve

Radar retrievals.

Radar retrievals
================

.. autosummary::
    :toctree: generated/

    kdp_maesaka
    calculate_snr_from_reflectivity
    calculate_velocity_texture
    compute_snr
    compute_l
    compute_cdr
    compute_noisedBZ
    fetch_radar_time_profile
    map_profile_to_gates
    steiner_conv_strat
    hydroclass_semisupervised
    get_freq_band
    texture_of_complex_phase
    grid_displacement_pc
    grid_shift
    est_rain_rate_zpoly
    est_rain_rate_z
    est_rain_rate_kdp
    est_rain_rate_a
    est_rain_rate_zkdp
    est_rain_rate_za
    est_rain_rate_hydro
    velocity_azimuth_display
    quasi_vertical_profile

"""

from .kdp_proc import kdp_maesaka, kdp_schneebeli, kdp_vulpiani
from .echo_class import steiner_conv_strat, hydroclass_semisupervised
from .echo_class import get_freq_band
from .gate_id import map_profile_to_gates, fetch_radar_time_profile
from .simple_moment_calculations import calculate_snr_from_reflectivity
from .simple_moment_calculations import calculate_velocity_texture
from .simple_moment_calculations import compute_snr, compute_l, compute_cdr
from .simple_moment_calculations import compute_noisedBZ
from .advection import grid_displacement_pc, grid_shift
from .qpe import est_rain_rate_zpoly, est_rain_rate_z, est_rain_rate_kdp
from .qpe import est_rain_rate_a, est_rain_rate_zkdp, est_rain_rate_za
from .qpe import est_rain_rate_hydro
from .vad import velocity_azimuth_display
from .qvp import quasi_vertical_profile

__all__ = [s for s in dir() if not s.startswith('_')]
