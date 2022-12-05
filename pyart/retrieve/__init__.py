"""
Radar retrievals.

"""

from .advection import grid_displacement_pc, grid_shift
from .comp_z import composite_reflectivity
from .echo_class import (
    conv_strat_yuter,
    get_freq_band,
    hydroclass_semisupervised,
    steiner_conv_strat,
)
from .gate_id import fetch_radar_time_profile, map_profile_to_gates
from .kdp_proc import kdp_maesaka, kdp_schneebeli, kdp_vulpiani
from .qpe import (
    est_rain_rate_a,
    est_rain_rate_hydro,
    est_rain_rate_kdp,
    est_rain_rate_z,
    est_rain_rate_za,
    est_rain_rate_zkdp,
    est_rain_rate_zpoly,
)
from .qvp import quasi_vertical_profile
from .simple_moment_calculations import (
    calculate_snr_from_reflectivity,
    calculate_velocity_texture,
    compute_cdr,
    compute_l,
    compute_noisedBZ,
    compute_snr,
)
from .spectra_calculations import dealias_spectra, spectra_moments
from .vad import vad_browning, vad_michelson

__all__ = [s for s in dir() if not s.startswith("_")]
