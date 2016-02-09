"""
=============================
Utilities (:mod:`pyart.util`)
=============================

Miscellaneous utility functions.

The location and names of these functions within Py-ART may change between
versions without depeciation, use with caution.

.. currentmodule:: pyart.util

Direction statistics
====================

.. autosummary::
    :toctree: generated/

    angular_mean
    angular_std
    angular_mean_deg
    angular_std_deg
    interval_mean
    interval_std
    mean_of_two_angles
    mean_of_two_angles_deg

Miscellaneous functions
=======================

.. autosummary::
    :toctree: generated/

    cross_section_ppi
    estimate_noise_hs74
    is_vpt
    to_vpt
    join_radar

"""

from .circular_stats import angular_mean, angular_std
from .circular_stats import angular_mean_deg, angular_std_deg
from .circular_stats import interval_mean, interval_std
from .circular_stats import mean_of_two_angles, mean_of_two_angles_deg
from .xsect import cross_section_ppi
from .hildebrand_sekhon import estimate_noise_hs74
from .radar_utils import is_vpt, to_vpt, join_radar

__all__ = [s for s in dir() if not s.startswith('_')]
