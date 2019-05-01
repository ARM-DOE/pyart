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
    cross_section_rhi
    datetime_from_radar
    datetimes_from_radar
    datetime_from_dataset
    datetimes_from_dataset
    datetime_from_grid
    estimate_noise_hs74
    is_vpt
    to_vpt
    join_radar
    simulated_vel_from_profile
    texture_along_ray
    texture
    rolling_window
    angular_texture_2d

"""

from .circular_stats import angular_mean, angular_std
from .circular_stats import angular_mean_deg, angular_std_deg
from .circular_stats import interval_mean, interval_std
from .circular_stats import mean_of_two_angles, mean_of_two_angles_deg
from .datetime_utils import datetime_from_radar, datetimes_from_radar
from .datetime_utils import datetime_from_dataset, datetimes_from_dataset
from .datetime_utils import datetime_from_grid
from .xsect import cross_section_ppi, cross_section_rhi
from .hildebrand_sekhon import estimate_noise_hs74
from .radar_utils import is_vpt, to_vpt, join_radar
from .simulated_vel import simulated_vel_from_profile
from .sigmath import texture_along_ray, rolling_window
from .sigmath import texture, angular_texture_2d

__all__ = [s for s in dir() if not s.startswith('_')]
