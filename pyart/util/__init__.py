"""
Miscellaneous utility functions.

The location and names of these functions within Py-ART may change between
versions without depreciation, use with caution.

"""

from .circular_stats import (
    angular_mean,
    angular_mean_deg,
    angular_std,
    angular_std_deg,
    interval_mean,
    interval_std,
    mean_of_two_angles,
    mean_of_two_angles_deg,
)
from .columnsect import (
    for_azimuth,
    get_column_rays,
    get_field_location,
    sphere_distance,
)
from .datetime_utils import (
    datetime_from_dataset,
    datetime_from_grid,
    datetime_from_radar,
    datetimes_from_dataset,
    datetimes_from_radar,
)
from .hildebrand_sekhon import estimate_noise_hs74
from .radar_utils import image_mute_radar, is_vpt, join_radar, subset_radar, to_vpt
from .sigmath import angular_texture_2d, rolling_window, texture, texture_along_ray
from .simulated_vel import simulated_vel_from_profile
from .xsect import cross_section_ppi, cross_section_rhi

__all__ = [s for s in dir() if not s.startswith("_")]
