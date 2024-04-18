"""
Miscellaneous utility functions.

The location and names of these functions within Py-ART may change between
versions without depreciation, use with caution.

"""

from .circular_stats import angular_mean  # noqa
from .circular_stats import angular_mean_deg  # noqa
from .circular_stats import angular_std  # noqa
from .circular_stats import angular_std_deg  # noqa
from .circular_stats import interval_mean  # noqa
from .circular_stats import interval_std  # noqa
from .circular_stats import mean_of_two_angles  # noqa
from .circular_stats import mean_of_two_angles_deg  # noqa
from .columnsect import for_azimuth  # noqa
from .columnsect import get_column_rays  # noqa
from .columnsect import get_field_location  # noqa
from .columnsect import sphere_distance  # noqa
from .columnsect import column_vertical_profile  # noqa
from .datetime_utils import datetime_from_dataset  # noqa
from .datetime_utils import datetime_from_grid  # noqa
from .datetime_utils import datetime_from_radar  # noqa
from .datetime_utils import datetimes_from_dataset  # noqa
from .datetime_utils import datetimes_from_radar  # noqa
from .hildebrand_sekhon import estimate_noise_hs74  # noqa
from .radar_utils import (  # noqa
    image_mute_radar,
    is_vpt,
    join_radar,
    determine_sweeps,
    subset_radar,
    to_vpt,
)
from .sigmath import (  # noqa
    angular_texture_2d,
    rolling_window,
    texture,
    texture_along_ray,
)
from .simulated_vel import simulated_vel_from_profile  # noqa
from .xsect import cross_section_ppi, cross_section_rhi  # noqa

__all__ = [s for s in dir() if not s.startswith("_")]
