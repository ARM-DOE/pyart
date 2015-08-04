"""
=============================
Utilities (:mod:`pyart.util`)
=============================

Miscellaneous utility functions.

The location and names of these functions within Py-ART may change between
versions without depeciation, use with caution.

.. currentmodule:: pyart.util

.. autosummary::
    :toctree: generated/

    angular_mean
    angular_std
    angular_mean_deg
    angular_std_deg
    interval_mean
    interval_std

"""

from .circular_stats import angular_mean, angular_std
from .circular_stats import angular_mean_deg, angular_std_deg
from .circular_stats import interval_mean, interval_std

__all__ = [s for s in dir() if not s.startswith('_')]
