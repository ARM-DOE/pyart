"""
==============================
Filters (:mod:`pyart.filters`)
==============================

.. currentmodule:: pyart.filters

Classes for specifying what gates are included and excluded from routines.

Filtering radar data
====================

.. autosummary::
    :toctree: generated/

    GateFilter
    moment_based_gate_filter

"""

from .gatefilter import GateFilter, moment_based_gate_filter
from .gatefilter import moment_based_gate_filter2

__all__ = [s for s in dir() if not s.startswith('_')]
