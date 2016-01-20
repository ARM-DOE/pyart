"""
========================
Core (:mod:`pyart.core`)
========================

.. currentmodule:: pyart.core

Core Py-ART classes and function for interacting with weather radar data.

Core classes
============

.. autosummary::
    :toctree: generated/

    Radar
    Grid

"""

from .radar import Radar
from .grid import Grid

# Depreciated function names in this name space
from ..exceptions import _depreciated_alias
from ..util import radar_utils as _radar_utils
is_vpt =  _depreciated_alias(
    _radar_utils.is_vpt,
    'pyart.core.radar.is_vpt',
    'pyart.util.radar_utils.is_vpt')
to_vpt =  _depreciated_alias(
    _radar_utils.to_vpt,
    'pyart.core.radar.to_vpt',
    'pyart.util.radar_utils.to_vpt')

__all__ = [s for s in dir() if not s.startswith('_')]
