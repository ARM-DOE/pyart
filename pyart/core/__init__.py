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

Coordinate transformations
==========================

.. autosummary::
    :toctree: generated/

    antenna_to_cartesian
    antenna_vectors_to_cartesian
    cartesian_to_geographic
    cartesian_vectors_to_geographic
    cartesian_to_geographic_aeqd
    geographic_to_cartesian
    geographic_to_cartesian_aeqd

"""

from .radar import Radar
from .grid import Grid

from .transforms import antenna_to_cartesian
from .transforms import antenna_vectors_to_cartesian
from .transforms import cartesian_to_geographic
from .transforms import cartesian_vectors_to_geographic
from .transforms import cartesian_to_geographic_aeqd
from .transforms import geographic_to_cartesian
from .transforms import geographic_to_cartesian_aeqd

# Deprecated function names in this name space
from ..exceptions import _deprecated_alias
from ..util import radar_utils as _radar_utils
is_vpt = _deprecated_alias(
    _radar_utils.is_vpt,
    'pyart.core.radar.is_vpt',
    'pyart.util.radar_utils.is_vpt')
to_vpt = _deprecated_alias(
    _radar_utils.to_vpt,
    'pyart.core.radar.to_vpt',
    'pyart.util.radar_utils.to_vpt')

__all__ = [s for s in dir() if not s.startswith('_')]
