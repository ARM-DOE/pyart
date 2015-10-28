"""
================================================
Bridging to other toolkits (:mod:`pyart.bridge`)
================================================

.. currentmodule:: pyart.bridge

Py-ART can act as bridge to other community software projects.

The functionality in this namespace is available in other pyart namespaces.

Phase functions
===============

.. autosummary::
    :toctree: generated/

    texture_of_complex_phase

"""

from .wradlib_bridge import texture_of_complex_phase
from .. import retrieve as _retrieve
_retrieve.texture_of_complex_phase = texture_of_complex_phase
