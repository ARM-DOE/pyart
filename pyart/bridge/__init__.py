"""
================================================
bridging to other toolkits (:mod:`pyart.bridge`)
================================================

.. currentmodule:: pyart.bridge

Py-ART  as a bridge to other community software projetcs.

.. autosummary::
    :toctree: generated/

"""

try:
    from .wradlib_bridge import texture_of_complex_phase
    _WRADLIB_AVAILABLE = True
except ImportError:
    _WRADLIB_AVAILABLE = False

if _WRADLIB_AVAILABLE:
    from .. import retrieve
    retrieve.texture_of_complex_phase = texture_of_complex_phase


