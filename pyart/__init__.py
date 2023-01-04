"""
Py-ART: The Python ARM Radar Toolkit
=====================================

"""

# print information on citing Py-ART, this message can be suppressed by
# setting the PYART_QUIET environment variable
_citation_text = """
## You are using the Python ARM Radar Toolkit (Py-ART), an open source
## library for working with weather radar data. Py-ART is partly
## supported by the U.S. Department of Energy as part of the Atmospheric
## Radiation Measurement (ARM) Climate Research Facility, an Office of
## Science user facility.
##
## If you use this software to prepare a publication, please cite:
##
##     JJ Helmus and SM Collis, JORS 2016, doi: 10.5334/jors.119
"""
from os import environ as _environ

if "PYART_QUIET" not in _environ:
    print(_citation_text)

# Make sure that deprecation warnings get printed by default
import warnings as _warnings

_warnings.simplefilter("always", DeprecationWarning)

from pkg_resources import DistributionNotFound, get_distribution

# import subpackages
# print out helpful message if build fails or importing from source tree
from . import __check_build  # noqa
from . import aux_io  # noqa
from . import bridge  # noqa
from . import config  # noqa
from . import core  # noqa
from . import correct  # noqa
from . import filters  # noqa
from . import graph  # noqa
from . import io  # noqa
from . import map  # noqa
from . import retrieve  # noqa
from . import testing  # noqa
from . import util  # noqa
from ._debug_info import _debug_info  # noqa

# root level functions
from .config import load_config  # noqa

# Get the version
try:
    __version__ = get_distribution("arm_pyart").version
except DistributionNotFound:
    # package is not installed
    __version__ = "0.0.0"
