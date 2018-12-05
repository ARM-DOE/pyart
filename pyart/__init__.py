"""
Py-ART: The Python ARM Radar Toolkit
=====================================

"""

from __future__ import print_function

# Detect if we're being called as part of Py-ART's setup procedure
try:
    __PYART_SETUP__
except NameError:
    __PYART_SETUP__ = False

if __PYART_SETUP__:
    import sys as _sys
    _sys.stderr.write("Running from Py-ART source directory.\n")
    del _sys
else:

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
    if 'PYART_QUIET' not in _environ:
        print(_citation_text)

    # Make sure that deprecation warnings get printed by default
    import warnings as _warnings
    _warnings.simplefilter("always", DeprecationWarning)

    # print out helpful message if build fails or importing from source tree
    from . import __check_build

    # versioning
    from .version import git_revision as __git_revision__
    from .version import version as __version__

    # import subpackages
    from . import core
    from . import io
    from . import correct
    from . import graph
    from . import map
    from . import filters
    from . import util
    from . import testing
    from . import config
    from . import aux_io
    from . import retrieve
    from . import bridge

    # root level functions
    from .config import load_config
    from ._debug_info import _debug_info

    # test function setup based on scikit-image test function
    import os.path as _osp
    import functools as _functools
    import sys as _sys

    try:
        if _sys.version_info[:2] >= (3, 4):
            import importlib as _importlib
            specs = _importlib.util.find_spec('pytest')
            specs.loader.load_module()
        else:
            import imp as _imp
            _imp.find_module('pytest')
    except (AttributeError, ImportError) as error:
        def _test(verbose=False):
            """
            This would invoke the Py-ART test suite, but pytest couldn't
            be imported so the test suite can not run.
            """
            raise ImportError("Could not load pytest. Unit tests not available."
                              " To run unit tests, please install pytest.")
    else:
        def _test(verbose=False):
            """
            Invoke the Py-ART test suite.
            """
            import pytest
            pkg_dir = _osp.abspath(_osp.dirname(__file__))
            args = [pkg_dir, '--pyargs', 'pyart']
            if verbose:
                args.extend(['-v', '-s'])
            pytest.main(args=args)

    # Do not use `test` as function name as this leads to a recursion problem
    # with the pytest test suite.
    test = _test
    test_verbose = _functools.partial(test, verbose=True)
    test_verbose.__doc__ = test.__doc__
