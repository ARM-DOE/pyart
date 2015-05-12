"""
Py-ART: The Python ARM Radar Toolkit
=====================================

"""

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
    from . import testing
    from . import config
    from . import aux_io
    from . import retrieve
    from . import bridge

    # root level functions
    from .config import load_config
    from ._debug_info import _debug_info

    # test function setup based on scikit-image test function
    import imp as _imp
    import os.path as _osp
    import functools as _functools

    try:
        _imp.find_module('nose')
    except ImportError:
        def _test(verbose=False):
            """
            This would invoke the Py-ART test suite, but nose couldn't be
            imported so the test suite can not run.
            """
            raise ImportError("Could not load nose. Unit tests not available.")
    else:
        def _test(verbose=False):
            """
            Invoke the Py-ART test suite.
            """
            import nose
            pkg_dir = _osp.abspath(_osp.dirname(__file__))
            args = ['', pkg_dir, '--exe']
            if verbose:
                args.extend(['-v', '-s'])
            nose.run('pyart', argv=args)

    # do not use `test` as function name as this leads to a recursion problem
    # with the nose test suite
    test = _test
    test_verbose = _functools.partial(test, verbose=True)
    test_verbose.__doc__ = test.__doc__
