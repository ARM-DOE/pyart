"""
Py-ART: The Python ARM Radar Toolkit
=====================================

"""
__all__ = ['sounding', 'io']

# versioning
from .version import git_revision as __git_revision__
from .version import version as __version__

# import subpackages
from . import io
from . import correct
from . import graph
from . import map
from . import metadata
from . import testing

# test function setup based on scikit-image test function
import imp as _imp
import os.path as _osp
import functools as _functools

pkg_dir = _osp.abspath(_osp.dirname(__file__))

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
        args = ['', pkg_dir, '--exe']
        if verbose:
            args.extend(['-v', '-s'])
        nose.run('pyart', argv=args)

# do not use `test` as function name as this leads to a recursion problem with
# the nose test suite
test = _test
test_verbose = _functools.partial(test, verbose=True)
test_verbose.__doc__ = test.__doc__
