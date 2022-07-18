""" Unit Tests for Py-ART's _debug_info module. """

import sys
import warnings

import pyart
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


def test_debug_info():
    # test to see that something is written when _debug_info is called
    # we don't care what is written, just that something is.
    buf = StringIO()
    pyart._debug_info(buf)
    assert len(buf.getvalue()) > 0


def test_debug_stdout():
    # lack of error is assumed to mean that call succeed
    pyart._debug_info()


# Class for Mocking ImportErrors from
# http://stackoverflow.com/questions/2481511/mocking-importerror-in-python
class DisableModules(object):

    def __init__(self, modules):
        self.modules = modules

    def find_module(self, fullname):
        if fullname in self.modules:
            raise ImportError('Debug import failure for %s' % fullname)


def test_debug_info_all_disabled():
    modules = ['numpy', 'scipy', 'matplotlib', 'netCDF4', 'cylp', 'glpk',
               'cvxopt', 'mpl_toolkits', 'platform']
    for module in modules:
        if module in sys.modules:
            del sys.modules[module]
    fail_loader = DisableModules(modules)
    sys.meta_path.append(fail_loader)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # test to see that something is written when _debug_info is called
        # we don't care what is written, just that something is.
        buf = StringIO()
        pyart._debug_info(buf)
        assert len(buf.getvalue()) > 0
    # remove the Mocked ImportErrors
    sys.meta_path.remove(fail_loader)
