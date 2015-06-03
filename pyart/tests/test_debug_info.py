""" Unit Tests for Py-ART's config.py module. """

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
