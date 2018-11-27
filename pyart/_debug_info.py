"""
Print out Py-ART version information.

This file can also be run as a script to report on dependencies before a
build: python pyart/_debug_info.py

"""

from __future__ import print_function
import os
import sys


def _debug_info(stream=None):
    """
    Print out version and status information for debugging.

    This file can be run as a script from the source directory to report on
    dependecies before a build using: **python pyart/_debug_info.py**.

    Parameters
    ----------
    stream : file-like object
        Stream to print the information to, None prints to sys.stdout.

    """
    if stream is None:
        stream = sys.stdout

    # remove the current path from the import search path
    # if this is not done ./io is found and not the std library io module.
    current_dir = os.path.dirname(os.path.abspath(__file__))
    if current_dir in sys.path:
        sys.path.remove(current_dir)

    try:
        import pyart
        pyart_version = pyart.__version__
    except:
        pyart_version = "MISSING"

    try:
        import platform
        python_version = platform.python_version()
    except:
        python_version = "MISSING"

    try:
        import numpy
        numpy_version = numpy.__version__
    except:
        numpy_version = "MISSING"

    try:
        import numpy
        numpy_version = numpy.__version__
    except:
        numpy_version = "MISSING"

    try:
        import scipy
        scipy_version = scipy.__version__
    except:
        scipy_version = "MISSING"

    try:
        import matplotlib
        matplotlib_version = matplotlib.__version__
    except:
        matplotlib_version = "MISSING"

    try:
        import netCDF4
        netCDF4_version = netCDF4.__version__
    except:
        netCDF4_version = "MISSING"

    try:
        rsl_version = pyart.io._rsl_interface._RSL_VERSION_STR
    except:
        rsl_version = "MISSING"

    try:
        import cylp
        cylp_available = "Available"
    except:
        cylp_available = "MISSING"

    try:
        import glpk
        glpk_version = "%i.%i" % (glpk.env.version)
    except:
        glpk_version = "MISSING"

    try:
        import cvxopt.info
        cvxopt_version = cvxopt.info.version
    except:
        cvxopt_version = "MISSING"

    try:
        import cartopy
        cartopy_version = cartopy.__version__
    except:
        cartopy_version = "MISSING"

    try:
        import pytest
        pytest_version = pytest.__version__
    except:
        pytest_version = "MISSING"

    print("Py-ART version:", pyart_version, file=stream)
    print("", file=stream)

    print("---- Dependencies ----", file=stream)
    print("Python version:", python_version, file=stream)
    print("NumPy version:", numpy_version, file=stream)
    print("SciPy version:", scipy_version, file=stream)
    print("matplotlib version:", matplotlib_version, file=stream)
    print("netCDF4 version:", netCDF4_version, file=stream)
    print("", file=stream)

    print("---- Optional dependencies ----", file=stream)
    print("TRMM RSL version:", rsl_version, file=stream)
    print("CyLP:", cylp_available, file=stream)
    print("PyGLPK version:", glpk_version, file=stream)
    print("CVXOPT version:", cvxopt_version, file=stream)
    print("Cartopy version:", cartopy_version, file=stream)
    print("pytest version:", pytest_version, file=stream)

if __name__ == "__main__":

    _debug_info()
