
import os
import sys

from numpy import get_include

RSL_MISSING_WARNING = """
==============================================================================
WARNING: RSL LIBS AND HEADERS COULD NOT BE FOUND AT THE PROVIDED LOCATION.

Py-ART will be build without bindings to the NASA TRMM RSL library but some
functionality will not be available.  If this functionality is desired please
rebuild and reinstall Py-ART after verifying:

    1. The NASA TRMM RSL library is installed and accessable.  This package
       can be obtained from:
            http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/.

    2. The RSL_PATH environmental variable points to where RSL is installed
       (the libs in $RSL_PATH/lib, the headers in $RSL_PATH/include).
       Currently the RSL_PATH variable is set to: %s

==============================================================================
"""


def guess_rsl_path():
    return {'darwin': '/usr/local/trmm',
            'linux2': '/usr/local/trmm',
            'linux': '/usr/local/trmm',
            'win32': 'XXX'}[sys.platform]


def check_rsl_path(rsl_lib_path, rsl_include_path):
    """ check if the rsl path is valid, return False if invalid. """
    ext = {'darwin': 'dylib',
           'linux2': 'so',
           'linux': 'so',
           'win32': 'DLL'}[sys.platform]
    lib_file = os.path.join(rsl_lib_path, 'librsl.' + ext)
    if os.path.isfile(lib_file) is False:
        return False

    inc_file = os.path.join(rsl_include_path, 'rsl.h')
    if os.path.isfile(inc_file) is False:
        return False
    return True


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info, BlasNotFoundError
    config = Configuration('io', parent_package, top_path)
    config.add_data_dir('tests')

    # determine and verify the at RSL location
    rsl_path = os.environ.get('RSL_PATH')
    if rsl_path is None:
        rsl_path = guess_rsl_path()
    rsl_lib_path = os.path.join(rsl_path, 'lib')
    rsl_include_path = os.path.join(rsl_path, 'include')

    # build the RSL interface if RSL is installed
    if check_rsl_path(rsl_lib_path, rsl_include_path):
        config.add_extension(
            '_rsl_interface',
            sources=['_rsl_interface.c'],
            libraries=['rsl'],
            library_dirs=[rsl_lib_path],
            include_dirs=[rsl_include_path] + [get_include()],
            runtime_library_dirs=[rsl_lib_path])
    else:
        import warnings
        warnings.warn(RSL_MISSING_WARNING % (rsl_path))

    config.add_extension('_sigmetfile',
                         sources=['_sigmetfile.c'],
                         include_dirs=[get_include()])

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
