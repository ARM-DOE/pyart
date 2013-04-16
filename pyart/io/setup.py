
import os
import sys

from numpy import get_include


def guess_rsl_path():
    return {'darwin': '/usr/local/trmm',
            'linux2': '/usr/local/trmm',
            'win32': 'XXX'}[sys.platform]


def check_rsl_path(rsl_lib_path, rsl_include_path):

    ext = {'darwin': 'dylib', 'linux2': 'so', 'win32': 'DLL'}[sys.platform]
    lib_file = os.path.join(rsl_lib_path, 'librsl.' + ext)
    if os.path.isfile(lib_file) is False:
        raise ValueError('rsl library not found in ' + rsl_lib_path)

    inc_file = os.path.join(rsl_include_path, 'rsl.h')
    if os.path.isfile(inc_file) is False:
        raise ValueError('rsl.h not found in ' + rsl_include_path)
    return


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info, BlasNotFoundError
    config = Configuration('io', parent_package, top_path)

    config.add_data_dir('tests')

    # determine and verify the at RSL location
    rsl_path = os.environ.get('RSL_PATH')
    if rsl_path is None:
        rsl_path = guess_rsl_path()
    print "rsl_path set to", rsl_path
    rsl_lib_path = os.path.join(rsl_path, 'lib')
    rsl_include_path = os.path.join(rsl_path, 'include')
    check_rsl_path(rsl_lib_path, rsl_include_path)

    # Cython wrapper around RSL
    config.add_extension(
        '_rsl_interface',
        sources=['_rsl_interface.c'],
        libraries=['rsl'],
        library_dirs=[rsl_lib_path],
        include_dirs=[rsl_include_path] + [get_include()],
        runtime_library_dirs=[rsl_lib_path])

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
