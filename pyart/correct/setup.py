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
        return False

    inc_file = os.path.join(rsl_include_path, 'rsl.h')
    if os.path.isfile(inc_file) is False:
        return False
    return True


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('correct', parent_package, top_path)
    config.add_data_dir('tests')

    # determine and verify the at RSL location
    rsl_path = os.environ.get('RSL_PATH')
    if rsl_path is None:
        rsl_path = guess_rsl_path()
    rsl_lib_path = os.path.join(rsl_path, 'lib')
    rsl_include_path = os.path.join(rsl_path, 'include')

    # build the FourDD dealiaser if RSL is installed
    if check_rsl_path(rsl_lib_path, rsl_include_path):
        fourdd_sources = ['src/findRay.c', 'src/firstGuess_noread.c',
                          'src/firstGuess.c', 'src/FourDD.c',
                          'src/prepVolume.c', 'src/unfoldVolume.c',
                          'src/window.c']
        # Cython wrapper around FourDD
        config.add_extension(
            '_fourdd_interface',
            sources=['_fourdd_interface.c'] + fourdd_sources,
            libraries=['rsl'],
            library_dirs=[rsl_lib_path],
            include_dirs=[rsl_include_path, 'src'] + [get_include()],
            runtime_library_dirs=[rsl_lib_path])

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
