import os
import sys

from numpy import get_include


def guess_rsl_path():
    return {'darwin': '/usr/local/trmm',
            'linux2': '/usr/local/trmm',
            'linux': '/usr/local/trmm',
            'win32': 'XXX'}[sys.platform]


def check_rsl_path(rsl_lib_path, rsl_include_path):

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
    config = Configuration('correct', parent_package, top_path)
    config.add_data_dir('tests')

    # determine and verify the RSL library and include location
    rsl_path = os.environ.get('RSL_PATH')
    if rsl_path is None:
        rsl_path = guess_rsl_path()
    rsl_lib_path = os.path.join(rsl_path, 'lib')
    rsl_include_path = os.path.join(rsl_path, 'include')

    # build the FourDD dealiaser if RSL is installed
    if check_rsl_path(rsl_lib_path, rsl_include_path):
        fourdd_sources = ['src/dealias_fourdd.c', 'src/sounding_to_volume.c',
                          'src/helpers.c']
        # Cython wrapper around FourDD
        config.add_extension(
            '_fourdd_interface',
            sources=['_fourdd_interface.c'] + fourdd_sources,
            libraries=['rsl'],
            library_dirs=[rsl_lib_path],
            include_dirs=[rsl_include_path, 'src'] + [get_include()],
            runtime_library_dirs=[rsl_lib_path])

    # phase unwrap extensions
    config.add_extension('_unwrap_1d', sources=['_unwrap_1d.c'],
                         include_dirs=[get_include()])
    unwrap_sources_2d = ['_unwrap_2d.c', 'unwrap_2d_ljmu.c']
    config.add_extension('_unwrap_2d', sources=unwrap_sources_2d,
                         include_dirs=[get_include()])
    unwrap_sources_3d = ['_unwrap_3d.c', 'unwrap_3d_ljmu.c']
    config.add_extension('_unwrap_3d', sources=unwrap_sources_3d,
                         include_dirs=[get_include()])

    # _fast_edge_finder extension
    config.add_extension('_fast_edge_finder', sources=['_fast_edge_finder.c'],
                         include_dirs=[get_include()])
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
