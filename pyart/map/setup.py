import os
from os.path import join
import warnings
import numpy


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info, BlasNotFoundError
    config = Configuration('map', parent_package, top_path)
    config.add_data_dir('tests')
    libraries = []
    if os.name == 'posix':
        libraries.append('m')

    config.add_extension('ckdtree',
                         sources=['ckdtree.c'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension('_load_nn_field_data',
                         sources=['_load_nn_field_data.c'],
                         include_dirs=[numpy.get_include()])

    config.add_extension('_gate_to_grid_map',
                         sources=['_gate_to_grid_map.c'],
                         libraries=libraries)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
