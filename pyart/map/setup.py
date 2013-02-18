from os.path import join
import warnings
import numpy


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info, BlasNotFoundError

    config = Configuration('map', parent_package, top_path)
    config.add_extension('ball_tree',
                         sources=[join('src', 'ball_tree.cpp')],
                         depends=[join('src', 'BallTree.h'),
                                  join('src', 'BallTreePoint.h')],
                         include_dirs=[numpy.get_include()])
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
