from os.path import join


def configuration(parent_package='', top_path=None):
    global config
    from numpy.distutils.misc_util import Configuration
    config = Configuration('retrieve', parent_package, top_path)
    config.add_data_dir('tests')

    # KDP processing Cython extension
    config.add_extension('_kdp_proc', sources=['_kdp_proc.c'])
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
