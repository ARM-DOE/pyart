

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('pyart', parent_package, top_path)
    config.add_subpackage('io')     # io first to detect if RSL is missing.
    config.add_subpackage('__check_build')
    config.add_subpackage('core')
    config.add_subpackage('correct')
    config.add_subpackage('graph')
    config.add_subpackage('map')
    config.add_subpackage('retrieve')
    config.add_subpackage('filters')
    config.add_subpackage('testing')
    config.add_subpackage('util')
    config.add_subpackage('aux_io')
    config.add_subpackage('bridge')

    config.add_data_dir('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
