

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('graph', parent_package, top_path)
    config.add_data_dir('tests')
    config.add_data_files('balance-rgb.txt')
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
