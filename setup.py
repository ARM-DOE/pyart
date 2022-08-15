#!/usr/bin/env python
"""Py-ART: Python ARM Radar Toolkit
The Python ARM Radar Toolkit, Py-ART, is an open source Python module containing
a growing collection of weather radar algorithms and utilities build on top of
the Scientific Python stack and distributed under the 3-Clause BSD license.
Py-ART is used by the Atmospheric Radiation Measurement (ARM) Climate Research
Facility for working with data from a number of precipitation and cloud radars,
but has been designed so that it can be used by others in the radar and
atmospheric communities to examine, processes, and analyse data from many types
of weather radars.
"""


DOCLINES = __doc__.split("\n")

import os
from os import path
import sys
import subprocess
import glob
from numpy import get_include

from setuptools import find_packages, setup, Extension
from Cython.Build import cythonize
import Cython

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

CLASSIFIERS = """\
    Development Status :: 5 - Production/Stable
    Intended Audience :: Science/Research
    Intended Audience :: Developers
    License :: OSI Approved :: BSD License
    Programming Language :: Python
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3.6
    Programming Language :: Python :: 3.7
    Programming Language :: Python :: 3.8
    Programming Language :: Python :: 3.9
    Programming Language :: Python :: 3.10
    Programming Language :: C
    Programming Language :: Cython
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Atmospheric Science
    Operating System :: POSIX :: Linux
    Operating System :: MacOS :: MacOS X
    Operating System :: Microsoft :: Windows
    Framework :: Matplotlib
"""

NAME = 'arm_pyart'
AUTHOR = "Scott Collis, Jonathan Helmus"
AUTHOR_EMAIL = "scollis@anl.gov"
MAINTAINER = "Py-ART Developers"
MAINTAINER_EMAIL = "zsherman@anl.gov, scollis@anl.gov, mgrover@anl.gov"
DESCRIPTION = DOCLINES[0]
LONG_DESCRIPTION = "\n".join(DOCLINES[2:])
URL = "https://github.com/ARM-DOE/pyart"
DOWNLOAD_URL = "https://github.com/ARM-DOE/pyart"
LICENSE = 'BSD'
CLASSIFIERS = list(filter(None, CLASSIFIERS.split("\n")))
PLATFORMS = ["Linux", "Mac OS-X", "Unix", "Windows"]
SCRIPTS = glob.glob('scripts/*')


# This is a bit hackish: we are setting a global variable so that the main
# pyart __init__ can detect if it is being loaded by the setup routine, to
# avoid attempting to load components that aren't built yet. While ugly, it's
# a lot more robust than what was previously being used.


# NOTE: This file must remain Python 2 compatible for the foreseeable future,
# to ensure that we error out properly for people with outdated setuptools
# and/or pip.
min_version = (3, 6)
if sys.version_info < min_version:
    error = """
act does not support Python {}.{}.
Python {}.{} and above is required. Check your Python version like so:
python3 --version
This may be due to an out-of-date pip. Make sure you have pip >= 9.0.1.
Upgrade pip like so:
pip install --upgrade pip
""".format(
        *sys.version_info[:2], *min_version
    )
    sys.exit(error)

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as readme_file:
    readme = readme_file.read()

with open(path.join(here, 'requirements.txt')) as requirements_file:
    # Parse requirements.txt, ignoring any commented-out lines.
    requirements = [
        line for line in requirements_file.read().splitlines() if not line.startswith('#')
    ]


extensions = []

# RSL Path if present
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

rsl_path = os.environ.get('RSL_PATH')
if rsl_path is None:
    rsl_path = guess_rsl_path()
rsl_lib_path = os.path.join(rsl_path, 'lib')
rsl_include_path = os.path.join(rsl_path, 'include')

# build the RSL IO and FourDD dealiaser if RSL is installed
if check_rsl_path(rsl_lib_path, rsl_include_path):
    fourdd_sources = [
        'pyart/correct/src/dealias_fourdd.c',
        'pyart/correct/src/sounding_to_volume.c',
        'pyart/correct/src/helpers.c'
  ]

    # Cython wrapper around FourDD
    extension_4dd = Extension(
        'pyart.correct._fourdd_interface',
        sources=['pyart/correct/_fourdd_interface.pyx',] + fourdd_sources,
        libraries=['rsl'],
        library_dirs=[rsl_lib_path],
        include_dirs=[
            rsl_include_path, 'pyart/correct/src'] + [get_include()],
        runtime_library_dirs=[rsl_lib_path])

    # Cython wrapper around RSL io
    extension_rsl = Extension(
        'pyart.io._rsl_interface',
        sources=['pyart/io/_rsl_interface.pyx'],
        libraries=['rsl'],
        library_dirs=[rsl_lib_path],
        include_dirs=[
            rsl_include_path] + [get_include()],
        runtime_library_dirs=[rsl_lib_path],)

    extensions.append(extension_rsl)
    extensions.append(extension_4dd)
else:
    import warnings
    warnings.warn(RSL_MISSING_WARNING % (rsl_path))

libraries = []
if os.name == 'posix':
    libraries.append('m')

# Check build pyx extensions
extension_check_build = Extension(
    'pyart.__check_build._check_build', sources=['pyart/__check_build/_check_build.pyx'],
    include_dirs=[get_include()])

extensions.append(extension_check_build)

# Correct pyx extensions
extension_edge_finder = Extension(
    'pyart.correct._fast_edge_finder', sources=['pyart/correct/_fast_edge_finder.pyx'],
    include_dirs=[get_include()])

extension_1d = Extension(
    'pyart.correct._unwrap_1d', sources=['pyart/correct/_unwrap_1d.pyx'],
    include_dirs=[get_include()])

unwrap_sources_2d = [
    'pyart/correct/_unwrap_2d.pyx', 'pyart/correct/unwrap_2d_ljmu.c']
extension_2d = Extension('pyart.correct._unwrap_2d', sources=unwrap_sources_2d,
                         include_dirs=[get_include()])

unwrap_sources_3d = [
    'pyart/correct/_unwrap_3d.pyx', 'pyart/correct/unwrap_3d_ljmu.c']
extension_3d = Extension('pyart.correct._unwrap_3d', sources=unwrap_sources_3d,
                         include_dirs=[get_include()])

extensions.append(extension_edge_finder)
extensions.append(extension_1d)
extensions.append(extension_2d)
extensions.append(extension_3d)

# IO pyx extensions
extension_sigmet = Extension(
    'pyart.io._sigmetfile', sources=['pyart/io/_sigmetfile.pyx'],
    include_dirs=[get_include()])

extension_nexrad = Extension(
    'pyart.io.nexrad_interpolate', sources=['pyart/io/nexrad_interpolate.pyx'],
    include_dirs=[get_include()])

extensions.append(extension_sigmet)
extensions.append(extension_nexrad)

# Map pyx extensions
extension_ckd = Extension(
    'pyart.map.ckdtree', sources=['pyart/map/ckdtree.pyx'],
    include_dirs=[get_include()],
    libraries=libraries)

extension_load_nn = Extension(
    'pyart.map._load_nn_field_data', sources=['pyart/map/_load_nn_field_data.pyx'],
    include_dirs=[get_include()])

extension_gate_to_grid = Extension(
    'pyart.map._gate_to_grid_map', sources=['pyart/map/_gate_to_grid_map.pyx'],
    libraries=libraries)

extensions.append(extension_ckd)
extensions.append(extension_load_nn)
extensions.append(extension_gate_to_grid)


# Retrieve pyx extensions
extension_kdp = Extension(
    'pyart.retrieve._kdp_proc', sources=['pyart/retrieve/_kdp_proc.pyx'])

extensions.append(extension_kdp)

setup(
    name='arm_pyart',
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    maintainer=MAINTAINER,
    maintainer_email=MAINTAINER_EMAIL,
    url=URL,
    packages=find_packages(include=['pyart'], exclude=['docs']),
    include_package_data=True,
    scripts=SCRIPTS,
    install_requires=requirements,
    license=LICENSE,
    platforms=PLATFORMS,
    classifiers=CLASSIFIERS,
    zip_safe=False,
    use_scm_version={
        'version_scheme': 'post-release',
        'local_scheme': 'dirty-tag',
    },
    ext_modules=cythonize(
        extensions, compiler_directives={'language_level' : "3"}),
)
