.. -*- mode: rst -*-
The Python ARM Radar Toolkit (Py-ART)
=====================================

|Travis| |AppVeyor|

|AnacondaCloud| |CondaDownloads|

|DocsDev| |DocsUsers| |DocsGuides|

|ARM| |Tweet|

.. |Travis| image:: https://api.travis-ci.org/ARM-DOE/pyart.png?branch=master
    :target: https://travis-ci.org/ARM-DOE/pyart

.. |AppVeyor| image:: https://ci.appveyor.com/api/projects/status/9do57qycha65j4v9/branch/master?svg=true
    :target: https://ci.appveyor.com/project/JonathanHelmus/pyart-l711v/branch/master

.. |AnacondaCloud| image:: https://anaconda.org/conda-forge/arm_pyart/badges/version.svg
    :target: https://anaconda.org/conda-forge/arm_pyart

.. |CondaDownloads| image:: https://anaconda.org/conda-forge/arm_pyart/badges/downloads.svg
    :target: https://anaconda.org/conda-forge/arm_pyart/files

.. |DocsDev| image:: https://img.shields.io/badge/docs-developers-4088b8.svg
    :target: http://arm-doe.github.io/pyart-docs-travis/dev_reference/index.html

.. |DocsUsers| image:: https://img.shields.io/badge/docs-users-4088b8.svg
    :target: http://arm-doe.github.io/pyart-docs-travis/user_reference/index.html

.. |DocsGuides| image:: https://img.shields.io/badge/docs-guides-4088b8.svg
    :target: https://github.com/ARM-DOE/pyart/tree/master/guides/

.. |ARM| image:: https://img.shields.io/badge/Sponsor-ARM-blue.svg?colorA=00c1de&colorB=00539c
    :target: https://www.arm.gov/

.. |Tweet| image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social
    :target: https://twitter.com/Py_ART

The Python ARM Radar Toolkit, Py-ART, is an open source Python module 
containing a growing collection of weather radar algorithms and utilities
build on top of the Scientific Python stack and distributed under the
3-Clause BSD license. Py-ART is used by the 
`Atmospheric Radiation Measurement (ARM) Climate Research Facility 
<http://www.arm.gov>`_ for working with data from a number of precipitation
and cloud radars, but has been designed so that it can be used by others in
the radar and atmospheric communities to examine, processes, and analyze
data from many types of weather radars. 


Important Links
===============

- Official source code repository: https://github.com/ARM-DOE/pyart
- HTML documentation: http://arm-doe.github.io/pyart-docs-travis/
- Examples: http://arm-doe.github.io/pyart/dev/auto_examples/index.html
- Mailing List: http://groups.google.com/group/pyart-users/
- Issue Tracker: https://github.com/ARM-DOE/pyart/issues


Citing
======

If you use the Python ARM Radar Toolkit (Py-ART) to prepare a publication
please cite:

    Helmus, J.J. & Collis, S.M., (2016). The Python ARM Radar Toolkit
    (Py-ART), a Library for Working with Weather Radar Data in the Python
    Programming Language. Journal of Open Research Software. 4(1), p.e25.
    DOI: http://doi.org/10.5334/jors.119

Py-ART implements many published scientific methods which should *also* be
cited if you make use of them. Refer to the **References** section in the
documentation of the functions used for information on these citations.


Install
=======

The easiest method for installing Py-ART is to use the conda packages from
the latest release. To do this you must download and install 
`Anaconda <https://www.anaconda.com/download/#>`_ or 
`Miniconda <https://conda.io/miniconda.html>`_.
With Anaconda or Miniconda install, it is recommended to create a new conda
environment when using Py-ART or even other packages. To create a new
environment based on the `environment.yml <https://github.com/ARM-DOE/pyart/blob/master/environment.yml>`_::

    conda env create -f environment.yml

Or for a basic environment and downloading optional dependencies as needed::

    conda create -n pyart_env -c conda-forge python=3.7 arm_pyart

Basic command in a terminal or command prompt to install the latest version of
Py-ART::

    conda install -c conda-forge arm_pyart

To update an older version of Py-ART to the latest release use::

    conda update -c conda-forge arm_pyart

If you do not wish to use Anaconda or Miniconda as a Python environment or want
to use the latest, unreleased version of Py-ART see the section below on 
**Installing from source**.


Configuration
=============

The configuration file in Py-ART specifies the default metadata, field names,
colormaps and plot limits. A custom configuration can be loaded
automatically be setting the environmental variable **PYART_CONFIG** to point
to a custom configuration file. For additional details on this process see the
documentation on the `pyart.load_config` function.


Extensions and related software
===============================

A number of projects are available which extend the functionality of Py-ART.
These include:

* `ARTView <https://github.com/nguy/artview>`_ : 
  Interactive radar viewing browser.

* `PyTDA <https://github.com/nasa/PyTDA>`_ : 
  Python Turbulence Detection Algorithm.

* `SingleDop <https://github.com/nasa/SingleDop>`_ : 
  Single Doppler Retrieval Toolkit.

* `DualPol <https://github.com/nasa/DualPol>`_ :
  Python Interface to Dual-Pol Radar Algorithms.

* `PyBlock <https://github.com/nasa/PyBlock>`_:
  Python Polarimetric Radar Beam Blockage Calculation


Other related open source software for working with weather radar data:

* `wradlib <https://wradlib.org>`_ :
  An open source library for weather radar data processing.
  
* `BALTRAD <https://baltrad.eu/>`_ : Community-based weather radar networking.

* `MMM-Py <https://github.com/nasa/MMM-Py>`_ : 
  Marshall MRMS Mosaic Python Toolkit.

* `CSU_RadarTools <https://github.com/CSU-Radarmet/CSU_RadarTools>`_ : 
  Colorado State University Radar Tools.

* `TRMM RSL <https://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/>`_ :
  TRMM Radar Software Library.

* `RadX <https://www.ral.ucar.edu/projects/titan/docs/radial_formats/radx.html>`_: 
  Radx C++ Software Package for Radial Radar Data.


Dependencies
============

Py-ART is tested to work under Python 2.7, 3.4, 3.5, and 3.6.

The required dependencies to install Py-ART in addition to Python are:

* `NumPy <https://www.numpy.org/>`_
* `SciPy <https://www.scipy.org>`_
* `matplotlib <https://matplotlib.org/>`_
* `netCDF4 <https://github.com/Unidata/netcdf4-python>`_

A working C/C++ compiler is required for some optional modules. An easy method
to install these dependencies is by using a 
`Scientific Python distributions <http://scipy.org/install.html>`_.
`Anaconda <https://www.anaconda.com/distribution/>`_ will install all of
the above packages by default on Windows, Linux and Mac computers and is
provided free of charge by Anaconda.


Optional Dependences
====================

The above Python modules are require before installing Py-ART, additional
functionality is available of the following modules are installed.

* `TRMM Radar Software Library (RSL) 
  <https://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/>`_.  
  If installed Py-ART will be able to read in radar data in a number of 
  additional formats (Lassen, McGill, Universal Format, and RADTEC) and 
  perform automatic dealiasing of Doppler velocities.  RSL should be
  install prior to installing Py-ART. The environmental variable `RSL_PATH`
  should point to the location where RSL was installed if RSL was not
  installed in the default location (/usr/local/trmm), such as a anaconda path
  (/usr/anaconda3/envs/pyart_env/.

* In order to read files which are stored in HDF5 files the
  `h5py <https://www.h5py.org/>`_ package and related libraries must be
  installed.

* A linear programming solver and Python wrapper to use the LP phase
  processing method. `CyLP <https://github.com/mpy/CyLP>`_ is recommended as
  it gives the fastest results, but 
  `PyGLPK <https://tfinley.net/software/pyglpk/>`_ and 
  `CVXOPT <https://cvxopt.org/>`_ are also supported. The underlying LP 
  solvers `CBC <https://projects.coin-or.org/Cbc>`_ or 
  `GLPK <https://www.gnu.org/software/glpk/>`_ will also be required depending
  on which wrapper is used. When using `CyLP <https://github.com/mpy/CyLP>`_
  a path to coincbc is needed by setting the `COIN_INSTALL_DIR` path, such as
  (/usr/anaconda3/envs/pyart_env/).

* `Cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_. If installed,
  the ability to plot grids on geographic maps is available.

* `xarray <https://xarray.pydata.org/en/stable/>`_. If installed, gives the
  ability to work with the grid dataset used in grid plotting.

* `Basemap <https://matplotlib.org/basemap/>`_. If installed, also gives the
  ability to plot grids on geographic maps, but Cartopy is recommended over
  Basemap.

* `wradlib <https://docs.wradlib.org/en/latest/>`_.  Needed to calculate the texture
  of a differential phase field.

* `pytest <https://docs.pytest.org/en/latest/>`_.
  Required to run the Py-ART unit tests.

* `gdal <https://pypi.python.org/pypi/GDAL/>`_.
  Required to output GeoTIFFs from `Grid` objects.
 
Installing from source
======================

Installing Py-ART from source is the only way to get the latest updates and
enhancement to the software that have not yet made it into a release.
The latest source code for Py-ART can be obtained from the GitHub repository,
https://github.com/ARM-DOE/pyart. Either download and unpack the 
`zip file <https://github.com/ARM-DOE/pyart/archive/master.zip>`_ of 
the source code or use git to checkout the repository::

    git clone https://github.com/ARM-DOE/pyart.git

To install in your home directory, use::

    python setup.py install --user

To install for all users on Unix/Linux::

    python setup.py build
    sudo python setup.py install


Development
===========

Py-ART is an open source, community software project. Contributions to
the package are welcomed from all users.

Code
----
The latest source code can be obtained with the command::
 
    git clone https://github.com/ARM-DOE/pyart.git

If you are planning on making changes that you would like included in Py-ART,
forking the repository is highly recommended.

Contributing
-------------

We welcome contributions for all uses of Py-ART, provided the code can be
distributed under the BSD 3-clause license. A copy of this license is
available in the **LICENSE.txt** file in this directory. For more on
contributing, see the `contributor's guide. <https://github.com/ARM-DOE/pyart/blob/master/guides/contributors_guide.rst>`_

Testing
-------

After installation, you can launch the test suite from outside the
source directory (you will need to have pytest installed)::

   $ pytest --pyargs pyart

In-place installs can be tested using the `pytest` command from within
the source directory.
