.. -*- mode: rst -*-

|Travis|_ |AppVeyor|_

.. |Travis| image:: https://api.travis-ci.org/ARM-DOE/pyart.png?branch=master
.. _Travis: https://travis-ci.org/ARM-DOE/pyart

.. |AppVeyor| image:: https://ci.appveyor.com/api/projects/status/9do57qycha65j4v9/branch/master?svg=true
.. _AppVeyor: https://ci.appveyor.com/project/JonathanHelmus/pyart-l711v/branch/master


The Python ARM Radar Toolkit (Py-ART)
=====================================

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
cited if you make use of them.  Refer to the **References** section in the
documentation of the functions used for information on these citations.


Install
=======

The easiest method for installing Py-ART is to use the conda packages from
the latest release.  To do this you must download and install 
`Anaconda <http://continuum.io/downloads>`_ or 
`Miniconda <http://continuum.io/downloads>`_.  
Then use the following command in a terminal or command prompt to install
the latest version of Py-ART::

    conda install -c https://conda.anaconda.org/jjhelmus pyart

To update an older version of Py-ART to the latest release use::

    conda update -c https://conda.anaconda.org/jjhelmus pyart

If you do not wish to use Anaconda or Miniconda as a Python environment or want
to use the latest, unreleased version of Py-ART see the section below on 
**Installing from source**.


Configuration
=============

The configuration file in Py-ART specifies the default metadata, field names,
colormaps and plot limits.  A custom configuration can be loaded
automatically be setting the environmental variable **PYART_CONFIG** to point
to a custom configuration file.  For additional details on this process see the
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

* `wradlib <http://wradlib.bitbucket.org/>`_ : 
  An open source library for weather radar data processing.
  
* `BALTRAD <http://baltrad.eu/>`_ : Community-based weather radar networking.

* `MMM-Py <https://github.com/nasa/MMM-Py>`_ : 
  Marshall MRMS Mosaic Python Toolkit.

* `CSU_RadarTools <https://github.com/CSU-Radarmet/CSU_RadarTools>`_ : 
  Colorado State University Radar Tools.

* `TRMM RSL <http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/>`_ :
  TRMM Radar Software Library.

* `RadX <http://www.ral.ucar.edu/projects/titan/docs/radial_formats/radx.html>`_: 
  Radx C++ Software Package for Radial Radar Data.


Dependencies
============

Py-ART is tested to work under Python 2.7, 3.4, and 3.5.

The required dependencies to install Py-ART in addition to Python are:

* `NumPy <http://www.scipy.org>`_
* `SciPy <http://www.scipy.org>`_
* `matplotlib <http://matplotlib.org/>`_
* `netCDF4 <https://github.com/Unidata/netcdf4-python>`_

As well as a working C/C++ compiler.  A Fortran compiler is required for some
optional modules. An easy method to install these dependencies is by using a 
`Scientific Python distributions <http://scipy.org/install.html>`_.
`Anaconda <https://store.continuum.io/cshop/anaconda/>`_ will install all of
the above packages by default on Windows, Linux and Mac computers and is
provided free of charge by Continuum Analytics.


Optional Dependences
====================

The above Python modules are require before installing Py-ART, additional
functionality is available of the following modules are installed.

* `TRMM Radar Software Library (RSL) 
  <http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/>`_.  
  If installed Py-ART will be able to read in radar data in a number of 
  additional formats (Lassen, McGill, Universal Format, and RADTEC) and 
  perform automatic dealiasing of Doppler velocities.  RSL should be
  install prior to installing Py-ART. The environmental variable `RSL_PATH`
  should point to the location where RSL was installed if RSL was not
  installed in the default location (/usr/local/trmm).

* In order to read files which are stored in HDF5 files the
  `h5py <http://www.h5py.org/>`_ package and related libraries must be
  installed.

* A linear programming solver and Python wrapper to use the LP phase
  processing method. `CyLP <https://github.com/mpy/CyLP>`_ is recommended as
  it gives the fastest results, but 
  `PyGLPK <http://tfinley.net/software/pyglpk/>`_ and 
  `CVXOPT <http://cvxopt.org/>`_ are also supported. The underlying LP 
  solvers `CBC <https://projects.coin-or.org/Cbc>`_ or 
  `GLPK <http://www.gnu.org/software/glpk/>`_ will also be required depending
  on which wrapper is used.

* `Basemap <http://matplotlib.org/basemap/>`_. If installed the ability to 
  plot grids on geographic maps is available.

* `wradlib <http://wradlib.bitbucket.org/>`_.  Needed to calculate the texture
  of a differential phase field.

* `nose <http://nose.readthedocs.org/en/latest/>`_.  
  Required to run the Py-ART unit tests.

* `gdal <https://pypi.python.org/pypi/GDAL/>`_.
  Required to output GeoTIFFs from `Grid` objects.
 
Installing from source
======================

Installing Py-ART from source is the only way to get the latest updates and
enhancement to the software that have not yet made it into a release.
The latest source code for Py-ART can be obtained from the GitHub repository,
https://github.com/ARM-DOE/pyart.  Either download and unpack the 
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

Py-ART is an open source, community software project.  Contributions to
the package are welcomed from all users.  

Code
----
The latest source code can be obtained with the command::
 
    git clone https://github.com/ARM-DOE/pyart.git

If you are planning on making changes that you would like included in Py-ART,
forking the repository is highly recommended.

Contributing
-------------

We welcome contributions for all used of Py-ART provided the code can be
distributed under the BSD 3-clause license.  A copy of this license is
available in the **LICENSE.txt** file in this directory.  

Testing
-------

After installation, you can launch the test suite from outside the
source directory (you will need to have nosetests installed)::

   $ nosetests --exe pyart

In-place installs can be tested using the `nosetest` command from within
the source directory.
