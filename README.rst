.. -*- mode: rst -*-

|Travis|_

.. |Travis| image:: https://api.travis-ci.org/ARM-DOE/pyart.png?branch=master
.. _Travis: https://travis-ci.org/ARM-DOE/pyart

The Python ARM Radar Toolkit (Py-ART)
=====================================

The Python ARM Radar Toolkit, Py-ART, is an open source Python module 
containing a growing collection of weather radar algorithms and utilities
build on top of the Scientific Python stack and distributed under the
3-Clause BSD license. Py-ART is used by the 
`Atmospheric Radiation Measurement (ARM) Climate Research Facility 
<http://www.arm.gov>`_ for working with data from a number of precipitation
and cloud radars, but has been designed so that it can be used by others in
the radar and atmospheric communities to examine, processes, and analyse
data from many types of weather radars. 

Important Links
===============

- Official source code repository: https://github.com/ARM-DOE/pyart
- HTML documentation: http://arm-doe.github.io/pyart-docs-travis/
- Examples: http://arm-doe.github.io/pyart/dev/auto_examples/index.html
- Mailing List: http://groups.google.com/group/pyart-users/
- Issue Tracker: https://github.com/ARM-DOE/pyart/issues

Dependencies
============

Py-ART is tested to work under Python 2.6, 2.7, 3.3, and 3.4.

The required dependencies to install Py-ART in addition to Python are:

* `NumPy <http://www.scipy.org>`_ 1.6+
* `SciPy <http://www.scipy.org>`_ 0.10+
* `matplotlib <http://matplotlib.org/>`_ 1.1.0+
* `netCDF4 <https://github.com/Unidata/netcdf4-python>`_ 1.0.2+ 

As well as a working C/C++ compiler.  An easy method to install these
dependencies is by using a 
`Scientific Python distributions <http://scipy.org/install.html>`_.
`Anaconda <https://store.continuum.io/cshop/anaconda/>`_ will install all of
the above packages by default on Linux and Mac computers and is provided
free of charge by Continuum Analytics.

Optional Dependences
====================

The above Python modules are require before installing Py-ART, additional
functionality is available of the the following modules are installed.

* `TRMM Radar Software Library (RSL) 
  <http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/>`_.  
  If installed Py-ART will be able to read in radar data in a number of 
  additional formats (Lassen, McGill, Universal Format, and RADTEC) and 
  perform automatic dealiasing of doppler velocities.  RSL should be
  install prior to installing Py-ART. The environmental variable `RSL_PATH`
  should point to the location where RSL was installed if RSL was not
  installed in the default location (/usr/local/trmm).

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
 
Install
=======

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

Additional install instructions can be found in the **INSTALL.rst** file in
this directory.

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
