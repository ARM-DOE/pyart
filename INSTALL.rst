==============================
Building and installing Py-ART
==============================

Required Dependencies
=====================

Py-ART requires the following software.

* Python__ 2.7.x, 3.6.x or 3.7.x

__ http://www.python.org

* NumPy__

__ http://www.scipy.org

* SciPy__

__ http://www.scipy.org

* matplotlib__

__ http://matplotlib.org/

* netCDF4__

__ https://github.com/Unidata/netcdf4-python


Optional Dependencies
=====================

The following packages are recommended for a fully-functional Py-ART
installation, but the package will install and work with reduced functionality
without these packages.

* `TRMM RSL <https://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/>`_

* `CyLP <https://github.com/mpy/CyLP>`_ or 
  `PyGLPK <https://tfinley.net/software/pyglpk/>`_ or
  `CVXOPT <https://cvxopt.org/>`_ and their dependencies.

* `Cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_ or
* `Basemap <https://matplotlib.org/basemap/>`_ But Cartopy is recommended as
  basemap will no longer have support.

* `xarray <https://xarray.pydata.org/en/stable/`_
* `pyproj <https://code.google.com/p/pyproj/`_

* `pytest <https://docs.pytest.org/en/latest/`_

Obtaining the latest source
===========================

The latest source code for Py-ART can be obtained from the GitHub repository,
https://github.com/ARM-DOE/pyart.

The latest source can be checked out using

::

    $ git clone https://github.com/ARM-DOE/pyart.git


Installing from Source
======================

The path to the TRMM RSL library must be provided during install. This can
either be done by setting the ``RSL_PATH`` environmentation variable. In bash
this can be done using ``export RSL_PATH=/path/to/rsl/``. If this location is
not specified, some common locations will be searched. Note that the location
provided should be the root TRMM RSL path, under which both a `lib` and
`include` directory are contained, the default location is ``/lib/local/trmm``.
If using CyLP, a path for the coincbc directory is needed. This can be done
using ``export COIN_INSTALL_DIR=/path/to/coincbc/``.

After specifying the TRMM RSL path Py-ART can be installed globally using

::

    $ python setup.py install

of locally using

::

    $ python setup.py install --user

If you prefer to use Py-ART without installing, simply add the this path to
your ``PYTHONPATH`` (directory or with a .pth file) and compile the extension
in-place.

::

    $ python setup.py build_ext -i

You can also install Py-ART in development mode by using

::

    $ pip install -e .

Frequently asked questions
==========================

* I'm getting a no 'io' module after installing pyart with pip.

  There is a pyart on pip that is a different package. Make sure to do::

      pip install arm_pyart

  and not::

      pip install pyart

* I'm getting a segfault or another error in python when using 
  ``pyart.io.read_rsl()`` with IRIS/other files.
  
  This is due to a bug in RSL, and can be remedied by adding
  ``-fno-stack-protector -D_FORTIFY_SOURCE=0`` to the CFLAGS parameter of the
  makefile of RSL.  This issue has been fixed with the release of rsl-v1.44.

* I'm having trouble getting PyGLPK to compile on my 64-bit operating system.
  
  Change the line in the setup.py file from
  
  ::
  
      define_macros = macros, extra_compile_args=['-m32'], extra_link_args=['-m32'],
  
  to
  
  ::
  
      define_macros = macros, extra_compile_args=['-m64'], extra_link_args=['-m64'],

  Then build and install PyGLPK as recommended in the PYGLPK README.txt file.
