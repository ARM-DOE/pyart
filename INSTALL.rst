============
Installation
============

Required Dependencies
=====================

Py-ART requires the following software.

* Python__ 3.11x, 3.12x or 3.13x

__ http://www.python.org

* NumPy__

__ http://www.scipy.org

* SciPy__

__ http://www.scipy.org

* matplotlib__

__ http://matplotlib.org/

* netCDF4__

__ https://github.com/Unidata/netcdf4-python

* Xarray__

__ https://docs.xarray.dev/en/stable/

* pooch__

__ https://pypi.org/project/pooch/

* Cython__

__ https://cython.readthedocs.io/en/latest/

* setuptools__

__ https://setuptools.pypa.io/en/latest/index.html

* cartopy__

__ https://scitools.org.uk/cartopy/docs/latest/

* cmweather__

__ https://cmweather.readthedocs.io/en/latest/

* xradar__

__ https://docs.openradarscience.org/projects/xradar/en/stable/

* xarray__

__ https://docs.xarray.dev/en/stable/

* mda-xdrlib__

__ https://github.com/MDAnalysis/mda-xdrlib

* fsspec__

__ https://filesystem-spec.readthedocs.io/en/latest/

* s3fs__

__ https://s3fs.readthedocs.io/en/latest/

* pandas__

__ https://pandas.pydata.org/

Optional Dependencies
=====================

The following packages are recommended for a fully-functional Py-ART
installation, but the package will install and work with reduced functionality
without these packages.

* `CyLP <https://github.com/mpy/CyLP>`_ or
  `PyGLPK <https://tfinley.net/software/pyglpk/>`_ or
  `CVXOPT <https://cvxopt.org/>`_ and their dependencies.

* `Cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_ or
* `Basemap <https://matplotlib.org/basemap/>`_ But Cartopy is recommended as
  basemap will no longer have support.

* `pyproj <https://github.com/pyproj4/pyproj>`_

* `pytest <https://docs.pytest.org/en/latest/>`_

* `metpy <https://unidata.github.io/MetPy/latest/>`_

Instructions for Installing
===========================

You can install using pip, conda, or pixi! We recommend using pixi, using the instructions below

1. Install pixi if you have not done so already https://pixi.sh/latest/installation/

2. Create a Py-ART environment using the following commands

::

    $ pixi init pyart_dev
    $ cd pyart_dev
    $ pixi add arm_pyart

This will build your environment and install all the required dependencies!


Obtaining the latest source
===========================

The latest source code for Py-ART can be obtained from the GitHub repository,
https://github.com/ARM-DOE/pyart.

The latest source can be checked out using

::

    $ git clone https://github.com/ARM-DOE/pyart.git

Installing from Source
======================

NOTE: TRMM RSL is deprecated, please consider using radx to convert files such as DORADE.

The path to the TRMM RSL library must be provided during install. This can
either be done by setting the ``RSL_PATH`` environmentation variable. In bash
this can be done using ``export RSL_PATH=/path/to/rsl/``. If this location is
not specified, some common locations will be searched. Note that the location
provided should be the root TRMM RSL path, under which both a `lib` and
`include` directory are contained, the default location is ``/lib/local/trmm``.
If using CyLP, a path for the coincbc directory is needed. This can be done
using ``export COIN_INSTALL_DIR=/path/to/coincbc/``. When using CyLP, on some
systems, installing the Anaconda compilers is needed. These can be found here:
https://docs.conda.io/projects/conda-build/en/latest/resources/compiler-tools.html

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

* I'm getting a segfault or compile error with CyLP in newer Python versions
  when installing in an environment.

  Anaconda has its own compilers now on conda-forge. Theres can be found here:
  https://docs.conda.io/projects/conda-build/en/latest/resources/compiler-tools.html
  Once the proper compilers are installed, reinstall CyLP.

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

* When running basemap, I get an error 'KeyError: PROJ_LIB'.

  Basemap is not being supported beyond 2020, some of these errors relate
  to it not playing nicely with newer versions of packages. We recommend using
  Cartopy instead, but some users have been able to use:
  import os
  os.environ['PROJ_LIB'] = 'C:/Users/xx Username xxx/Anaconda3/Lib/site-packages/mpl_toolkits/basemap'
  To get basemap working, but again Cartopy should be used instead of Basemap.
