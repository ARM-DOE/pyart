==============================
Building and installing Py-ART
==============================

Required Dependencies
=====================

Py-ART requires the following software.

* Python__ 2.7.x, 3.4.x or 3.5.x

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

* TRMM RSL__

__ http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/

* `CyLP <https://github.com/mpy/CyLP>`_ or 
  `PyGLPK <http://tfinley.net/software/pyglpk/>`_ or
  `CVXOPT <http://cvxopt.org/>`_ and their dependencies.

* Basemap__

__ http://matplotlib.org/basemap/

* pyproj__

__ http://code.google.com/p/pyproj/

* nose__

__ http://nose.readthedocs.org/en/latest/

Obtaining the latest source
===========================

The latest source code for Py-ART can be obtained from the GitHub repository,
https://github.com/ARM-DOE/pyart.

The latest source can be checked out using

    $ git clone https://github.com/ARM-DOE/pyart.git


Installing from Source
======================

The path to the TRMM RSL library must be provided during install.  This can
either be done by setting the `RSL_PATH` environmentation variable.  In bash
this can be done using `export RSL_PATH=/path/to/rsl/`. If this location is
not specified, some common locations will be searched.  Note that the location
provided should be the root TRMM RSL path, under which both a `lib` and
`include` directory are contained, the default location is `/lib/local/trmm`.

After specifying the TRMM RSL path Py-ART can be installed globally using

    $ python setup.py install

of locally using

    $ python setup.py install --user

If you prefer to use Py-ART without installing, simply add the this path to
your PYTHONPATH (directory or with a .pth file) and compile the extension
in-place.

    $ python setup.py build_ext -i

Frequently asked questions
==========================

* I'm getting a segfault or another error in python when using 
  pyart.io.read_rsl() with IRIS/other files.
  
  This is due to a bug in RSL, and can be remedied by adding
  -fno-stack-protector -D_FORTIFY_SOURCE=0 to the CFLAGS parameter of the
  makefile of RSL.  This issue has been fixed with the release of rsl-v1.44.

* I'm having trouble getting PyGLPK to compile on my 64-bit operating system.
  Change the line in the setup.py file from
  
  define_macros = macros, extra_compile_args=['-m32'], extra_link_args=['-m32'],
  
  to
  
  define_macros = macros, extra_compile_args=['-m64'], extra_link_args=['-m64'],

  Then build and install PyGLPK as recommended in the PYGLPK README.txt file.
