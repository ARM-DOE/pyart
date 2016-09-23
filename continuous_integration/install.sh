#!/bin/bash
# This script is adapted from the install.sh script from the scikit-learn
# project: https://github.com/scikit-learn/scikit-learn

# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

set -e
# use next line to debug this script
#set -x

# Use Miniconda to provide a Python environment.  This allows us to perform
# a conda based install of the SciPy stack on multiple versions of Python
# as well as use conda and binstar to install additional modules which are not
# in the default repository.
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh \
    -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b
export PATH=/home/travis/miniconda2/bin:$PATH
conda update --yes conda
conda update --yes conda

# Create a testenv with the correct Python version
conda create -n testenv --yes pip python=$PYTHON_VERSION
source activate testenv

# Install Py-ART dependencies
conda install --yes numpy scipy matplotlib netcdf4 nose
conda install --yes basemap
conda install --yes -c jjhelmus trmm_rsl

if [[ $PYTHON_VERSION == '2.7' ]]; then
    conda install --yes -c http://conda.anaconda.org/jjhelmus cbc cylp
    conda install --yes -c http://conda.anaconda.org/jjhelmus glpk pyglpk
    conda install --yes -c http://conda.anaconda.org/jjhelmus cvxopt_glpk

    # wradlib and dependencies
    conda install --yes h5py
    # KLUDGE libgdal does not report its version dependency on geos which
    # causes either gdal or basemap to break, force the exact libgdal version
    # see: https://github.com/ContinuumIO/anaconda-issues/issues/584
    conda install --yes gdal basemap libgdal=2.0.0=0 krb5 proj4
    conda install --no-deps --yes -c conda-forge wradlib
fi

# install coverage modules
pip install nose-cov
if [[ "$COVERALLS" == "true" ]]; then
    pip install python-coveralls
fi

# install Py-ART
export RSL_PATH=~/miniconda2/envs/testenv

if [[ "$FROM_RECIPE" == "true" ]]; then
    source deactivate
    conda install --yes conda-build
    conda install --yes jinja2 setuptools
    conda config --add channels http://conda.anaconda.org/jjhelmus
    source activate testenv
    conda build --no-test -q conda_recipe/
   
    export CONDA_PACKAGE=`conda build --output conda_recipe/ | grep bz2`
    conda install --yes $CONDA_PACKAGE
    mkdir foo   # required so source directory not picked up during tests
    cd foo
else
    python setup.py build_ext --inplace
fi

# KLUDGE
# cylp and cvxopt_glpk depend on BLAS and LAPACK which are provided by the
# system and depend on the system libgfortran.  The conda libgfortran does not
# export the symbols required for the system packages, so it must be removed.
conda install --yes libgfortran
conda remove --yes --force libgfortran
