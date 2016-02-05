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
conda install --yes -c http://conda.anaconda.org/jjhelmus trmm_rsl

if [[ $PYTHON_VERSION == '2.7' ]]; then
    conda install --yes basemap 
    conda install --yes -c http://conda.anaconda.org/jjhelmus cbc cylp
    conda install --yes -c http://conda.anaconda.org/jjhelmus glpk pyglpk
    conda install --yes -c http://conda.anaconda.org/jjhelmus cvxopt_glpk

    # wradlib and dependencies
    # KLUDGE libgdal does not report its version dependency on geos which
    # causes either gdal or basemap to break, force the exact libgdal version
    # see: https://github.com/ContinuumIO/anaconda-issues/issues/584
    conda install --yes sphinx gdal numpydoc h5py basemap libgdal=2.0.0=0
    conda install --yes sphinx_rtd_theme
    pip install sphinxcontrib-bibtex
    pip install xmltodict
    pip install wradlib
fi
if [[ $PYTHON_VERSION == '3.4' ]]; then
    conda install --yes basemap 
fi
if [[ $PYTHON_VERSION == '3.5' ]]; then
    conda install --yes basemap 
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
   
    export CONDA_PACKAGE=`conda build --output conda_recipe/`
    conda install --yes $CONDA_PACKAGE
    conda update --yes libnetcdf   # KLUDGE to upgrade downgraded libnetcdf
    mkdir foo   # required so source directory not picked up during tests
    cd foo
else
    python setup.py build_ext --inplace
fi
