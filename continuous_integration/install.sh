#!/bin/bash

set -e
# use next line to debug this script
#set -x

# Install Miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b
export PATH=/home/travis/miniconda3/bin:$PATH
conda config --set always_yes yes
conda config --set show_channel_urls true
conda update -q conda

## Create a testenv with the correct Python version
conda env create -f continuous_integration/environment-$PYTHON_VERSION.yml
source activate testenv
# TODO install cylp and  glpk in Python 2.7

# upgrade pip
pip install --upgrade pip

# install coverage modules
conda install -c conda-forge -q pytest-cov
if [[ "$COVERALLS" == "true" ]]; then
    conda install -c conda-forge -q coveralls
fi

# install Py-ART
export RSL_PATH=~/miniconda3/envs/testenv

if [[ "$FROM_RECIPE" == "true" ]]; then
    source deactivate
    conda install -q conda-build
    conda install -q jinja2 setuptools
    conda config --add channels conda-forge
    conda config --add channels jjhelmus
    conda build --no-test --python $PYTHON_VERSION -q conda_recipe/
   
    export CONDA_PACKAGE=`conda build --python $PYTHON_VERSION --output conda_recipe/ | grep bz2`
    source activate testenv
    conda install -q $CONDA_PACKAGE
    mkdir foo   # required so source directory not picked up during tests
    cd foo
else
    pip install -e .
fi
