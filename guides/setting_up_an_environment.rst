Setting up an Environment
=========================


Anaconda
--------

Creating environments using Anaconda is recommended due to the ability to
create more than one environment. It is also recommended because you can
keep dependencies separate from one another that might conflict if you had
them all in your root environment. For example, if you had all the dependencies
for a Pandas environment and all the dependencies for a Py-ART environment in
your root environment, there might be conflicts between channels and packages.
So Anaconda allows you to create multiple environments to avoid these issues.

To download and install `Anaconda <https://www.anaconda.com/download/#>`_.

While Anaconda is downloading, it will ask if you want to set a path to it, or
let Anaconda set a default path. After choosing, Anaconda should finish
downloading. After it is done, exit the terminal and open a new one to make
sure the environment path is set. If conda command is not found, there is help
on running conda and fixing the environment path, found here:

* `How to Run Conda <https://stackoverflow.com/questions/18675907/how-to-run-conda>`_

Setting a Channel
-----------------

Anaconda has a cloud that stores many of its packages. It is recommended, at
times, to use the conda-forge channel instead. Conda-Forge is a community led
collection of packages, and typically contains the most recent versions of the
packages required for Py-ART. Also Py-ART is on Conda-Forge. Having packages in
an environment, within the same channel, helps avoid conflict issues. To add
conda-forge as the priority channel, simply do::

        conda config --add channels conda-forge

You can also just flag the channel when conda install packages such as::

        conda install -c conda-forge numpy

More on managing channels can be found here:

* `Managing Channels <https://conda.io/docs/user-guide/tasks/manage-channels.html>`_

Creating an Environment
-----------------------

There are a few ways to create a conda environment for using Py-ART or other
packages. One way is to use the environment file, found here:

* https://github.com/ARM-DOE/pyart/blob/master/environment.yml

To create an environment using this file, use the command::

        conda env create -f environment.yml

This will then create an environment called pyart_env that can be activated
by::

        source activate pyart_env

or deactivated after use::

        source deactivate pyart_env

Once the environment is created and activated, you can install more packages
into the environment by simply conda installing them. An example of this is,
if you want Jupyter Notebook to run in that enviroment with those packages::

        conda install -c conda-forge jupyter notebook

while that environment is activated. Another way to create a conda environment
is by doing it from scratch using the conda create command. An example of this::

        conda create -n pyart_env -c conda-forge python=3.8 arm_pyart netCDF4
        cartopy scipy numpy matplotlib

This will also create an environment called pyart_env that can be activate the
same way, as mentioned above. To then run your coding editor within the
environment, run in the command line::

        python

or::

        ipython

or::

        jupyter notebook

or even::

        spyder

depending on what you installed in your environment and want to use for coding.

Adding Optional Dependencies with setting Paths
-----------------------------------------------

There are other optional dependencies that can enhance the use of Py-ART. One,
such package is `CyLP <https://github.com/jjhelmus/CyLP>`_. To get CyLP to work,
installing of the package `coincbc <https://projects.coin-or.org/Cbc>`_ is
needed as a dependency for CyLP. Simply do::

        conda install -c conda-forge coincbc

within your pyart_env. After that though, the coincbc path needs to be exported
so CyLP knows where to find it during its install. To do this::

        export COIN_INSTALL_DIR=/Users/yourusername/youranacondadir/envs/pyart_env

or real example on a Linux machine::

        export COIN_INSTALL_DIR=/home/zsherman/anaconda3/envs/pyart_env

CyLP was actually adapted by Jonathan Helmus to be Python 3 compatible, so we
will install a specific CyLP branch after doing the export path step above.
GitHub repositories can actually be pip installed within your environment. So
to install the CyLP version we want::

        pip install git+https://github.com/jjhelmus/CyLP.git@py3

This will install a Python 3 compatible version of CyLP found on GitHub.

More Information
----------------

For more an conda and help with conda:

* https://conda.io/docs/
* https://gitter.im/conda/conda
