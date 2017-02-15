Contributor's Guide
===================

Note: This guide is still a work in progress.
Note: When an example shows filename as such::
	
	nosetests filename

filename is the filename and location, such as::

	nosetests /home/zsherman/pyart/pyart/retrieve/tests/test_vad.py

Py-ART Introduction and Information
===================================

The Python ARM Radar Toolkit (Py-ART)
-------------------------------------

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
---------------

- Official source code repository: https://github.com/ARM-DOE/pyart
- HTML documentation: http://arm-doe.github.io/pyart-docs-travis/
- Examples: http://arm-doe.github.io/pyart/dev/auto_examples/index.html
- Mailing List: http://groups.google.com/group/pyart-users/
- Issue Tracker: https://github.com/ARM-DOE/pyart/issues


Citing
------

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
-------

The easiest method for installing Py-ART is to use the conda packages from
the latest release.  To do this you must download and install 
`Anaconda <http://continuum.io/downloads>`_ or 
`Miniconda <http://continuum.io/downloads>`_.  
Then use the following command in a terminal or command prompt to install
the latest version of Py-ART::

    conda install -c conda-forge arm_pyart

To update an older version of Py-ART to the latest release use::

    conda update -c conda-forge arm_pyart

Resources
---------

Pyart:

- https://github.com/openradar/AMS-Short-Course-on-Open-Source-Radar-Software
- https://github.com/EVS-ATMOS/pyart_short_course
- https://www.youtube.com/watch?v=diiP-Q3bKZw
- http://arm-doe.github.io/pyart/dev/auto_examples/index.html

Git:

- https://git-scm.com/book/en/v2
- https://docs.scipy.org/doc/numpy/dev/gitwash/development_workflow.html



Proper Code Style
=================

Py-ART follows pep8 coding standards. To make sure your code follows the pep8
style, you can use a variety of modules that will check for you. Two popular
pep8 check modules are pycodestyle and pylint.

For more on pep8 style:

- https://www.python.org/dev/peps/pep-0008/

To install pycodestyle::

        pip install pycodestyle

To use pycodestyle::

        pycodestyle filename

pylint is also useful to check for pep8 compliance and should already be
installed with python 3.5.

To get a detailed pylint report::

        pylint filename

If you want to just see what line number and the issue, just use::

        pylint -r n filename


Python File Setup
=================

In a new .py file, the top of the code should have the function's location,
sphinx comments for template configuration, and public and private functions
within the .py file. Public functions are listed first and then private
functions with and underscore in front of the function. A space is needed
between the last function and the closing docstring quotation.

For example:

.. code-block:: python

        """
	pyart.retrieve.velocity_azimuth_display
	=======================================
	
	Retrieval of VADs from a radar object.

	.. autosummary::
    	    :toctreeL generated/
    	    :template: dev_template.rst

	    velocity_azimuth_display
            _inverse_dist_squared
            _Average1D

        """

Following the introduction code, modules are then added. Main imports come
first, followed by 'from imports'.

For example:

.. code-block:: python

        import numpy as np
        import numpy.ma as ma
        from scipy.interpolate import interp1d

        from ..core import HorizontalWindProfile

Following the main function def line, but before the code within it, a doc
string is needed to explain parameters, returns, references if needed, and
other helpful information.

For example:

.. code-block:: python

        def velocity_azimuth_display(
            radar, velocity=None, z_want=None, valid_ray_min=16,
            gatefilter=False, window=2):
	
            """
   	    Velocity azimuth display.

            Note: If a specific sweep is desired, before using the
            velocity_azimuth_display function, use:
            radar = radar.extract_sweeps([0])

            Parameters
            ----------
            radar : Radar
                Radar object used.
            velocity : string
                Velocity field to use for VAD calculation.
                If None, the default velocity field will be used.

            Other Parameters
            ----------------
            z_want : array
                Height array user would like for the VAD
                calculation. None will result in a z_want of
        	np.linspace and use of _inverse_dist_squared
        	and _Average1D functions. Note, height must have
        	same shape as expected u_wind and v_wind if user
        	provides z_want.
    	    valid_ray_min : int
        	Amount of rays required to include that level in
        	the VAD calculation.
            gatefilter : GateFilter
        	Used to correct the velocity field before its use
        	in the VAD calculation. Uses Py-ART's region dealiaser.
    	    window : int
        	Value to use for window calculation in _Averag1D
        	function.

            Returns
            -------
    	    height : array
        	Heights in meters above sea level at which horizontal winds were
        	sampled.
    	    speed : array
        	Horizontal wind speed in meters per second at each height.
    	    direction : array
        	Horizontal wind direction in degrees at each height.
    	    u_wind : array
        	U-wind mean in meters per second.
    	    v_wind : array
        	V-wind mean in meters per second.

    	    Reference
    	    ----------
    	    K. A. Browning and R. Wexler, 1968: The Determination
    	    of Kinematic Properties of a Wind Field Using Doppler
	    Radar. J. Appl. Meteor., 7, 105–113

    	    """
            
As seen, each variable has what type of object it is, an explaination of what
it is, mention of units, and if a variable has a default value, statement of
what that default value is and why.

When adding a new function to pyart it is important to add your function to
the _init.py file under the corresponding pyart folder.

Create a test for your function and have assert from numpy test the known
values to the calculated values. If changes are made in the future to pyart,
nose will use the test created to see if the function is still valid and
produces the same values. 

An example:

.. code-block:: python

        def test_vad():
            test_radar = pyart.testing.make_target_radar()
            height = np.arange(0, 1000, 200)
            speed = np.ones_like(height) * 5
            direction = np.array([0, 90, 180, 270, 45])
            profile = pyart.core.HorizontalWindProfile(height, speed, direction)
            sim_vel = pyart.util.simulated_vel_from_profile(test_radar, profile)
            
            test_radar.add_field('velocity', sim_vel,
                                 replace_existing=True)

            velocity = 'velocity'
            z_start = 0
            z_end = 10
            z_count = 5

            vad_height = ([0., 2.5, 5., 7.5, 10.])
            vad_speed = ([4.98665725, 4.94020686, 4.88107152,
                          4.81939374, 4.75851962])
            vad_direction = ([359.84659496, 359.30240553, 358.58658589,
                              357.81073051, 357.01353486])
            u_wind = ([0.01335138, 0.06014712, 0.12039762,
                       0.18410404, 0.24791911])
            v_wind = ([-4.98663937, -4.9398407, -4.87958641,
                       -4.81587601, -4.75205693])

            vad = pyart.retrieve.velocity_azimuth_display(test_radar,
                                                          velocity,
                                                          z_start, z_end,
                                                          z_count)

            assert_almost_equal(vad.height, vad_height, 8)
            assert_almost_equal(vad.speed, vad_speed, 8)
            assert_almost_equal(vad.direction, vad_direction, 8)
            assert_almost_equal(vad.u_wind, u_wind, 8)
            assert_almost_equal(vad.v_wind, v_wind, 8)

To install nose::

   		conda install nose

To run all tests in pyart with nose::

		nosetests --exe pyart

All test with in depth details::

		nosetests -v -s

Just one file::

		nosetests filename


GitHub
======

When contributing to pyart, the changes created should be in a new branch
under your forked repository. Let’s say your adding a new map display.
Instead of creating that new function in your master branch. Create a new
branch called ‘cartopy_map’. If everything checks out and the admin
accepts the pull request, you can then merge the master branch and
cartopy_map branch. 

To delete a branch both locally and remotely, if done with it::

		git push origin --delete <branch_name>
		git branch -d <branch_name>

or in this case::
		
		git push origin --delete cartopy_map
		git branch -d cartopy_map


To create a new branch, the command is `git checkout -b <branch_name>`.
If you type `git status` it will inform you of the branch you are in.

To switch between branches, simply type::

		git checkout <branch_name>

When commiting to GitHub, start the statement with a acronym such as
‘ADD:’ depending on what your commiting, could be ‘MAINT:’ or
‘BUG:’ or more. Then following should be a short statement such as
“ADD: Adding cartopy map display.”, but after the short statement, before
finishing the quotations, hit enter and in your terminal you can then type
a more in depth description on what your commiting. 

If you would like to type your commit in the terminal and skip the default
editor::

	git commit -m "PEP: Removing whitespace from vad.py."

To use the default editor(in Linux, usually VIM), simply type::

	git commit

One thing to keep in mind is before doing a pull request, update your
branches with the original upstream repository.

This could be done by::

	git fetch upstream

After creating a pull request through GitHub, two outside code checkers,
Appveyor and TravisCI will determine if the code past all checks. If the code
fails either tests, as the pull request sits, make changes to fix the code
and when pushed to GitHub, the pull request will automatically update and
TravisCI and Appveyor will automatically rerun.


GitLab
======


                
