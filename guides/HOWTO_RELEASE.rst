How to make a Py-ART release
============================

This file gives an overview of the steps needed to make a release of Py-ART.
including making binary and source files available on PyPI and conda-forge.


Make the tagged git release
---------------------------

Change the ``ISRELEASE`` line in the **setup.py** file in the root directory to 
True. Check this change into git and create a tagged commit with the version
string. The commit message should be "REL: Release X.X.X" with the proper
version filled in. This commit and the tag should be pushed directly to the
ARM-DOE repository on GitHub. This can be done using the following, assuming
upstream points to https://github.com/ARM-DOE/pyart.git.

::

    git add setup.py
    git commit -m "REL: Release 1.8.0"
    git tag -a v1.8.0 -m "Py-ART version 1.8.0"
    git push --tags upstream master

As an example see the `commit`_ for the 1.8.0 release.

.. _commit : https://github.com/ARM-DOE/pyart/commit/d1d08c1d2136051d865d7a6269326328035db7c0


Create and upload a source distribution to PyPI
-----------------------------------------------

To create a release and the associated source distribution on PyPI, begin by 
cloning the ARM-DOE repository into a temporary location and checking out the
tagged release. Starting from a fresh clone is recommended to prevent any
extra files from accidentally being included in the source distribution. 

With the tagged version checked out, generate the source distribution file
using ``python setup.py sdist --formats=gztar``. This will create a
arm_pyart-X.X.X.tar.gz file in the dist directory. Tests should be run to check
that this file can be used to install a working version of Py-ART.

Once tested the tar.gz source distribution file can be uploaded to PyPI using
`twine <https://pypi.python.org/pypi/twine>`_, installing if needed, with the 
command ``twine upload dist/*``. This command will prompt the user for their
PyPI credentials.  After uploading check that this file is available for
download from 
`pypi.python.org/pypi/arm_pyart <https://pypi.python.org/pypi/arm_pyart>`_.

Some older releases also included a zip file, these are no longer included as
PyPI only allows a single source distribution to be uploaded.

A typical workflow for creating and uploading the source distribution to PyPI
is as follows:

::

    git clone git@github.com:ARM-DOE/pyart.git
    cd pyart
    git checkout tags/v1.8.0
    # test the .tar.gz file in the dist directory
    twine upload dist/*


Prepare the release notes
-------------------------

Prepare a short description of the changes, improvements and deprecations in 
the new release. The format used in previous releases makes a good template
for this text. The following git commands can be used to list all merge commits
and the contributors to a particular version which may be helpful when drafting
these notes.

::

    git log --merges v1.7.0..v1.8.0  # list all merge commits
    git shortlog v1.7.0..v1.8.0      # log of all commit by author


Update the GitHub Release
-------------------------

From the ARM-DOE Github page, click on the Releases tab and make a new release.
Use the existing tag pushed for the version. Include the release notes. Upload
the tar.gz file from PyPI in the release. Some older releases also included a
zip file, these are no longer included as PyPI only allows a single source
distribution to be uploaded.


Create and upload version specific documentation
------------------------------------------------

The Py-ART documentation page contains links to documentation and examples from
the latest release. These pages must be updated when a new release is made.
To build the documentation for a release, go to the doc directory and run the
``rebuild_full_docs.sh`` shell script. Note that this should be done while the
release commit is checked out and on a system with all optional requirements
for the examples and sphinx install. Additionally, the files needed for the
examples must be copied into the correct locations in the examples directory.

Once the documentation has been built, the contents of the doc/build/html
directory needs to be committed to the gh-pages branch of the ARM-DOE repository
in the dev directory. Typically this process is best done in a directory where
the ARM-DOE repository has been cloned and the gh-pages branch checked out.

First, the current dev directory should be moved to a versioned directory using
``git mv dev v1_7`` and a commit made. Then create an empty ``dev`` directory 
and copy all the files from the doc/build/html directory. Check these new 
files into git and create a commit. Finally, push these changes to the ARM-DOE
repository. After a few minutes check that the new documentation is available
at `http://arm-doe.github.io/pyart/dev/ <http://arm-doe.github.io/pyart/dev/>`_.


Prepare the repository for the next version
-------------------------------------------

To begin the development of the next version of Py-ART the ``MINOR`` variable
should be incremented and ``ISRELEASED`` should be changed back to False.  
Make this change to the ``setup.py`` file in the root directory and create a 
commit with the message ``MAINT: start v1.9+dev``. See the commit for the 
`start of the 1.8+dev`_ as an example. Push this commit directly to the
ARM-DOE repository.

.. _start of the 1.8+dev : https://github.com/ARM-DOE/pyart/commit/e3d01232ac3787da38456152da6b9db3ab46d4d7

If any features are to be depreciated in this version, this is a good time
to remove them.


Create conda packages on conda-forge
------------------------------------

Conda package for Py-ART are created using the conda-forge feedstock. To
create a new release first fork the `arm_pyart-feedstock`_ repository or update
your fork. Then clone the fork and checkout a new branch. Update the first
few lines of the ``recipe/meta.yaml`` with the new version number and a sha256
hash of the file on PyPI. Check this change in as a new commit. Next rerender
the feedstock using ``conda smithy rerender``, installing ``conda-smithy``
if needed.  

Finally, push the branch up to your fork and create a pull request against the
conda-forge feedstock. The CI services will test the recipe by building conda
packages on all three platforms and Python+NumPy combinations. Once these
pass the pull request can be merged and the CI services will again build the
packages and uploaded them to the conda-forge Anaconda cloud channel.

The `arm_pyart-feedstock`_ repository has additional details on how this
process works.

.. _arm_pyart-feedstock : https://github.com/conda-forge/arm_pyart-feedstock


Announce the release on the pyart-user mailing list
---------------------------------------------------

Announce the new release on the `Py-ART mailing list`_. Use a previous
announcement as a template and include the release notes written above.

.. _Py-ART mailing list : http://groups.google.com/group/pyart-users/
