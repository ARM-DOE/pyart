How to make a Py-ART release
============================

This file gives an overview of the steps needed to make a release of Py-ART.
including making binary and source files available on PyPI and conda-forge.


Make the tagged git release
---------------------------

First you will create a new tag and message associated with that tag.
After that you will push the new tag to Py-ART's upstream. Upstream should
point to https://github.com/ARM-DOE/pyart.git. Make sure when creating a new tag
to not commit anything to the tag until the release process is done. Otherwise
the automatic version creator will create a .post extension that will break PyPI.

::

    git tag -a v1.8.0 -m "Py-ART version 1.8.0"
    git push --tags upstream master

Prepare the release notes
-------------------------

There are two ways to prepare release notes for the new release. Manually and
automatically.

If manually, prepare a short description of the changes, improvements and deprecations in
the new release. The format used in previous releases makes a good template
for this text. The following git commands can be used to list all merge commits
and the contributors to a particular version which may be helpful when drafting
these notes.

::

    git log --merges v1.7.0..v1.8.0  # list all merge commits
    git shortlog v1.7.0..v1.8.0      # log of all commit by author

Automatically, when you click draft a new release found here
https://github.com/ARM-DOE/pyart/releases/new there is a button to click
called Generate release notes. Which will generate release notes between
the new tag and the one prior.

Update the GitHub Release and Trigger the PyPI Upload
-----------------------------------------------------

From the ARM-DOE Github page, click on the Releases tab and make a new release.
https://github.com/ARM-DOE/pyart/releases Use the existing tag pushed for the
version. Include the release notes. Once you select release, a Github action
will be triggered, generating the source files and uploding them to PyPI.


Documentation
-------------
Py-ART's documentation will update on its own via successful build doc
pull requests.

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

Note, the conda-forge bot, will sometimes automatically create a PR on the
arm_pyart-feedstock when a new version is uploaded to PyPI. Thus avoiding
the conda-forge process above. Be aware, bot will not change pinning of run and
host dependencies for you. Recommended following the method above if the recipe
changes.


Announce the release on the Open Radar Forum
--------------------------------------------

Announce the new release on the `Open Radar Forum`_. Use a previous
announcement as a template and include the release notes written above.

.. _Open Radar Forum : https://openradar.discourse.group/
