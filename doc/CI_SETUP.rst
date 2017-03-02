Py-ART Continuous Integration Setup
===================================

Py-ART makes use of continuous integration (CI) to run the library's unit tests
against every pull request and change made to the repository. This document
gives a brief explanation of the details of this setup and hints on how to fix
the CI when it breaks.


Services
--------

Py-ART's CI setup is hosted on two services, `Travis CI`_ and `AppVeyor`_.
Both of these are hosted continuous integration services that are free for open
source projects. Travis CI runs the tests on a Linux platform and AppVeyor
runs on Windows.  Currently no tests are run on macOS although Travis CI does
support this platform.

.. _Travis CI : https://travis-ci.org/
.. _AppVeyor : https://www.appveyor.com/


Travis CI
---------

The CI setup on Travis CI is controlled by the ``.travis.yml`` file in the
repository root.  This file contains multiple sections which define the steps
taken in a CI run.  The `Travis CI documentation`_ provide details on these
sections and Travis CI in general.

Briefly the sections in the file include:

- An encrypted environment variable used to push the documentation built
  upon each merge to the master branch to the `pyart-docs-travis`_ repo.
- A matrix section which defines the various combinations of Python versions
  and other parameters which will be tested.  Currently Python 2.7, 3.4, 3.5
  and 3.6 are tested with a second Python 2.7 run also testing that the
  package can be built from a conda recipe and that the docs build.
- An install section which executes the ``install.sh`` script in the
  continuous_integration directory.
- A script section which runs the tests using ``nosetests``
- A after_success section which reports test coverage to `Coveralls`_.
  Currently nothing is done with this information but it may be useful later
  to view how test coverage has changed over time.
- The after_success section also runs the ``build_docs.sh`` script in the
  continuous_integration directory for the appropriate matrix element.  On
  non-merge commits this tests that the documentation can be build.  On merge
  commits this builds the docs and uploads them to the `pyart-docs-travis`_
  repository which makes them available on the `Py-ART documentation page`_.

The bulk of the test environment setup is done by the ``install.sh`` script in
the continuous_integration directory.  This script downloads, installs and
configures the latest version of Miniconda3, creates a conda environment for
the Python version being tested with the various requirements specified in the
environment-XXX.yml file, and builds Py-ART in that environment.  For the
matrix element where Py-ART is build from a conda recipe, the if/then block
installs conda-build, builds Py-ART from the recipe in the repository and
installs this conda package into the test conda environment.

.. _Travis CI documentation : https://docs.travis-ci.com/
.. _pyart-docs-travis : https://github.com/ARM-DOE/pyart-docs-travis
.. _Coveralls : https://coveralls.io/
.. _Py-ART documentation page : http://arm-doe.github.io/pyart-docs-travis/


Fixing issues with Travis CI
----------------------------

Most of the issues that occur with Travis CI which require updates to the
setup are caused by conflicts and issues with the conda packages used to
create the test environment. These packages are specified in the
enviroment-X.X.yml files in the continuous_integration directory. Over time
these may need to be updated.

To debug Travis CI issue use the logs provided for the builds themselves to see
where the issues are occurring.  Un-commenting the ``set -x`` line in the
``install.sh`` script may be helpful when debugging Travis CI issue.  With this
line uncommented the commands being executed are echoed in the log.

Additionally, it can be helpful to run the CI builds locally to quickly test
out changes to the setup and try out specific commands.  Travis CI provides
the docker `images`_ used by their service and instructions on how to
`run these locally`_.

.. _images : https://quay.io/organization/travisci
.. _run these locally : https://docs.travis-ci.com/user/common-build-problems/#Running-a-Container-Based-Docker-Image-Locally


AppVeyor
--------

The second CI service used by Py-ART is `AppVeyor_`.  This setup is controlled
by the ``appveyor.yml`` file in the repository and builds and runs the unit
tests on a Windows platform.  This services is slower than Travis CI so only
two Python versions are tested (2.7 and 3.5).  All setup and configuration is
done in the ``appveyor.yml`` file besides a setup function used by many
projects which can be found in the ``run_with_env.cmd`` file in the
continuous_integration/appveyor directory.  This file contains comments which
explain the steps being performed and some debugging details.  Refer to the
`AppVeyor documentation`_ for details.

.. _AppVeyor : https://www.appveyor.com/
.. _AppVeyor documentation : https://www.appveyor.com/docs/
