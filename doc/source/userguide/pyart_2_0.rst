.. _pyart_2_0:

==========
Py-ART 2.0
==========

In preparation for version 2.0.0 of Py-ART, codes were standardized for consistency purposes as further defined in the `Contributor's Guide <https://arm-doe.github.io/pyart/userguide/CONTRIBUTING.html>`_.  These changes will break some users code as the API has changed.  This guide will detail the changes for each module.

How to Try Py-ART 2.0
=====================

The Py-ART 2.0 release candidate can be installed directly from github - this is still a work in progress, feedback is welcome!::

    pip install git+https://github.com/ARM-DOE/pyart@release/2.0

Input/Output (IO)
=================
We now offer the option to use xradar for IO. See :ref:`xradar_integration`
for the full xradar-first reading guide, the supported-API table for
:class:`pyart.xradar.Xradar`/:class:`pyart.xradar.Xgrid`, and a worked
gridding example.

Correct
=======
The `dealias_fourdd <https://arm-doe.github.io/pyart/API/generated/pyart.correct.dealias_fourdd.html>`_ algorithm has been removed given the now unsupported RSL library.

It is recommended that users move to the `region-based dealiasing algorithm <https://arm-doe.github.io/pyart/API/generated/pyart.correct.dealias_region_based.html>`_.

Graph
=====
Colormaps have been moved to a dedicated package outside Py-ART, `cmweather <https://cmweather.readthedocs.io/>`_.

For example, visualizing our grid mentioned previously, it is recommended to install/import cmweather and change the colormap name from `pyart_ChaseSpectral` to `ChaseSpectral`

.. code-block:: python

    import cmweather

    display = pyart.graph.GridMapDisplay(grid)
    display.plot_grid(
        "reflectivity_horizontal", level=0, vmin=-20, vmax=60, cmap="ChaseSpectral"
    )
