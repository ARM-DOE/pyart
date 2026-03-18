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
We now offer the option to use xradar for IO, with the following interface (a typical gridding workflow is shown below):

.. code-block:: python

    import xradar as xd
    import pyart

    # Access sample cfradial1 data from Py-ART and read using xradar
    filename = get_test_data("swx_20120520_0641.nc")
    tree = xd.io.open_cfradial1_datatree(filename)

    # Add the associated pyart methods - ensuring compatibility with Py-ART functionality
    radar = tree.pyart.to_radar()

    # Grid using 11 vertical levels, and 101 horizontal grid cells at a resolution on 1 km
    grid = pyart.map.grid_from_radars(
        (radar,),
        grid_shape=(11, 101, 101),
        grid_limits=(
            (0.0, 10_000),
            (-50_000.0, 50_000.0),
            (-50_000, 50_000.0),
        ),
    )

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
