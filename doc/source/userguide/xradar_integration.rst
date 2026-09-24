.. _xradar_integration:

=================================
Reading Data: xradar and Py-ART
=================================

This page answers the question new users ask first: **which IO path do I
use, and what works on the object it gives me?**

Reading data in Py-ART today
=============================

The recommended path for reading radar data is **xradar-first**: open the
file with `xradar <https://docs.openradarscience.org/projects/xradar/en/stable/>`_
into an xarray ``DataTree``, then attach Py-ART's ``pyart`` accessor to get
a :class:`pyart.xradar.Xradar` object that behaves like a
:class:`pyart.core.Radar` for the purposes of running Py-ART algorithms,
correction routines, and displays.

.. code-block:: python

    import xradar as xd
    import pyart

    # Access sample cfradial1 data from Py-ART and read using xradar
    filename = pyart.testing.get_test_data("swx_20120520_0641.nc")
    tree = xd.io.open_cfradial1_datatree(filename)

    # Attach the Py-ART accessor -- ``radar`` now exposes Radar-like
    # attributes/methods (fields, azimuth, elevation, get_elevation, ...)
    radar = tree.pyart.to_radar()

xradar ships ``open_*_datatree`` readers for CfRadial1, CfRadial2/FM301,
ODIM_H5, Furuno, Iris/Sigmet, NEXRAD Level 2, and others -- see the
`xradar documentation <https://docs.openradarscience.org/projects/xradar/en/stable/>`_
for the full list. Any algorithm in Py-ART that accepts a ``radar``
argument can be handed the resulting :class:`pyart.xradar.Xradar` object
directly, in place of a :class:`pyart.core.Radar`. Code that needs to
accept either interchangeably (for example, a helper function) can
normalize the input with :func:`pyart.xradar.to_pyart_radar`, which
returns a :class:`~pyart.core.Radar` or :class:`~pyart.xradar.Xradar`
unchanged and wraps a bare ``xarray.DataTree`` for you:

.. code-block:: python

    import pyart

    radar = pyart.xradar.to_pyart_radar(tree)  # Radar, Xradar, or DataTree in

When to use ``pyart.io.read`` instead
--------------------------------------

Py-ART's built-in readers (:func:`pyart.io.read` and friends, e.g.
:func:`pyart.io.read_cfradial`, :func:`pyart.io.read_nexrad_archive`) remain
available and return a classic :class:`pyart.core.Radar`. Reach for them
when:

- You need a format xradar does not (yet) read, or a Py-ART-specific reader
  such as MDV, UF, CHL, or RSL-backed files (see :mod:`pyart.io` and
  :mod:`pyart.aux_io`).
- You depend on behavior specific to :class:`pyart.core.Radar` that is not
  (yet) mirrored on :class:`pyart.xradar.Xradar` -- check the
  :ref:`supported-API table <xradar-supported-api>` below first.
- You are maintaining existing code and are not ready to migrate.

.. code-block:: python

    import pyart

    radar = pyart.io.read(filename)  # returns a pyart.core.Radar

Reader migration: the roadmap direction
=========================================

New radar-format IO support is expected to land in `xradar
<https://docs.openradarscience.org/projects/xradar/en/stable/>`_ going
forward rather than as new readers in :mod:`pyart.io`. Py-ART's own legacy
readers will be deprecated over time as xradar's format coverage grows,
matching "Moderate #1" of the rc3.0 roadmap
(``roadmaps/pyart-roadmap-rc3.0.pdf``). No legacy reader is being removed
today -- this is a direction to plan for, not an immediate break.

.. _xradar-supported-api:

Supported-API table
=====================

:class:`pyart.xradar.Xradar` wraps an xradar ``DataTree`` and duck-types
:class:`pyart.core.Radar`: most attributes and methods used by Py-ART's
correction, retrieval, mapping, and graphing subpackages work unchanged.
The table below lists what has been verified to work, sourced from
``pyart/xradar/accessor.py`` and the test suite in ``tests/xradar/``
(notably ``tests/xradar/test_compat_matrix.py``, which runs the same
algorithm against a :class:`~pyart.core.Radar` and the equivalent
:class:`~pyart.xradar.Xradar` and asserts numerically equal results).

``Xradar`` (radar-level)
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Member
     - Notes
   * - ``fields``
     - ``numpy.ma.MaskedArray`` per field, masked where the source
       variable is NaN (or was already masked)
   * - ``azimuth``, ``elevation``, ``range``, ``time``, ``fixed_angle``
     - Standard Radar-style ``dict(data=..., **attrs)``
   * - ``latitude``, ``longitude``, ``altitude``
     - Populated from the DataTree root
   * - ``sweep_number``, ``sweep_mode``, ``sweep_start_ray_index``,
       ``sweep_end_ray_index``
     - Derived from the combined per-sweep data
   * - ``rays_per_sweep``
     - Lazily derived from the sweep start/end indices
   * - ``gate_x``/``y``/``z``, ``gate_longitude``/``latitude``,
       ``gate_altitude``
     - Computed the same way as on ``Radar``
   * - ``get_elevation(sweep)``
     - Per-sweep elevation slice
   * - ``get_azimuth(sweep)``
     - Per-sweep azimuth slice
   * - ``get_gate_area(sweep)``
     - Per-gate area, same formula as ``Radar``
   * - ``get_gate_x_y_z``, ``get_gate_lat_lon_alt``
     - Per-sweep gate locations, computed via ``xradar.georeference()``
       under the hood
   * - ``get_nyquist_vel``
     - From ``instrument_parameters`` when present
   * - ``info(level=...)``
     - Same three levels as ``Radar.info``
   * - ``add_field``, ``add_field_like``, ``add_filter``
     - Mutates both ``self.fields`` and the underlying ``DataTree``
       sweep datasets
   * - ``extract_sweeps``
     - Returns a new ``Xradar``
   * - ``instrument_parameters``
     - Populated from a ``radar_parameters`` subgroup and/or
       root-level variables when present, else ``{}``
   * - ``radar_calibration``
     - Populated from a ``radar_calibration`` subgroup when present in
       the DataTree, else ``None``
   * - ``antenna_transition``, ``rays_are_indexed``, ``ray_angle_res``,
       ``target_scan_rate``, ``scan_rate``, ``altitude_agl``,
       ``rotation``, ``tilt``, ``roll``, ``drift``, ``heading``,
       ``pitch``, ``georefs_applied``
     - ``None`` when the source file has no corresponding
       variable/group -- **never fabricated**. Populated when present.
   * - ``scan_type``
     - ``"ppi"`` by default, or pass ``scan_type="rhi"`` to
       ``tree.pyart.to_radar(scan_type="rhi")`` -- RHI sweeps are
       combined on elevation rather than azimuth

``Xgrid`` (grid-level)
------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Member
     - Notes
   * - ``fields``, ``x``/``y``/``z``, ``origin_latitude`` etc.
     - Same shape and layout as ``pyart.core.Grid``
   * - ``get_point_longitude_latitude``, ``add_field``, ``to_xarray``
     - Same as ``Grid``
   * - ``projection_proj``
     - ``pyproj.Proj`` built from the ``projection`` attrs; raises
       ``ValueError`` for the ``pyart_aeqd`` projection (no
       ``pyproj.Proj`` equivalent) and
       :class:`~pyart.exceptions.MissingOptionalDependency` if
       ``pyproj`` is not installed
   * - ``write(filename)``
     - Writes the grid via :func:`pyart.io.grid_io.write_grid`

Known caveats
==============

- **Fields are masked, not filled.** ``Xradar.fields[name]["data"]`` is a
  ``numpy.ma.MaskedArray``; any NaN in the source variable (xarray's
  decoded ``_FillValue``) is masked, mirroring how :class:`pyart.core.Radar`
  masks missing gates.
- **RHI sweeps are supported** by passing ``scan_type="rhi"`` to
  ``tree.pyart.to_radar()``; internally rays are combined per-sweep on
  elevation (which varies per ray in an RHI) instead of azimuth (which is
  ~constant across an RHI sweep).
- **``radar_calibration``** is populated only when the source DataTree has a
  ``radar_calibration`` child group (mirroring how ``instrument_parameters``
  is sourced from a ``radar_parameters`` subgroup); otherwise it is ``None``,
  matching :class:`~pyart.core.Radar`'s handling of files that lack this
  information.
- **Comparing a file read both ways is not bit-for-bit.** ``pyart.io.read``
  and xradar's ``open_*_datatree`` readers can disagree on ray bookkeeping
  for messy real-world files (e.g. which ray straddling the 0/360 azimuth
  wrap belongs to which sweep); this is a property of two independent
  parsers, not an ``Xradar`` accessor defect. See
  ``tests/xradar/test_compat_matrix.py`` for how the compatibility test
  suite works around this to get an apples-to-apples numerical comparison.

Worked example
================

.. code-block:: python

    import xradar as xd
    import pyart

    filename = pyart.testing.get_test_data("swx_20120520_0641.nc")
    tree = xd.io.open_cfradial1_datatree(filename)
    radar = tree.pyart.to_radar()

    # Grid using 11 vertical levels, and 101 horizontal grid cells at a
    # resolution of 1 km -- pyart.map works directly on the Xradar object
    grid = pyart.map.grid_from_radars(
        (radar,),
        grid_shape=(11, 101, 101),
        grid_limits=(
            (0.0, 10_000),
            (-50_000.0, 50_000.0),
            (-50_000, 50_000.0),
        ),
    )

More end-to-end examples are in the example gallery:

- :ref:`sphx_glr_examples_xradar_plot_xradar.py`
- :ref:`sphx_glr_examples_xradar_plot_grid_xradar.py`
- :ref:`sphx_glr_examples_xradar_plot_dealias_xradar.py`
- :ref:`sphx_glr_examples_xradar_plot_rhi_xradar.py`

See also :ref:`pyart_2_0` for the rest of the Py-ART 2.0 API changes, and
:ref:`overview` for how the ``Radar``/xradar data models compare.
