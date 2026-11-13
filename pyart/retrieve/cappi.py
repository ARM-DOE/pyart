"""
Constant Altitude Plan Position Indicator
"""

import numpy as np
from netCDF4 import num2date
from pandas import to_datetime
from scipy.interpolate import RectBivariateSpline

from pyart.core import Radar


def create_cappi(
    radar,
    fields=None,
    height=2000,
    gatefilter=None,
    vel_field="velocity",
    same_nyquist=True,
    nyquist_vector_idx=0,
    max_height_diff=None,
):
    """
    Create a Constant Altitude Plan Position Indicator (CAPPI) from radar data.

    Parameters
    ----------
    radar : Radar
        Py-ART Radar object containing the radar data.
    fields : list of str, optional
        List of radar fields to be used for creating the CAPPI.
        If None, all available fields will be used. Default is None.
    height : float, optional
        The altitude at which to create the CAPPI. Default is 2000 meters.
    gatefilter : GateFilter, optional
        A GateFilter object to apply masking/filtering to the radar data.
        Default is None.
    vel_field : str, optional
        The name of the velocity field to be used for determining the Nyquist velocity.
        Default is 'velocity'.
    same_nyquist : bool, optional
        Whether to only stack sweeps with the same Nyquist velocity.
        Default is True.
    nyquist_vector_idx : int, optional
        Index for the Nyquist velocity vector if `same_nyquist` is True.
        Default is 0.
    max_height_diff : float, optional
        Maximum allowed difference (meters) between the requested CAPPI height
        and the closest available gate height. Gates farther than this
        threshold are masked out. If None, a data-driven tolerance equal to
        twice the median vertical gate spacing is used. Default is None.

    Returns
    -------
    Radar
        A Py-ART Radar object containing the CAPPI at the specified height.

    Notes
    -----
    CAPPI (Constant Altitude Plan Position Indicator) is a radar visualization
    technique that provides a horizontal view of meteorological data at a fixed altitude.
    Reference: https://glossary.ametsoc.org/wiki/Cappi

    Author
    ------
    Hamid Ali Syed (@syedhamidali)
    """

    if fields is None:
        fields = list(radar.fields.keys())

    # Initialize the first sweep as the reference
    first_sweep = 0

    # Initialize containers for the stacked data and nyquist velocities
    data_stack = []
    nyquist_stack = []
    gate_z_stack = []

    # Process each sweep individually
    for sweep in range(radar.nsweeps):
        sweep_slice = radar.get_slice(sweep)
        try:
            nyquist = radar.get_nyquist_vel(sweep=sweep)
            nyquist = np.round(nyquist)
        except LookupError:
            print(
                f"Nyquist velocity unavailable for sweep {sweep}. Estimating using maximum velocity."
            )
            nyquist = radar.fields[vel_field]["data"][sweep_slice].max()

        # Extract and sort azimuth angles (used across all fields)
        azimuth = radar.azimuth["data"][sweep_slice]
        azimuth_sorted_idx = np.argsort(azimuth)
        azimuth_sorted = azimuth[azimuth_sorted_idx]
        time = radar.time["data"][sweep_slice][azimuth_sorted_idx]

        # Gate heights from radar geometry already include the 4/3 Earth-radius refractivity model
        sweep_gate_z = radar.gate_z["data"][sweep_slice][azimuth_sorted_idx]

        if sweep == first_sweep:
            azimuth_final = azimuth_sorted
            time_final = time
        else:
            # Interpolate gate heights to the reference azimuth grid (keeps refractivity-aware heights aligned)
            gate_interpolator = RectBivariateSpline(
                azimuth_sorted, radar.range["data"], sweep_gate_z
            )
            sweep_gate_z = gate_interpolator(azimuth_final, radar.range["data"])

        sweep_data = {}

        for field in fields:
            data = radar.get_field(sweep, field)

            # Apply gatefilter if provided
            if gatefilter is not None:
                data = np.ma.masked_array(
                    data, gatefilter.gate_excluded[sweep_slice, :]
                )

            data = data[azimuth_sorted_idx]

            if sweep != first_sweep:
                # Interpolate data for consistent azimuth ordering across sweeps
                interpolator = RectBivariateSpline(
                    azimuth_sorted, radar.range["data"], data
                )
                data = interpolator(azimuth_final, radar.range["data"])

            sweep_data[field] = data[np.newaxis, :, :]

        data_stack.append(sweep_data)
        nyquist_stack.append(nyquist)
        gate_z_stack.append(sweep_gate_z[np.newaxis, :, :])

    nyquist_stack = np.array(nyquist_stack)

    # Filter for sweeps with similar Nyquist velocities
    if same_nyquist:
        nyquist_range = nyquist_stack[nyquist_vector_idx]
        nyquist_mask = np.abs(nyquist_stack - nyquist_range) <= 1
        data_stack = [
            sweep_data for i, sweep_data in enumerate(data_stack) if nyquist_mask[i]
        ]
        gate_z_stack = [
            sweep_gate_z
            for i, sweep_gate_z in enumerate(gate_z_stack)
            if nyquist_mask[i]
        ]

    z_3d = np.concatenate(gate_z_stack, axis=0)

    # Find nearest gates to requested height using refractivity/curvature-adjusted geometry
    height_idx = np.argmin(np.abs(z_3d - height), axis=0)
    selected_gate_z = np.take_along_axis(
        z_3d, height_idx[np.newaxis, :, :], axis=0
    ).squeeze(axis=0)

    # Derive a height tolerance if none provided. Use the median positive gate-to-gate
    # spacing across the volume to approximate vertical resolution.
    if max_height_diff is None:
        z_sorted = np.sort(z_3d.ravel())
        dz = np.diff(z_sorted)
        dz = dz[dz > 0]
        if dz.size == 0:
            derived_tol = 1000.0  # fallback to 1 km
        else:
            derived_tol = 2.0 * np.median(dz)
        height_tol = max(500.0, derived_tol)  # enforce a sensible floor
    else:
        height_tol = float(max_height_diff)

    # Mask gates that are too far from the requested height (e.g., asking for 20 km
    # when the volume tops out near 14 km). Also warn when the request is outside
    # the available volume.
    max_available = z_3d.max()
    min_available = z_3d.min()
    if height > max_available + height_tol or height < min_available - height_tol:
        print(
            f"Warning: requested CAPPI height {height} m is outside available gates "
            f"({min_available:.0f}-{max_available:.0f} m). Output will be fully masked."
        )

    height_mask = np.abs(selected_gate_z - height) > height_tol
    # Generate CAPPI for each field using data_stack
    fields_data = {}
    for field in fields:
        data_3d = np.concatenate(
            [sweep_data[field] for sweep_data in data_stack], axis=0
        )

        # Sort azimuth for all sweeps
        dim0 = data_3d.shape[1:]

        # Extract the data slice corresponding to the requested height using gate geometry (curved Earth)
        CAPPI = np.take_along_axis(
            data_3d, height_idx[np.newaxis, :, :], axis=0
        ).squeeze(axis=0)
        CAPPI = np.ma.array(CAPPI, mask=height_mask)

        # Retrieve units and handle case where units might be missing
        units = radar.fields[field].get("units", "").lower()

        # Determine valid_min and valid_max based on units
        if units == "dbz":
            valid_min, valid_max = -10, 80
        elif units in ["m/s", "meters per second"]:
            valid_min, valid_max = -100, 100
        elif units == "db":
            valid_min, valid_max = -7.9, 7.9
        else:
            # If units are not found or don't match known types, set default values or skip masking
            valid_min, valid_max = None, None

        # If valid_min or valid_max are still None, set them to conservative defaults or skip
        if valid_min is None:
            print(f"Warning: valid_min not set for {field}, using default of -1e10")
            valid_min = -1e10  # Conservative default
        if valid_max is None:
            print(f"Warning: valid_max not set for {field}, using default of 1e10")
            valid_max = 1e10  # Conservative default

        # Apply valid_min and valid_max masking
        if valid_min is not None:
            CAPPI = np.ma.masked_less(CAPPI, valid_min)
        if valid_max is not None:
            CAPPI = np.ma.masked_greater(CAPPI, valid_max)

        # Convert to masked array with the specified fill value
        CAPPI.set_fill_value(radar.fields[field].get("_FillValue", np.nan))
        CAPPI = np.ma.masked_invalid(CAPPI)
        CAPPI = np.ma.masked_outside(CAPPI, valid_min, valid_max)

        fields_data[field] = {
            "data": CAPPI,
            "units": radar.fields[field]["units"],
            "long_name": f"CAPPI {field} at {height} meters",
            "comment": f"CAPPI {field} calculated at a height of {height} meters",
            "_FillValue": radar.fields[field].get("_FillValue", np.nan),
        }

    # Set the elevation to zeros for CAPPI
    elevation_final = np.zeros(dim0[0], dtype="float32")

    # Since we are using the whole volume scan, report mean time
    try:
        dtime = to_datetime(
            num2date(radar.time["data"], radar.time["units"]).astype(str),
            format="ISO8601",
        )
    except ValueError:
        dtime = to_datetime(
            num2date(radar.time["data"], radar.time["units"]).astype(str)
        )
    dtime = dtime.mean()

    time = radar.time.copy()
    time["data"] = time_final
    time["mean"] = dtime

    # CAPPI represents a single synthesized sweep on the reference azimuth grid
    n_rays = azimuth_final.size
    sweep_number = {"data": np.array([0], dtype="int32")}
    sweep_mode = {"data": np.array(["ppi"], dtype="object")}
    sweep_start_ray_index = {"data": np.array([0], dtype="int32")}
    sweep_end_ray_index = {"data": np.array([n_rays - 1], dtype="int32")}
    fixed_angle = {"data": np.array([0.0], dtype="float32")}

    # Create the Radar object with the new CAPPI data.  Keep site metadata but
    # reset sweep descriptors and azimuth to match the CAPPI grid so display
    # libraries do not mis-align fields (which showed up as a 90° rotation).
    return Radar(
        time=time,
        _range=radar.range.copy(),
        fields=fields_data,
        metadata=radar.metadata.copy(),
        scan_type="ppi",
        latitude=radar.latitude.copy(),
        longitude=radar.longitude.copy(),
        altitude=radar.altitude.copy(),
        sweep_number=sweep_number,
        sweep_mode=sweep_mode,
        fixed_angle=fixed_angle,
        sweep_start_ray_index=sweep_start_ray_index,
        sweep_end_ray_index=sweep_end_ray_index,
        azimuth={"data": azimuth_final},
        elevation={"data": elevation_final},
        instrument_parameters=radar.instrument_parameters,
    )
