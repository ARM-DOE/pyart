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

        sweep_data = {}

        for field in fields:
            data = radar.get_field(sweep, field)

            # Apply gatefilter if provided
            if gatefilter is not None:
                data = np.ma.masked_array(
                    data, gatefilter.gate_excluded[sweep_slice, :]
                )
            time = radar.time["data"][sweep_slice]

            # Extract and sort azimuth angles
            azimuth = radar.azimuth["data"][sweep_slice]
            azimuth_sorted_idx = np.argsort(azimuth)
            azimuth = azimuth[azimuth_sorted_idx]
            data = data[azimuth_sorted_idx]

            # Store initial lat/lon for reordering
            if sweep == first_sweep:
                azimuth_final = azimuth
                time_final = time
            else:
                # Interpolate data for consistent azimuth ordering across sweeps
                interpolator = RectBivariateSpline(azimuth, radar.range["data"], data)
                data = interpolator(azimuth_final, radar.range["data"])

            sweep_data[field] = data[np.newaxis, :, :]

        data_stack.append(sweep_data)
        nyquist_stack.append(nyquist)

    nyquist_stack = np.array(nyquist_stack)

    # Filter for sweeps with similar Nyquist velocities
    if same_nyquist:
        nyquist_range = nyquist_stack[nyquist_vector_idx]
        nyquist_mask = np.abs(nyquist_stack - nyquist_range) <= 1
        data_stack = [
            sweep_data for i, sweep_data in enumerate(data_stack) if nyquist_mask[i]
        ]

    # Generate CAPPI for each field using data_stack
    fields_data = {}
    for field in fields:
        data_3d = np.concatenate(
            [sweep_data[field] for sweep_data in data_stack], axis=0
        )

        # Sort azimuth for all sweeps
        dim0 = data_3d.shape[1:]
        azimuths = np.linspace(0, 359, dim0[0])
        elevation_angles = radar.fixed_angle["data"][: data_3d.shape[0]]
        ranges = radar.range["data"]

        theta = (450 - azimuths) % 360
        THETA, PHI, R = np.meshgrid(theta, elevation_angles, ranges)
        Z = R * np.sin(PHI * np.pi / 180)

        # Extract the data slice corresponding to the requested height
        height_idx = np.argmin(np.abs(Z - height), axis=0)
        CAPPI = np.array(
            [
                data_3d[height_idx[j, i], j, i]
                for j in range(dim0[0])
                for i in range(dim0[1])
            ]
        ).reshape(dim0)

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

    # Create the Radar object with the new CAPPI data
    return Radar(
        time=radar.time.copy(),
        _range=radar.range.copy(),
        fields=fields_data,
        metadata=radar.metadata.copy(),
        scan_type=radar.scan_type,
        latitude=radar.latitude.copy(),
        longitude=radar.longitude.copy(),
        altitude=radar.altitude.copy(),
        sweep_number=radar.sweep_number.copy(),
        sweep_mode=radar.sweep_mode.copy(),
        fixed_angle=radar.fixed_angle.copy(),
        sweep_start_ray_index=radar.sweep_start_ray_index.copy(),
        sweep_end_ray_index=radar.sweep_end_ray_index.copy(),
        azimuth=radar.azimuth.copy(),
        elevation={"data": elevation_final},
        instrument_parameters=radar.instrument_parameters,
    )
