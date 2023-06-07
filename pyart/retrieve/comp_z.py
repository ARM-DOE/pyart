"""
Calculate the composite reflectivity

"""

import copy

import numpy as np
from netCDF4 import num2date
from pandas import to_datetime
from scipy.interpolate import interp2d

from pyart.core import Radar


def composite_reflectivity(radar, field="reflectivity", gatefilter=None):
    """
    Composite Reflectivity

    Often a legacy product, composite reflectivity is:
    "A display or mapping of the maximum radar reflectivity factor at any
    altitude as a function of position on the ground." - AMS Glossary
    This is more useful for the dry regions of the world, where maximum
    reflectivity values are found aloft, as opposed to the lowest scan.
    Alternatively this is useful for comparing to NWP since composite Z
    is a standard output of NWP.

    Why this is not as easy as one would think: Turns out the data are
    not natively stored with index 0 being azimuth 0. Likely due to the
    physical spinning of the radar antenna.

    Author: Randy J. Chase (@dopplerchase)

    Parameters
    ----------
    radar : Radar
        Radar object used.
    field : str
        Reflectivity field name to use to look up reflectivity data. In the
        radar object. Default field name is 'reflectivity'.
    gatefilter : GateFilter
        GateFilter instance. None will result in no gatefilter mask being
        applied to data.

    Returns
    -------
    radar : Radar
        The radar object containing the radar dimensions, metadata and
        composite field.

    References
    ----------
    American Meteorological Society, 2022: "Composite reflectivity". Glossary of Meteorology,
    http://glossary.ametsoc.org/wiki/Composite_reflectivity

    """

    # Determine the lowest sweep (used for metadata and such)
    minimum_sweep = np.min(radar.sweep_number["data"])

    # loop over all measured sweeps
    for sweep in sorted(radar.sweep_number["data"]):
        # get start and stop index numbers
        sweep_slice = radar.get_slice(sweep)

        # grab radar data
        z = radar.get_field(sweep, field)
        z_dtype = z.dtype

        # Use gatefilter
        if gatefilter is not None:
            mask_sweep = gatefilter.gate_excluded[sweep_slice, :]
            z = np.ma.masked_array(z, mask_sweep)

        # extract lat lons
        lon = radar.gate_longitude["data"][sweep_slice, :]
        lat = radar.gate_latitude["data"][sweep_slice, :]

        # get the range and time
        ranges = radar.range["data"]
        time = radar.time["data"]

        # get azimuth
        az = radar.azimuth["data"][sweep_slice]
        # get order of azimuths
        az_ids = np.argsort(az)

        # reorder azs so they are in order
        az = az[az_ids]
        z = z[az_ids]
        lon = lon[az_ids]
        lat = lat[az_ids]
        time = time[az_ids]

        # if the first sweep, store re-ordered lons/lats
        if sweep == minimum_sweep:
            azimuth_final = az
            time_final = time
            lon_0 = copy.deepcopy(lon)
            lon_0[-1, :] = lon_0[0, :]
            lat_0 = copy.deepcopy(lat)
            lat_0[-1, :] = lat_0[0, :]

        else:
            # Configure the intperpolator
            z_interpolator = interp2d(ranges, az, z, kind="linear")

            # Apply the interpolation
            z = z_interpolator(ranges, azimuth_final)

        # if first sweep, create new dim, otherwise concat them up
        if sweep == minimum_sweep:
            z_stack = copy.deepcopy(z[np.newaxis, :, :])
        else:
            z_stack = np.concatenate([z_stack, z[np.newaxis, :, :]])

    # now that the stack is made, take max across vertical
    compz = z_stack.max(axis=0).astype(z_dtype)

    # since we are using the whole volume scan, report mean time
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

    # return dict, because this is was pyart does with lots of things
    fields = {}
    fields["composite_reflectivity"] = {
        "data": compz,
        "units": "dBZ",
        "long_name": "composite_reflectivity",
        "comment": "composite reflectivity computed from calculating the max radar value in each radar gate vertically after reordering",
    }

    time = radar.time.copy()
    time["data"] = time_final
    time["mean"] = dtime

    gate_longitude = radar.gate_longitude.copy()
    gate_longitude["data"] = lon_0
    gate_longitude["comment"] = "reordered longitude grid, [az,range]"
    gate_latitude = radar.gate_latitude.copy()
    gate_latitude["data"] = lat_0
    gate_latitude["comment"] = "reordered latitude grid, [az,range]"

    _range = radar.range.copy()
    metadata = radar.metadata.copy()
    scan_type = radar.scan_type
    latitude = radar.latitude.copy()
    longitude = radar.longitude.copy()
    altitude = radar.altitude.copy()
    instrument_parameters = radar.instrument_parameters

    sweep_number = radar.sweep_number.copy()
    sweep_number["data"] = np.array([0], dtype="int32")
    sweep_mode = radar.sweep_mode.copy()
    sweep_mode["data"] = np.array([radar.sweep_mode["data"][0]])
    ray_shape = compz.shape[0]
    azimuth = radar.azimuth.copy()
    azimuth["data"] = azimuth_final
    elevation = radar.elevation.copy()
    elevation["data"] = np.zeros(ray_shape, dtype="float32")
    fixed_angle = radar.fixed_angle.copy()
    fixed_angle["data"] = np.array([0.0], dtype="float32")
    sweep_start_ray_index = radar.sweep_start_ray_index.copy()
    sweep_start_ray_index["data"] = np.array([0], dtype="int32")
    sweep_end_ray_index = radar.sweep_end_ray_index.copy()
    sweep_end_ray_index["data"] = np.array([ray_shape - 1], dtype="int32")

    return Radar(
        time,
        _range,
        fields,
        metadata,
        scan_type,
        latitude,
        longitude,
        altitude,
        sweep_number,
        sweep_mode,
        fixed_angle,
        sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth,
        elevation,
        instrument_parameters=instrument_parameters,
    )
