"""
Reading and RadX Grid objects.

    read_radx_grid

"""

import warnings

import netCDF4

from ..config import get_metadata
from ..core import Grid
from ..io.cfradial import _ncvar_to_dict
from ..io.common import _test_arguments


def read_radx_grid(filename, exclude_fields=None, **kwargs):
    """
    Read a netCDF grid file produced by radx2grid within LROSE.

    Parameters
    ----------
    filename : str
        Filename of netCDF grid file to read.  This file must have been
        produced by :py:func:`write_grid` or have identical layout.

    Other Parameters
    ----------------
    exclude_fields : list
        A list of fields to exclude from the grid object.

    Returns
    -------
    grid : Grid
        Grid object containing gridded data.

    """
    # test for non empty kwargs
    _test_arguments(kwargs)

    if exclude_fields is None:
        exclude_fields = []

    reserved_variables = [
        "time",
        "x",
        "y",
        "z",
        "origin_latitude",
        "origin_longitude",
        "origin_altitude",
        "point_x",
        "point_y",
        "point_z",
        "projection",
        "point_latitude",
        "point_longitude",
        "point_altitude",
        "radar_latitude",
        "radar_longitude",
        "radar_altitude",
        "radar_name",
        "radar_time",
        "base_time",
        "time_offset",
        "ProjectionCoordinateSystem",
    ]

    dset = netCDF4.Dataset(filename, mode="r")

    # metadata
    metadata = {k: getattr(dset, k) for k in dset.ncattrs()}

    # required reserved variables
    time = _ncvar_to_dict(dset.variables["time"])

    # below is altered from original read_grid to read in origin info from radx file
    # Get dict information from lat array in radx file
    try:
        origin_latitude = _ncvar_to_dict(dset.variables["lat0"])
        # set data single value from Grid info origin from radx file
        origin_latitude["data"] = [
            dset.variables["grid_mapping_0"].latitude_of_projection_origin
        ]
        # Get dict information from lon array in radx file
        origin_longitude = _ncvar_to_dict(dset.variables["lon0"])
        # set data single value from Grid info origin from radx file
        try:
            origin_longitude["data"] = [
                dset.variables["grid_mapping_0"].longitude_of_central_meridian
            ]
        except AttributeError:
            origin_longitude["data"] = [
                dset.variables["grid_mapping_0"].longitude_of_projection_origin
            ]
    except KeyError:
        origin_latitude = get_metadata("latitude")
        origin_latitude["data"] = [
            dset.variables["grid_mapping_0"].latitude_of_projection_origin
        ]

        origin_longitude = get_metadata("longitude")
        try:
            origin_longitude["data"] = [
                dset.variables["grid_mapping_0"].longitude_of_central_meridian
            ]
        except AttributeError:
            origin_longitude["data"] = [
                dset.variables["grid_mapping_0"].longitude_of_projection_origin
            ]

    # only need first alt and it needs to be in meters
    origin_altitude = _ncvar_to_dict(dset.variables["z0"])
    origin_altitude["data"] = [origin_altitude["data"][0] * 1000]
    origin_altitude["units"] = "m"

    # Below is altered from rad_grid to read in the correct variables and
    # change the units of the data in the radx file. Get catersian grid
    # spacing information from radx file and set values and units to meters.
    x = _ncvar_to_dict(dset.variables["x0"])
    x["data"] = x["data"] * 1000
    x["units"] = "m"
    y = _ncvar_to_dict(dset.variables["y0"])
    y["data"] = y["data"] * 1000
    y["units"] = "m"
    z = _ncvar_to_dict(dset.variables["z0"])
    z["data"] = z["data"] * 1000
    z["units"] = "m"

    # Projection
    # Below has been significantly altered form the original rad_grid to
    # obtain the correct info from the radx file to be plotted with Cartopy.
    # Create necessary projection params dict for cartopy and pyart from
    # values in grid_mapping_0 variable in radx created file false easting and
    # northing need to be in meters
    projection = get_grid_projection_dict(dset)

    # map _include_lon_0_lat_0 key to bool type
    if "_include_lon_0_lat_0" in projection:
        v = projection["_include_lon_0_lat_0"]
        projection["_include_lon_0_lat_0"] = {"true": True, "false": False}[v]

    # read in the fields
    fields = {}

    # fields in the file has a shape of (1, nz, ny, nx) with the leading 1
    # indicating time but should shaped (nz, ny, nx) in the Grid object
    field_shape = tuple(len(dset.dimensions[d]) for d in ["z0", "y0", "x0"])
    field_shape_with_time = (1,) + field_shape

    # check all non-reserved variables, those with the correct shape
    # are added to the field dictionary, if a wrong sized field is
    # detected a warning is raised
    field_keys = [k for k in dset.variables if k not in reserved_variables]
    for field in field_keys:
        if field in exclude_fields:
            continue
        field_dic = _ncvar_to_dict(dset.variables[field])
        if field_dic["data"].shape == field_shape_with_time:
            field_dic["data"].shape = field_shape
            fields[field] = field_dic
        else:
            warnings.warn(f"Field {field} skipped due to incorrect shape")

    # radar_ variables
    if "radar_latitude" in dset.variables:
        radar_latitude = _ncvar_to_dict(dset.variables["radar_latitude"])
    else:
        radar_latitude = None

    if "radar_longitude" in dset.variables:
        radar_longitude = _ncvar_to_dict(dset.variables["radar_longitude"])
    else:
        radar_longitude = None

    if "radar_altitude" in dset.variables:
        radar_altitude = _ncvar_to_dict(dset.variables["radar_altitude"])
    else:
        radar_altitude = None

    if "radar_name" in dset.variables:
        radar_name = _ncvar_to_dict(dset.variables["radar_name"])
    else:
        radar_name = None

    if "radar_time" in dset.variables:
        radar_time = _ncvar_to_dict(dset.variables["radar_time"])
    else:
        radar_time = None

    dset.close()

    return Grid(
        time,
        fields,
        metadata,
        origin_latitude,
        origin_longitude,
        origin_altitude,
        x,
        y,
        z,
        projection=projection,
        radar_latitude=radar_latitude,
        radar_longitude=radar_longitude,
        radar_altitude=radar_altitude,
        radar_name=radar_name,
        radar_time=radar_time,
    )


def get_grid_projection_dict(dset):
    grid_map_name = dset.variables["grid_mapping_0"].grid_mapping_name
    if grid_map_name == "azimuthal_equidistant":
        projstr = "pyart_aeqd"
        projection = {
            "lat_0": dset.variables["grid_mapping_0"].latitude_of_projection_origin,
            "lon_0": dset.variables["grid_mapping_0"].longitude_of_projection_origin,
            "proj": projstr,
            "units": "m",
            "vunits": "m",
            "x_0": dset.variables["grid_mapping_0"].false_easting * 1000,
            "y_0": dset.variables["grid_mapping_0"].false_northing * 1000,
        }

    elif grid_map_name == "transverse_mercator":
        projstr = "tmerc"
        projection = {
            "+k_0": dset.variables["grid_mapping_0"].scale_factor_at_central_meridian,
            "lat_0": dset.variables["grid_mapping_0"].latitude_of_projection_origin,
            "lon_0": dset.variables["grid_mapping_0"].longitude_of_central_meridian,
            "proj": projstr,
            "units": "m",
            "vunits": "m",
            "x_0": dset.variables["grid_mapping_0"].false_easting * 1000,
            "y_0": dset.variables["grid_mapping_0"].false_northing * 1000,
        }

    return projection
