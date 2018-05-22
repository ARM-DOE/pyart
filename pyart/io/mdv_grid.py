"""
pyart.io.mdv_grid
==================

Utilities for reading and writing of MDV grid files.

.. autosummary::
    :toctree: generated/

    write_grid_mdv
    read_grid_mdv
    _time_dic_to_datetime

"""


import datetime
import warnings

from netCDF4 import num2date, date2num
import numpy as np

from ..config import FileMetadata, get_fillvalue, get_metadata
from ..core.grid import Grid
from .common import make_time_unit_str, _test_arguments, prepare_for_read
from ..lazydict import LazyLoadDict
from . import mdv_common


def write_grid_mdv(filename, grid, mdv_field_names=None,
                   field_write_order=None):
    """
    Write grid object to MDV file.

    Create a MDV file containing data from the provided grid instance.

    The MDV file will contain parameters from the 'source' key if contained
    in grid.metadata.  If this key or parameters related to the radar location
    and name are not present in the grid a default or sentinel value.
    will be written in the MDV file in the place of the parameter.

    Grid fields will be saved in float32 unless the `_Write_as_dtype` key is
    present.

    Parameters
    ----------
    filename : str or file-like object.
        Filename of MDV file to create.  If a file-like object is specified
        data will be written using the write method.
    grid : Grid
        Grid object from which to create MDV file.
    mdv_field_names : dict or None, optional
        Mapping between grid fields and MDV data type names. Field names
        mapped to None or with no mapping will be excluded from
        writing.  If None, the same field names will be used.
    field_write_order : list or None, optional
        Order in which grid fields should be written out in the MDV file.
        None, the default, will determine a valid order automatically.

    Notes
    -----
    Do to limitations of the MDV format, not all grid objects are writable.
    To write a grid the following conditions must be satisfied:

        * XY grid must be regular (equal spacing), Z can be irregular.
        * The number of Z levels must not exceed 122.
        * Fields can be encoded in the file using the '_Write_as_dtype' key
          specifying one of 'uint8', 'uint16' or 'float32'.  Use the
          'scale_factor' and 'add_offset' keys to specify scaling.  Field
          data in the Grid object should be uncompressed, that is to say
          it has had the scaling applied.

    """
    # first of all firm field list
    if field_write_order is None:
        field_write_order = list(grid.fields.keys())

    # remove fields not in mdv_field_names
    if mdv_field_names is not None:
        for ifield, field in enumerate(field_write_order):
            if field not in mdv_field_names or mdv_field_names[field] is None:
                field_write_order.pop(ifield)

    # mount empty mdv file
    grid_shape = grid.fields[field_write_order[0]]['data'].shape
    nz, ny, nx = grid_shape
    nfields = len(field_write_order)
    if nz > mdv_common.MDV_MAX_VLEVELS:
        warnings.warn(('%i vlevels exceed MDV_MAX_VLEVELS = %i. Extra ' +
                       'levels will be ignored') %
                      (nz, mdv_common.MDV_MAX_VLEVELS))
        nz = mdv_common.MDV_MAX_VLEVELS
    mdv = mdv_common.MdvFile(None)
    mdv.field_headers = mdv._get_field_headers(nfields)
    mdv.vlevel_headers = mdv._get_vlevel_headers(nfields)
    mdv.fields_data = [None] * nfields
    mdv.projection = 'flat'

    # fill headers
    d = mdv.master_header
    grid_datetime = _time_dic_to_datetime(grid.time)
    mdv.times['time_centroid'] = grid_datetime
    mdv.times["time_begin"] = grid_datetime
    mdv.times["time_end"] = grid_datetime
    mdv._time_dict_into_header()

    d["data_dimension"] = 3  # XXX are grid's always 3d?
    # =DATA_SYNTHESIS, I don't realy know, so miscellaneous!
    d["data_collection_type"] = 3
    if (grid.z['units'] == 'm' or grid.z['units'] == 'meters'):
        d["native_vlevel_type"] = 4
    elif (grid.z['units'] == '\xc2' or grid.z['units'] == 'degree'):
        d["native_vlevel_type"] = 9
    d["vlevel_type"] = d["native_vlevel_type"]
    d["nfields"] = nfields
    d["max_nx"] = nx
    d["max_ny"] = ny
    d["max_nz"] = nz
    td = datetime.datetime.utcnow() - datetime.datetime(1970, 1, 1, 0, 0)
    d["time_written"] = int(round(td.microseconds + (td.seconds + td.days *
                                  24 * 3600) * 10**6) / 10**6)

    # sensor location using radar_ attribute or origin_ if not available
    if grid.radar_longitude is not None and grid.nradar != 0:
        d["sensor_lon"] = grid.radar_longitude['data'][0]
    else:
        d["sensor_lon"] = grid.origin_longitude['data'][0]

    if grid.radar_latitude is not None and grid.nradar != 0:
        d["sensor_lat"] = grid.radar_latitude['data'][0]
    else:
        d["sensor_lat"] = grid.origin_latitude['data'][0]

    if grid.radar_altitude is not None and grid.nradar != 0:
        d["sensor_alt"] = grid.radar_altitude['data'][0] / 1000.
    else:
        d["sensor_alt"] = grid.origin_altitude['data'][0] / 1000.

    for meta_key, mdv_key in mdv_common.MDV_METADATA_MAP.items():
        if meta_key in grid.metadata:
            d[mdv_key] = grid.metadata[meta_key].encode("ASCII")

    for ifield, field in enumerate(field_write_order):
        d = mdv.field_headers[ifield]
        l = mdv.vlevel_headers[ifield]

        # fill fields_header
        d["nx"] = nx
        d["ny"] = ny
        d["nz"] = nz
        d["proj_type"] = mdv_common.PROJ_FLAT
        try:
            dtype = np.dtype(grid.fields[field]['_Write_as_dtype'])
        except KeyError:    # default to float32 encoding
            dtype = np.float32
        if dtype == np.uint8:
            d["encoding_type"] = mdv_common.ENCODING_INT8
            d["data_element_nbytes"] = 1
        elif dtype == np.uint16:
            d["encoding_type"] = mdv_common.ENCODING_INT16
            d["data_element_nbytes"] = 2
        elif dtype == np.float32 or dtype == np.float64:
            d["encoding_type"] = mdv_common.ENCODING_FLOAT32
            d["data_element_nbytes"] = 4
        else:
            raise TypeError("Unsuported encoding %s, encoding must be "
                            "uint8, uint16 or float32 as specfied by"
                            "the '_Write_as_dtype key" % dtype)
        d["compression_type"] = 3   # zlib

        d["scaling_type"] = 4  # SCALING_SPECIFIED (by the user)
        d["native_vlevel_type"] = mdv.master_header["vlevel_type"]
        d["vlevel_type"] = mdv.master_header["vlevel_type"]
        d["data_dimension"] = 3
        d["proj_origin_lat"] = grid.origin_latitude['data'][0]
        d["proj_origin_lon"] = grid.origin_longitude['data'][0]
        d["grid_dx"] = (grid.x['data'][1] - grid.x['data'][0]) / 1000.
        d["grid_dy"] = (grid.y['data'][1] - grid.y['data'][0]) / 1000.
        d["grid_minx"] = (grid.x['data'][0]) / 1000.
        d["grid_miny"] = (grid.y['data'][0]) / 1000.
        if "scale_factor" in grid.fields[field].keys():
            d["scale"] = grid.fields[field]["scale_factor"]
        if "add_offset" in grid.fields[field].keys():
            d["bias"] = grid.fields[field]["add_offset"]

        # bad_data prioritise _FillValue
        if "_FillValue" in grid.fields[field].keys():
            d["bad_data_value"] = grid.fields[field]["_FillValue"]
        elif "missing_value" in grid.fields[field].keys():
            d["bad_data_value"] = grid.fields[field]["missing_value"]
        else:
            d["bad_data_value"] = get_fillvalue()
        # missing_data prioritise missing_value
        if "missing_value" in grid.fields[field].keys():
            d["missing_data_value"] = grid.fields[field]["missing_value"]
        elif "_FillValue" in grid.fields[field].keys():
            d["missing_data_value"] = grid.fields[field]["_FillValue"]
        else:
            d["missing_data_value"] = get_fillvalue()

        d["min_value"] = np.amax(grid.fields[field]['data'])
        d["max_value"] = np.amin(grid.fields[field]['data'])
        if "standard_name" in grid.fields[field].keys():
            d["field_name_long"] = (
                grid.fields[field]["standard_name"].encode("ASCII"))
        elif "long_name" in grid.fields[field].keys():
            d["field_name_long"] = (
                grid.fields[field]["long_name"].encode("ASCII"))
        if mdv_field_names is not None:
            d["field_name"] = mdv_field_names[field].encode("ASCII")
        else:
            d["field_name"] = field.encode("ASCII")
        if "units" in grid.fields[field].keys():
            d["units"] = grid.fields[field]["units"].encode("ASCII")
        d["transform"] = "none".encode("ASCII")  # XXX not implemented

        # fill vlevels_header
        typ = [0] * 122
        level = [0] * 122
        for iz in range(nz):
            typ[iz] = d["vlevel_type"]
            level[iz] = grid.z["data"][iz] / 1000.
        l["type"] = typ
        l["level"] = level

        # put data to field
        mdv.fields_data[ifield] = grid.fields[field]["data"]

    # write the file
    mdv.write(filename)


def read_grid_mdv(filename, field_names=None, additional_metadata=None,
                  file_field_names=False, exclude_fields=None,
                  delay_field_loading=False, **kwargs):
    """
    Read a MDV file to a Grid Object.

    Parameters
    ----------
    filename : str
        Name of MDV file to read or file-like object pointing to the
        beginning of such a file.
    field_names : dict, optional
        Dictionary mapping MDV data type names to radar field names. If a
        data type found in the file does not appear in this dictionary or has
        a value of None it will not be placed in the radar.fields dictionary.
        A value of None, the default, will use the mapping defined in the
        Py-ART configuration file.
    additional_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve metadata from during this read.
        This metadata is not used during any successive file reads unless
        explicitly included.  A value of None, the default, will not
        introduct any addition metadata and the file specific or default
        metadata as specified by the Py-ART configuration file will be used.
    file_field_names : bool, optional
        True to use the MDV data type names for the field names. If this
        case the field_names parameter is ignored. The field dictionary will
        likely only have a 'data' key, unless the fields are defined in
        `additional_metadata`.
    exclude_fields : list or None, optional
        List of fields to exclude from the grid object. This is applied
        after the `file_field_names` and `field_names` parameters.
    delay_field_loading : bool
        True to delay loading of field data from the file until the 'data'
        key in a particular field dictionary is accessed.  In this case
        the field attribute of the returned Radar object will contain
        LazyLoadDict objects not dict objects.

    Returns
    -------
    grid : Grid
        Grid object containing data from MDV file.

    Notes
    -----
    This function can only read cartesian MDV files with fields
    compressed with gzip or zlib. For polar files see
    :py:func:`pyart.io.read_mdv`

    MDV files and Grid object are not fully interchangeable.  Specific
    limitation include:

        * All fields must have the same shape and dimensions.
        * All fields must have the same projection.
        * Vlevels types must not vary.
        * Projection must not be PROJ_POLAR_RADAR (9) or PROJ_RHI_RADAR (13).
        * Correct unit in the Z axis are just availible for 'vlevel_type'
          equal to VERT_TYPE_Z(4), VERT_TYPE_ELEV(9), VERT_TYPE_AZ(17),
          VERT_TYPE_PRESSURE(3) and VERT_TYPE_THETA(7).
        * The behavior in cases of 2D data is unknown but most likely will not
          fail.

    """
    # test for non empty kwargs
    _test_arguments(kwargs)

    # XXX add tests for conversion limitations
    # create metadata retrieval object
    filemetadata = FileMetadata('mdv', field_names, additional_metadata,
                                file_field_names, exclude_fields)
    mdv = mdv_common.MdvFile(prepare_for_read(filename))

    # time dictionaries
    units = make_time_unit_str(mdv.times['time_begin'])
    time = get_metadata('grid_time')
    time['data'] = np.array([date2num(mdv.times['time_centroid'], units)])
    time['units'] = units

    # origin dictionaries
    origin_altitude = get_metadata('origin_altitude')
    origin_altitude['data'] = np.array(
        [mdv.master_header["sensor_alt"] * 1000.], dtype='float64')

    origin_latitude = get_metadata('origin_latitude')
    origin_latitude['data'] = np.array(
        [mdv.master_header["sensor_lat"]], dtype='float64')

    origin_longitude = get_metadata('origin_longitude')
    origin_longitude['data'] = np.array(
        [mdv.master_header["sensor_lon"]], dtype='float64')

    # grid coordinate dictionaries
    nz = mdv.master_header["max_nz"]
    ny = mdv.master_header["max_ny"]
    nx = mdv.master_header["max_nx"]
    z_line = mdv.vlevel_headers[0]["level"][0:nz]
    y_start = mdv.field_headers[0]["grid_miny"] * 1000.
    x_start = mdv.field_headers[0]["grid_minx"] * 1000.
    y_step = mdv.field_headers[0]["grid_dy"] * 1000.
    x_step = mdv.field_headers[0]["grid_dx"] * 1000.

    if mdv.field_headers[0]["proj_type"] == mdv_common.PROJ_LATLON:
        # x and y dimensions have units of degrees
        xunits = 'degree_E'
        yunits = 'degree_N'
        y_start = mdv.field_headers[0]["grid_miny"]
        x_start = mdv.field_headers[0]["grid_minx"]
        y_step = mdv.field_headers[0]["grid_dy"]
        x_step = mdv.field_headers[0]["grid_dx"]
    elif (mdv.field_headers[0]["proj_type"] != mdv_common.PROJ_POLAR_RADAR and
          mdv.field_headers[0]["proj_type"] != mdv_common.PROJ_RHI_RADAR):
        xunits = 'm'
        yunits = 'm'

    if mdv.field_headers[0]["vlevel_type"] == 4:  # VERT_TYPE_Z
        zunits = 'm'
        z_line = [e * 1000. for e in z_line]
    elif mdv.field_headers[0]["vlevel_type"] == 9:  # VERT_TYPE_ELEV
        zunits = 'degree'  # elevation
    elif mdv.field_headers[0]["vlevel_type"] == 17:  # VERT_TYPE_AZ
        zunits = 'degree'  # azimuth
    elif mdv.field_headers[0]["vlevel_type"] == 3:  # VERT_TYPE_PRESSURE
        zunits = 'mb'
    elif mdv.field_headers[0]["vlevel_type"] == 7:  # VERT_TYPE_THETA
        zunits = 'kelvin'
    else:
        warnings.warn(("While reading MDV found unexpected 'vlevel_type'" +
                      " (%i), units in the z axis set to 'unknown'") %
                      mdv.field_headers[0]["vlevel_type"])
        zunits = 'unknown'

    x = get_metadata('x')
    x['data'] = np.linspace(x_start, x_start + x_step * (nx-1), nx)
    x['units'] = xunits

    y = get_metadata('y')
    y['data'] = np.linspace(y_start, y_start + y_step * (ny-1), ny)
    y['units'] = yunits

    z = get_metadata('z')
    z['data'] = np.array(z_line, dtype='float64')
    z['units'] = zunits

    # metadata
    metadata = filemetadata('metadata')
    for meta_key, mdv_key in mdv_common.MDV_METADATA_MAP.items():
        metadata[meta_key] = mdv.master_header[mdv_key]

    # fields
    fields = {}
    mdv_fields = mdv._make_fields_list()
    for mdv_field in set(mdv_fields):
        field_name = filemetadata.get_field_name(mdv_field)
        if field_name is None:
            continue
        field_dic = filemetadata(field_name)

        field_dic['_FillValue'] = get_fillvalue()
        dataextractor = mdv_common._MdvVolumeDataExtractor(
            mdv, mdv.fields.index(mdv_field), get_fillvalue(), two_dims=False)
        if delay_field_loading:
            field_dic = LazyLoadDict(field_dic)
            field_dic.set_lazy('data', dataextractor)
        else:
            field_dic['data'] = dataextractor()
        fields[field_name] = field_dic

    if not delay_field_loading:
        mdv.close()
    return Grid(time, fields, metadata,
                origin_latitude, origin_longitude, origin_altitude, x, y, z)


# This function may be helpful in other cases and could be moved in common
def _time_dic_to_datetime(dic):
    """ Return a datetime for the first element in a time dictionary. """
    if 'calendar' in dic:
        calendar = dic['calendar']
    else:
        calendar = 'standard'
    return num2date(dic['data'][0], dic['units'], calendar)
