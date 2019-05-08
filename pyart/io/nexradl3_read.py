"""
pyart.io.nexradl3_read
======================

Functions for reading NEXRAD Level 3 products.

.. autosummary::
    :toctree: generated/

    read_nexrad_level3

"""

import numpy as np

from ..config import FileMetadata, get_fillvalue
from ..core.radar import Radar
from .common import make_time_unit_str, _test_arguments, prepare_for_read
from .nexrad_level3 import NEXRADLevel3File


def read_nexrad_level3(filename, field_names=None, additional_metadata=None,
                       file_field_names=False, exclude_fields=None,
                       include_fields=None, **kwargs):
    """
    Read a NEXRAD Level 3 product.

    Parameters
    ----------
    filename : str
        Filename of NEXRAD Level 3 product file. The files hosted by
        at the NOAA National Climate Data Center [1]_ as well as on the
        NWS WSR-88D Level III Data Collection and Distribution Network
        have been tests. Other NEXRAD Level 3 files may or may not work.
        A file-like object pointing to the beginning of such a file is also
        supported.
    field_names : dict, optional
        Dictionary mapping NEXRAD level 3 product number to radar field names.
        If the product number of the file does not appear in this dictionary
        or has a value of None it will not be placed in the radar.fields
        dictionary. A value of None, the default, will use the mapping
        defined in the metadata configuration file.
    additional_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve metadata from during this read.
        This metadata is not used during any successive file reads unless
        explicitly included. A value of None, the default, will not
        introduct any addition metadata and the file specific or default
        metadata as specified by the metadata configuration file will be used.
    file_field_names : bool, optional
        True to use the product number for the field name. In this case the
        field_names parameter is ignored. The field dictionary will likely
        only have a 'data' key, unless the fields are defined in
        `additional_metadata`.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters. Set
        to None to include all fields specified by include_fields.
    include_fields : list or None, optional
        List of fields to include from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters. Set
        to None to include all fields not specified by exclude_fields.

    Returns
    -------
    radar : Radar
        Radar object containing all moments and sweeps/cuts in the volume.
        Gates not collected are masked in the field data.

    References
    ----------
    .. [1] http://www.ncdc.noaa.gov/
    .. [2] http://www.roc.noaa.gov/wsr88d/Level_III/Level3Info.asp

    """
    # test for non empty kwargs
    _test_arguments(kwargs)

    # create metadata retrieval object
    filemetadata = FileMetadata('nexrad_level3', field_names,
                                additional_metadata, file_field_names,
                                exclude_fields, include_fields)

    # open the file
    nfile = NEXRADLevel3File(prepare_for_read(filename))
    nradials = nfile.packet_header['nradials']
    msg_code = nfile.msg_header['code']

    # time
    time = filemetadata('time')
    time_start = nfile.get_volume_start_datetime()
    time['units'] = make_time_unit_str(time_start)
    time['data'] = np.zeros((nradials, ), dtype='float64')

    # range
    _range = filemetadata('range')
    _range['data'] = nfile.get_range()
    _range['meters_to_center_of_first_gate'] = _range['data'][0]
    _range['meters_between_gates'] = _range['data'][1] - _range['data'][0]

    # fields
    fields = {}
    field_name = filemetadata.get_field_name(msg_code)
    if field_name is None:
        fields = {}
    else:
        dic = filemetadata(field_name)
        dic['_FillValue'] = get_fillvalue()
        dic['data'] = nfile.get_data()
        fields = {field_name: dic}

    # metadata
    metadata = filemetadata('metadata')
    metadata['original_container'] = 'NEXRAD Level 3'

    # scan_type
    scan_type = 'ppi'

    # latitude, longitude, altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')

    lat, lon, height = nfile.get_location()
    latitude['data'] = np.array([lat], dtype='float64')
    longitude['data'] = np.array([lon], dtype='float64')
    altitude['data'] = np.array([height], dtype='float64')

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index
    # sweep_end_ray_index
    sweep_number = filemetadata('sweep_number')
    sweep_mode = filemetadata('sweep_mode')
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')

    sweep_number['data'] = np.array([0], dtype='int32')
    sweep_mode['data'] = np.array(1 * ['azimuth_surveillance'], dtype='S')

    sweep_start_ray_index['data'] = np.array([0], dtype='int32')
    sweep_end_ray_index['data'] = np.array([nradials - 1], dtype='int32')

    # azimuth, elevation, fixed_angle
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')
    fixed_angle = filemetadata('fixed_angle')
    azimuth['data'] = nfile.get_azimuth()
    elev = nfile.get_elevation()
    elevation['data'] = np.ones((nradials, ), dtype='float32') * elev
    fixed_angle['data'] = np.array([elev], dtype='float32')

    nfile.close()
    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=None)
