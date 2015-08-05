"""
pyart.aux_io.pattern
====================

Routines for reading files from the X-band radar from the PATTERN_ project.

.. _PATTERN: http://www.mi.uni-hamburg.de/PATTERN-Pre.6763.0.html


.. autosummary::
    :toctree: generated/

    read_pattern

"""

import datetime

import numpy as np
import netCDF4

from ..config import FileMetadata
from ..io.common import make_time_unit_str, _test_arguments
from ..core.radar import Radar


def read_pattern(filename, **kwargs):
    """
    Read a netCDF file from a PATTERN project X-band radar.

    Parameters
    ----------
    filename : str
        Name of netCDF file to read data from.

    Returns
    -------
    radar : Radar
        Radar object.

    """
    # test for non empty kwargs
    _test_arguments(kwargs)

    # create metadata retrieval object
    filemetadata = FileMetadata('pattern')

    # read the data
    ncobj = netCDF4.Dataset(filename)
    ncvars = ncobj.variables

    # general parameters
    nrays = ncvars['Azimuth'].shape[0]
    scan_type = 'ppi'

    # time
    # interpolate between the first and last timestamps in the Time variable
    time = filemetadata('time')
    nctime = ncvars['Time']
    time['units'] = make_time_unit_str(
        datetime.datetime.utcfromtimestamp(nctime[0]))
    time['data'] = np.linspace(0, nctime[-1] - nctime[0], nrays)

    # range
    _range = filemetadata('range')
    _range['data'] = ncvars['Distance'][:]
    _range['meters_to_center_of_first_gate'] = _range['data'][0]
    # assuming the distance between all gates is constant, may not
    # always be true.
    _range['meters_between_gates'] = (_range['data'][1] - _range['data'][0])

    # fields
    # files contain a single corrected reflectivity field
    fields = {}
    field_name = filemetadata.get_field_name('corrected_reflectivity')
    field_dic = filemetadata(field_name)
    field_dic['_FillValue'] = ncvars['Corrected_Reflectivity']._FillValue
    field_dic['data'] = ncvars['Corrected_Reflectivity'][:]
    fields[field_name] = field_dic

    # metadata
    metadata = filemetadata('metadata')
    for k in ['institution', 'title', 'used_algorithms']:
        if k in ncobj.ncattrs():
            metadata[k] = ncobj.getncattr(k)

    # latitude, longitude, altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')
    latitude['data'] = np.array([ncobj.latitude[:-1]], dtype='float64')
    longitude['data'] = np.array([ncobj.longitude[:-1]], dtype='float64')
    altitude['data'] = np.array([ncobj.elevation], dtype='float64')

    # sweep parameters
    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
    # sweep_end_ray_index
    sweep_number = filemetadata('sweep_number')
    sweep_mode = filemetadata('sweep_mode')
    fixed_angle = filemetadata('fixed_angle')
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')

    sweep_number['data'] = np.arange(1, dtype='int32')
    sweep_mode['data'] = np.array(1 * ['azimuth_surveillance'])
    fixed_angle['data'] = np.array([0], dtype='float32')
    sweep_start_ray_index['data'] = np.array([0], dtype='int32')
    sweep_end_ray_index['data'] = np.array([nrays-1], dtype='int32')

    # azimuth, elevation
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')

    azimuth['data'] = ncvars['Azimuth'][:]
    elevation['data'] = np.array([0.], dtype='float32')

    # instrument parameters
    instrument_parameters = None

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)
