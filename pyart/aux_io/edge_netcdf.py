"""
Utilities for reading EDGE NetCDF files.

"""

import datetime

import numpy as np
import netCDF4

from ..config import FileMetadata, get_fillvalue
from ..io.common import make_time_unit_str, _test_arguments
from ..core.radar import Radar


def read_edge_netcdf(filename, **kwargs):
    """
    Read a EDGE NetCDF file.

    Parameters
    ----------
    filename : str
        Name of EDGE NetCDF file to read data from.

    Returns
    -------
    radar : Radar
        Radar object.

    """
    # test for non empty kwargs
    _test_arguments(kwargs)

    # create metadata retrieval object
    filemetadata = FileMetadata('edge_netcdf')

    # Open netCDF4 file
    dset = netCDF4.Dataset(filename)
    nrays = len(dset.dimensions['Azimuth'])
    nbins = len(dset.dimensions['Gate'])

    # latitude, longitude and altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')
    latitude['data'] = np.array([dset.Latitude], 'float64')
    longitude['data'] = np.array([dset.Longitude], 'float64')
    altitude['data'] = np.array([dset.Height], 'float64')

    # metadata
    metadata = filemetadata('metadata')
    metadata_mapping = {
        'vcp-value': 'vcp',
        'radarName-value': 'radar_name',
        'ConversionPlugin': 'conversion_software',
    }
    for netcdf_attr, metadata_key in metadata_mapping.items():
        if netcdf_attr in dset.ncattrs():
            metadata[metadata_key] = dset.getncattr(netcdf_attr)

    # sweep_start_ray_index, sweep_end_ray_index
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')
    sweep_start_ray_index['data'] = np.array([0], dtype='int32')
    sweep_end_ray_index['data'] = np.array([nrays-1], dtype='int32')

    # sweep number
    sweep_number = filemetadata('sweep_number')
    sweep_number['data'] = np.array([0], dtype='int32')

    # sweep_type
    scan_type = 'ppi'

    # sweep_mode, fixed_angle
    sweep_mode = filemetadata('sweep_mode')
    fixed_angle = filemetadata('fixed_angle')
    sweep_mode['data'] = np.array(1 * ['azimuth_surveillance'])
    fixed_angle['data'] = np.array([dset.Elevation], dtype='float32')

    # time
    time = filemetadata('time')
    start_time = datetime.datetime.utcfromtimestamp(dset.Time)
    time['units'] = make_time_unit_str(start_time)
    time['data'] = np.zeros((nrays, ), dtype='float64')

    # range
    _range = filemetadata('range')
    step = float(dset.getncattr('MaximumRange-value')) / nbins * 1000.
    _range['data'] = (np.arange(nbins, dtype='float32') * step + step / 2)
    _range['meters_to_center_of_first_gate'] = step / 2.
    _range['meters_between_gates'] = step

    # elevation
    elevation = filemetadata('elevation')
    elevation_angle = dset.Elevation
    elevation['data'] = np.ones((nrays, ), dtype='float32') * elevation_angle

    # azimuth
    azimuth = filemetadata('azimuth')
    azimuth['data'] = dset.variables['Azimuth'][:]

    # fields
    field_name = dset.TypeName

    field_data = np.ma.array(dset.variables[field_name][:])
    if 'MissingData' in dset.ncattrs():
        field_data[field_data == dset.MissingData] = np.ma.masked
    if 'RangeFolded' in dset.ncattrs():
        field_data[field_data == dset.RangeFolded] = np.ma.masked

    fields = {field_name: filemetadata(field_name)}
    fields[field_name]['data'] = field_data
    fields[field_name]['units'] = dset.variables[field_name].Units
    fields[field_name]['_FillValue'] = get_fillvalue()

    # instrument_parameters
    instrument_parameters = {}

    if 'PRF-value' in dset.ncattrs():
        dic = filemetadata('prt')
        prt = 1. / float(dset.getncattr('PRF-value'))
        dic['data'] = np.ones((nrays, ), dtype='float32') * prt
        instrument_parameters['prt'] = dic

    if 'PulseWidth-value' in dset.ncattrs():
        dic = filemetadata('pulse_width')
        pulse_width = dset.getncattr('PulseWidth-value') * 1.e-6
        dic['data'] = np.ones((nrays, ), dtype='float32') * pulse_width
        instrument_parameters['pulse_width'] = dic

    if 'NyquistVelocity-value' in dset.ncattrs():
        dic = filemetadata('nyquist_velocity')
        nyquist_velocity = float(dset.getncattr('NyquistVelocity-value'))
        dic['data'] = np.ones((nrays, ), dtype='float32') * nyquist_velocity
        instrument_parameters['nyquist_velocity'] = dic

    if 'Beamwidth' in dset.variables:
        dic = filemetadata('radar_beam_width_h')
        dic['data'] = dset.variables['Beamwidth'][:]
        instrument_parameters['radar_beam_width_h'] = dic

    dset.close()

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)
