"""
pyart.io.nexrad_archive
=======================

Functions for reading NEXRAD Level II Archive files.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    _NEXRADLevel2StagedField

.. autosummary::
    :toctree: generated/

    read_nexrad_archive

"""

import numpy as np

from ..config import FileMetadata, get_fillvalue
from ..core.radar import Radar
from .common import make_time_unit_str
from .nexrad_level2 import NEXRADLevel2File
from .lazydict import LazyLoadDict


def read_nexrad_archive(filename, field_names=None, additional_metadata=None,
                        file_field_names=False, exclude_fields=None,
                        delay_field_loading=False, bzip=None):
    """
    Read a NEXRAD Level 2 Archive file.

    Parameters
    ----------
    filename : str
        Filename of NEXRAD Level 2 Archive file.  The files hosted by
        at the NOAA National Climate Data Center [1]_ as well as on the
        UCAR THREDDS Data Server [2]_ have been tested.  Other NEXRAD
        Level 2 Archive files may or may not work.  Message type 1 file
        at not yet supported, only message type 31.
    field_names : dict, optional
        Dictionary mapping NEXRAD moments to radar field names. If a
        data type found in the file does not appear in this dictionary or has
        a value of None it will not be placed in the radar.fields dictionary.
        A value of None, the default, will use the mapping defined in the
        metadata configuration file.
    additional_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve metadata from during this read.
        This metadata is not used during any successive file reads unless
        explicitly included.  A value of None, the default, will not
        introduct any addition metadata and the file specific or default
        metadata as specified by the metadata configuration file will be used.
    file_field_names : bool, optional
        True to use the NEXRAD field names for the field names. If this
        case the field_names parameter is ignored. The field dictionary will
        likely only have a 'data' key, unless the fields are defined in
        `additional_metadata`.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters.
    delay_field_loading : bool
        True to delay loading of field data from the file until the 'data'
        key in a particular field dictionary is accessed.  In this case
        the field attribute of the returned Radar object will contain
        LazyLoadDict objects not dict objects.
    bzip : bool or None
        True if the file is compressed as a bzip2 file, False otherwise.
        None will examine the filename for a bzip extension.

    Returns
    -------
    radar : Radar
        Radar object containing all moments and sweeps/cuts in the volume.
        Gates not collected are masked in the field data.

    References
    ----------
    .. [1] http://www.ncdc.noaa.gov/
    .. [2] http://thredds.ucar.edu/thredds/catalog.html

    """
    # create metadata retrieval object
    filemetadata = FileMetadata('nexrad_archive', field_names,
                                additional_metadata, file_field_names,
                                exclude_fields)
    # parse bzip parameter
    if bzip is None:
        try:
            if filename.endswith('.bz2') or filename.endswith('bzip2'):
                bzip = True
            else:
                bzip = False
        except:
            bzip = False

    # open the file and retrieve scan information
    nfile = NEXRADLevel2File(filename, bzip)
    scan_info = nfile.scan_info()

    # time
    time = filemetadata('time')
    time_start, _time = nfile.get_times()
    time['data'] = _time
    time['units'] = make_time_unit_str(time_start)

    # range
    scan_max_gates = np.argmax([max(s['ngates']) for s in scan_info])
    i = np.argmax(scan_info[scan_max_gates]['ngates'])
    moment_max_gates = scan_info[scan_max_gates]['moments'][i]
    _range = filemetadata('range')
    _range['data'] = nfile.get_range(scan_max_gates, moment_max_gates)
    _range['meters_to_center_of_first_gate'] = _range['data'][0]
    _range['meters_between_gates'] = _range['data'][1] - _range['data'][0]

    # fields
    max_ngates = max([max(s['ngates']) for s in scan_info])
    available_moments = set([])
    for info in scan_info:
        for moment in info['moments']:
            available_moments.add(moment)

    fields = {}
    for moment in available_moments:
        field_name = filemetadata.get_field_name(moment)
        if field_name is None:
            continue
        dic = filemetadata(field_name)
        dic['_FillValue'] = get_fillvalue()
        if delay_field_loading:
            dic = LazyLoadDict(dic)
            data_call = _NEXRADLevel2StagedField(nfile, moment, max_ngates)
            dic.set_lazy('data', data_call)
        else:
            dic['data'] = nfile.get_data(moment, max_ngates)
        fields[field_name] = dic

    # metadata
    metadata = filemetadata('metadata')
    metadata['original_container'] = 'NEXRAD Level II'

    # scan_type
    scan_type = 'ppi'

    # latitude, longitude, altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')

    lat, lon, alt = nfile.location()
    latitude['data'] = np.array([lat], dtype='float64')
    longitude['data'] = np.array([lon], dtype='float64')
    altitude['data'] = np.array([alt], dtype='float64')

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index
    # sweep_end_ray_index
    sweep_number = filemetadata('sweep_number')
    sweep_mode = filemetadata('sweep_mode')
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')

    nsweeps = int(nfile.nscans)
    sweep_number['data'] = np.arange(nsweeps, dtype='int32')
    sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])

    rays_per_scan = [s['nrays'] for s in scan_info]
    sweep_end_ray_index['data'] = np.cumsum(rays_per_scan, dtype='int32') - 1

    rays_per_scan.insert(0, 0)
    sweep_start_ray_index['data'] = np.cumsum(rays_per_scan[:-1],
                                              dtype='int32')

    # azimuth, elevation, fixed_angle
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')
    fixed_angle = filemetadata('fixed_angle')
    azimuth['data'] = nfile.get_azimuth_angles()
    elevation['data'] = nfile.get_elevation_angles().astype('float32')
    fixed_angle['data'] = nfile.get_target_angles()

    # instrument_parameters
    nyquist_velocity = filemetadata('nyquist_velocity')
    unambiguous_range = filemetadata('unambiguous_range')
    nyquist_velocity['data'] = nfile.get_nyquist_vel().astype('float32')
    unambiguous_range['data'] = nfile.get_unambigous_range().astype('float32')

    instrument_parameters = {'unambiguous_range': unambiguous_range,
                             'nyquist_velocity': nyquist_velocity, }

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)


class _NEXRADLevel2StagedField(object):
    """
    A class to facilitate on demand loading of field data from a Level 2 file.
    """

    def __init__(self, nfile, moment, max_ngates):
        """ initialize. """
        self.nfile = nfile
        self.moment = moment
        self.max_ngates = max_ngates

    def __call__(self):
        """ Return the array containing the field data. """
        return self.nfile.get_data(self.moment, self.max_ngates)
