"""
pyart.io.nexrad_archive
=======================

Functions for reading NEXRAD Level II Archive files.

.. autosummary::
    :toctree: generated/

    read_nexrad_archive
    ARCHIVE_MAPPING
    NEXRAD_METADATA

"""

import numpy as np

from .radar import Radar
from .common import get_metadata, make_time_unit_str
from .nexrad_common import NEXRAD_METADATA
from .nexrad_level2 import NEXRADLevel2File


def read_nexrad_archive(filename, bzip=None, field_mapping=None,
                        field_metadata=None):
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
    bzip : bool or None
        True if the file is compressed as a bzip2 file, False otherwise.
        None will examine the filename for a bzip extension.
    field_mapping : None or dict, optional
        Dictionary mapping NEXRAD moments to the corresponding field names in
        the radar objects returned. None will use :data:`ARCHIVE_MAPPING`.
        which also served as an example of the format for this parameter.
        If a dictionary parameter is used it must have the same dictionary
        keys as ARCHIVE_MAPPING.  In addition, field_metadata must also be
        provided which contains the field metadata for the fields specified.
    field_metadata : None or dict, optional
        Metadata for the fields specified by field_mapping, None will use the
        field metadata provided in :data:`NEXRAD_METADATA`, which also serves
        as an example of the format for this parameter.  This metadata will
        be used for the field in the created radar objects returned.

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
    # parse the parameters
    if field_mapping is None:
        field_mapping = ARCHIVE_MAPPING.copy()
    if field_metadata is None:
        field_metadata = NEXRAD_METADATA.copy()

    if bzip is None:
        if filename.endswith('.bz2') or filename.endswith('bzip2'):
            bzip = True
        else:
            bzip = False

    # open the file and retrieve scan information
    nfile = NEXRADLevel2File(filename, bzip)
    scan_info = nfile.scan_info()

    # time
    time = get_metadata('time')
    time_start, _time = nfile.get_times()
    time['data'] = _time
    time['units'] = make_time_unit_str(time_start)

    # range
    scan_max_gates = np.argmax([max(s['ngates']) for s in scan_info])
    i = np.argmax(scan_info[scan_max_gates]['ngates'])
    moment_max_gates = scan_info[scan_max_gates]['moments'][i]
    _range = get_metadata('range')
    _range['data'] = nfile.get_range(scan_max_gates, moment_max_gates)
    _range['meters_to_center_of_first_gate'] = _range['data'][0]
    _range['meters_between_gates'] = _range['data'][1] - _range['data'][0]

    # fields

    # determine shape of field data
    ngates = max([max(s['ngates']) for s in scan_info])
    nrays = sum([s['nrays'] for s in scan_info])
    field_shape = (nrays, ngates)
    nexrad_moments = ['REF', 'VEL', 'SW', 'ZDR', 'PHI', 'RHO']

    # create empty field dictionary
    fields = {}
    for moment in nexrad_moments:
        field_name = field_mapping[moment]
        dic = field_metadata[field_name].copy()
        dic['_FillValue'] = -9999.0
        dic['data'] = np.ma.masked_all(field_shape, dtype='float32')
        fields[field_name] = dic

    # fill in data
    for i, msg in enumerate(nfile.msg31s):
        for moment in [m for m in msg.keys() if m in nexrad_moments]:
            offset = msg[moment]['offset']
            scale = msg[moment]['scale']
            ngates = msg[moment]['ngates']
            data = msg[moment]['data']
            mdata = (np.ma.masked_less_equal(data, 1) - offset) / (scale)
            field_name = field_mapping[moment]
            fields[field_name]['data'][i, :ngates] = mdata
            #fields[moment]['data'][i, :ngates] = mdata

    # metadata
    metadata = {'original_container': 'NEXRAD Level II'}
    # additional required CF/Radial metadata set to blank strings
    metadata['title'] = ''
    metadata['institution'] = ''
    metadata['references'] = ''
    metadata['source'] = ''
    metadata['comment'] = ''
    metadata['instrument_name'] = ''

    # scan_type
    scan_type = 'ppi'

    # latitude, longitude, altitude
    latitude = get_metadata('latitude')
    longitude = get_metadata('longitude')
    altitude = get_metadata('altitude')

    lat, lon, alt = nfile.location()
    latitude['data'] = np.array([lat], dtype='float64')
    longitude['data'] = np.array([lon], dtype='float64')
    altitude['data'] = np.array([alt], dtype='float64')

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index
    # sweep_end_ray_index
    sweep_number = get_metadata('sweep_number')
    sweep_mode = get_metadata('sweep_mode')
    sweep_start_ray_index = get_metadata('sweep_start_ray_index')
    sweep_end_ray_index = get_metadata('sweep_end_ray_index')

    nsweeps = int(nfile.nscans)
    sweep_number['data'] = np.arange(nsweeps, dtype='int32')
    sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])

    rays_per_scan = [s['nrays'] for s in scan_info]
    sweep_end_ray_index['data'] = np.cumsum(rays_per_scan, dtype='int32') - 1

    rays_per_scan.insert(0, 0)
    sweep_start_ray_index['data'] = np.cumsum(rays_per_scan[:-1],
                                              dtype='int32')

    # azimuth, elevation, fixed_angle
    azimuth = get_metadata('azimuth')
    elevation = get_metadata('elevation')
    fixed_angle = get_metadata('fixed_angle')
    azimuth['data'] = nfile.get_azimuth_angles()
    elevation['data'] = nfile.get_elevation_angles().astype('float32')
    fixed_angle['data'] = nfile.get_target_angles()

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=None)


# default mapping from Archive Level 2 moment to Radar object field names
ARCHIVE_MAPPING = {
    'REF': 'reflectivity',
    'VEL': 'velocity',
    'SW': 'spectrum_width',
    'ZDR': 'differential_reflectivity',
    'PHI': 'differential_phase',
    'RHO': 'correlation_coefficient'
}
