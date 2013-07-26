"""
pyart.io.nexrad_archive
=======================

Functions for reading NEXRAD Level II Archive files.

.. autosummary::
    :toctree: generated/

    read_nexrad_archive
    _radar_from_nexradl2
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
        UCAR THREDDS Data Server [2]_ have been tested.  Other NEXRAD Level
        2 Archibe files may or may not work.
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
    refl_hi : Radar or None
        Radar object containing the super resolution reflectivity field.
        If the file does not contain these fields None is returned, this is
        true for all returns.
    dopl_hi : Radar or None
        Radar object containing the super resolution doppler fields present.
        Possible fields in the radar are: velocity, spectrum_width,
        differential_reflectivity, differential_phase, and
        correlation_coefficient.
    refl_sd : Radar or None
        Radar object containing the standard resolution reflectivity field.
    dopl_sd : Radar or None
        Radar object containing the standard resolution doppler fields.

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

    # function for making radars res, 1: super resolution, 2: standard
    def _mr(moment, res):
        """Make a radar with given moment and res"""
        # check that at least some scans exist with the moments
        if len(scan_info[moment][res]['scans']) == 0:
            return None

        # determine scans, max_gates, dualpol_offset
        scans = scan_info[moment][res]['scans']
        max_gates = max(scan_info[moment][res]['ngates'])
        dualpol_offset = 0
        if moment == 'REF':
            moments = ['REF']
        else:
            moments = ['VEL', 'SW']
            if len(scan_info['ZDR'][res]['scans']) != 0:
                moments = ['VEL', 'SW', 'ZDR', 'PHI', 'RHO']
                dualpol_offset = (scan_info['ZDR'][res]['scans'][0] -
                                  scan_info['VEL'][res]['scans'][0])

        # read the radar from the archive file
        return _radar_from_nexradl2(nfile, moments, scans, max_gates,
                                    field_mapping, field_metadata,
                                    dualpol_offset)

    return _mr('REF', 1), _mr('VEL', 1), _mr('REF', 2), _mr('VEL', 2)


# default mapping from Archive Level 2 moment to Radar object field names
ARCHIVE_MAPPING = {
    'REF': 'reflectivity',
    'VEL': 'velocity',
    'SW': 'spectrum_width',
    'ZDR': 'differential_reflectivity',
    'PHI': 'differential_phase',
    'RHO': 'correlation_coefficient'
}


def _radar_from_nexradl2(nfile, moments, scans, max_gates, field_mapping,
                         field_metadata, dualpol_offset=0):
    """
    Create a radar object from a NEXRADLevel2File object.

    Parameters
    ----------
    nfile : NEXRADLevel2File
        NEXRADLevel2File to read radar data from.
    moments : list
        Moments to retrieve.
    scans : list
        List of scans (0 based) to retrieve.
    max_gates : int
        Maximum number of gates expected in scans.
    field_mapping : dict
        Mapping between moment and radar field names.  See
        :func:`read_nexrad_archive` for details.
    field_metadata
        Metadata for the field specified in field_mapping.  See
        :func:`read_nexrad_archive` for details.
    dualpol_offset : int, optional
        Offset of dualpol moment ('ZDR', 'RHO', and 'PHI') scans relative
        to those in the scans parameter.

    Returns
    -------
    radar : Radar
        Radar containing requested moments and scans.

    """

    # time
    time = get_metadata('time')
    time_start, _time = nfile.get_times(scans)
    time['data'] = _time
    time['units'] = make_time_unit_str(time_start)

    # _range
    _range = get_metadata('range')
    _range['data'] = nfile.get_range(scans[0], moments[0])
    _range['meters_to_center_of_first_gate'] = _range['data'][0]
    _range['meters_between_gates'] = _range['data'][1] - _range['data'][0]

    # fields
    fields = {}
    for moment in moments:
        field_name = field_mapping[moment]
        fields[field_name] = field_metadata[field_name].copy()
        fields[field_name]['_FillValue'] = -9999.0
        if moment in ['ZDR', 'RHO', 'PHI']:
            offset_scans = [s + dualpol_offset for s in scans]
            fdata = nfile.get_data(offset_scans, moment, max_gates)
        else:
            fdata = nfile.get_data(scans, moment, max_gates)
        fields[field_name]['data'] = fdata

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

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
    # sweep_end_ray_index
    sweep_number = get_metadata('sweep_number')
    sweep_mode = get_metadata('sweep_mode')
    fixed_angle = get_metadata('fixed_angle')
    sweep_start_ray_index = get_metadata('sweep_start_ray_index')
    sweep_end_ray_index = get_metadata('sweep_end_ray_index')

    nsweeps = len(scans)
    sweep_number['data'] = np.arange(len(scans), dtype='int32')
    sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])
    rays_per_scan = nfile.get_nrays(scans[0])
    ssri = np.arange(nsweeps, dtype='int32') * rays_per_scan
    sweep_start_ray_index['data'] = ssri
    sweep_end_ray_index['data'] = ssri + (rays_per_scan - 1)

    # azimuth, elevation
    azimuth = get_metadata('azimuth')
    elevation = get_metadata('elevation')
    azimuth['data'] = nfile.get_azimuth_angles(scans)
    elev = nfile.get_elevation_angles(scans)
    elevation['data'] = elev.astype('float32')
    fixed_angle['data'] = elev.astype('float32')

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=None)
