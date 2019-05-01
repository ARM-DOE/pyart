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
    _find_range_params
    _find_scans_to_interp
    _interpolate_scan

"""

import warnings

import numpy as np

from ..config import FileMetadata, get_fillvalue
from ..core.radar import Radar
from .common import make_time_unit_str, _test_arguments, prepare_for_read
from .nexrad_level2 import NEXRADLevel2File
from ..lazydict import LazyLoadDict
from .nexrad_common import get_nexrad_location
from .nexrad_interpolate import _fast_interpolate_scan


def read_nexrad_archive(filename, field_names=None, additional_metadata=None,
                        file_field_names=False, exclude_fields=None,
                        include_fields=None, delay_field_loading=False,
                        station=None, scans=None,
                        linear_interp=True, **kwargs):
    """
    Read a NEXRAD Level 2 Archive file.

    Parameters
    ----------
    filename : str
        Filename of NEXRAD Level 2 Archive file. The files hosted by
        at the NOAA National Climate Data Center [1]_ as well as on the
        UCAR THREDDS Data Server [2]_ have been tested. Other NEXRAD
        Level 2 Archive files may or may not work. Message type 1 file
        and message type 31 files are supported.
    field_names : dict, optional
        Dictionary mapping NEXRAD moments to radar field names. If a
        data type found in the file does not appear in this dictionary or has
        a value of None it will not be placed in the radar.fields dictionary.
        A value of None, the default, will use the mapping defined in the
        metadata configuration file.
    additional_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve metadata from during this read.
        This metadata is not used during any successive file reads unless
        explicitly included. A value of None, the default, will not
        introduct any addition metadata and the file specific or default
        metadata as specified by the metadata configuration file will be used.
    file_field_names : bool, optional
        True to use the NEXRAD field names for the field names. If this
        case the field_names parameter is ignored. The field dictionary will
        likely only have a 'data' key, unless the fields are defined in
        `additional_metadata`.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters. Set
        to None to include all fields specified by include_fields.
    include_fields : list or None, optional
        List of fields to include from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters. Set
        to None to include all fields not specified by exclude_fields.
    delay_field_loading : bool, optional
        True to delay loading of field data from the file until the 'data'
        key in a particular field dictionary is accessed. In this case
        the field attribute of the returned Radar object will contain
        LazyLoadDict objects not dict objects.
    station : str or None, optional
        Four letter ICAO name of the NEXRAD station used to determine the
        location in the returned radar object. This parameter is only
        used when the location is not contained in the file, which occur
        in older NEXRAD message 1 files.
    scans : list or None, optional
        Read only specified scans from the file. None (the default) will read
        all scans.
    linear_interp : bool, optional
        True (the default) to perform linear interpolation between valid pairs
        of gates in low resolution rays in files mixed resolution rays.
        False will perform a nearest neighbor interpolation. This parameter is
        not used if the resolution of all rays in the file or requested sweeps
        is constant.

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
    # test for non empty kwargs
    _test_arguments(kwargs)

    # create metadata retrieval object
    filemetadata = FileMetadata('nexrad_archive', field_names,
                                additional_metadata, file_field_names,
                                exclude_fields, include_fields)

    # open the file and retrieve scan information
    nfile = NEXRADLevel2File(prepare_for_read(filename))
    scan_info = nfile.scan_info(scans)

    # time
    time = filemetadata('time')
    time_start, _time = nfile.get_times(scans)
    time['data'] = _time
    time['units'] = make_time_unit_str(time_start)

    # range
    _range = filemetadata('range')
    first_gate, gate_spacing, last_gate = _find_range_params(
        scan_info, filemetadata)
    _range['data'] = np.arange(first_gate, last_gate, gate_spacing, 'float32')
    _range['meters_to_center_of_first_gate'] = float(first_gate)
    _range['meters_between_gates'] = float(gate_spacing)

    # metadata
    metadata = filemetadata('metadata')
    metadata['original_container'] = 'NEXRAD Level II'
    vcp_pattern = nfile.get_vcp_pattern()
    if vcp_pattern is not None:
        metadata['vcp_pattern'] = vcp_pattern
    if 'icao' in nfile.volume_header.keys():
        metadata['instrument_name'] = nfile.volume_header['icao'].decode()

    # scan_type
    scan_type = 'ppi'

    # latitude, longitude, altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')

    if nfile._msg_type == '1' and station is not None:
        lat, lon, alt = get_nexrad_location(station)
    else:
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

    if scans is None:
        nsweeps = int(nfile.nscans)
    else:
        nsweeps = len(scans)
    sweep_number['data'] = np.arange(nsweeps, dtype='int32')
    sweep_mode['data'] = np.array(
        nsweeps * ['azimuth_surveillance'], dtype='S')

    rays_per_scan = [s['nrays'] for s in scan_info]
    sweep_end_ray_index['data'] = np.cumsum(rays_per_scan, dtype='int32') - 1

    rays_per_scan.insert(0, 0)
    sweep_start_ray_index['data'] = np.cumsum(
        rays_per_scan[:-1], dtype='int32')

    # azimuth, elevation, fixed_angle
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')
    fixed_angle = filemetadata('fixed_angle')
    azimuth['data'] = nfile.get_azimuth_angles(scans)
    elevation['data'] = nfile.get_elevation_angles(scans).astype('float32')
    fixed_agl = []
    for i in nfile.get_target_angles(scans):
        if i > 180:
            i = i - 360.
            warnings.warn("Fixed_angle(s) greater than 180 degrees present."
                          + " Assuming angle to be negative so subtrating 360",
                          UserWarning)
        else:
            i = i
        fixed_agl.append(i)
    fixed_angles = np.array(fixed_agl, dtype='float32')
    fixed_angle['data'] = fixed_angles

    # fields
    max_ngates = len(_range['data'])
    available_moments = set([m for scan in scan_info for m in scan['moments']])
    interpolate = _find_scans_to_interp(
        scan_info, first_gate, gate_spacing, filemetadata)

    fields = {}
    for moment in available_moments:
        field_name = filemetadata.get_field_name(moment)
        if field_name is None:
            continue
        dic = filemetadata(field_name)
        dic['_FillValue'] = get_fillvalue()
        if delay_field_loading and moment not in interpolate:
            dic = LazyLoadDict(dic)
            data_call = _NEXRADLevel2StagedField(
                nfile, moment, max_ngates, scans)
            dic.set_lazy('data', data_call)
        else:
            mdata = nfile.get_data(moment, max_ngates, scans=scans)
            if moment in interpolate:
                interp_scans = interpolate[moment]
                warnings.warn(
                    "Gate spacing is not constant, interpolating data in " +
                    "scans %s for moment %s." % (interp_scans, moment),
                    UserWarning)
                for scan in interp_scans:
                    idx = scan_info[scan]['moments'].index(moment)
                    moment_ngates = scan_info[scan]['ngates'][idx]
                    start = sweep_start_ray_index['data'][scan]
                    end = sweep_end_ray_index['data'][scan]
                    _interpolate_scan(mdata, start, end, moment_ngates,
                                      linear_interp)
            dic['data'] = mdata
        fields[field_name] = dic

    # instrument_parameters
    nyquist_velocity = filemetadata('nyquist_velocity')
    unambiguous_range = filemetadata('unambiguous_range')
    nyquist_velocity['data'] = nfile.get_nyquist_vel(scans).astype('float32')
    unambiguous_range['data'] = (
        nfile.get_unambigous_range(scans).astype('float32'))

    instrument_parameters = {'unambiguous_range': unambiguous_range,
                             'nyquist_velocity': nyquist_velocity, }

    nfile.close()
    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)


def _find_range_params(scan_info, filemetadata):
    """ Return range parameters, first_gate, gate_spacing, last_gate. """
    min_first_gate = 999999
    min_gate_spacing = 999999
    max_last_gate = 0
    for scan_params in scan_info:
        ngates = scan_params['ngates'][0]
        for i, moment in enumerate(scan_params['moments']):
            if filemetadata.get_field_name(moment) is None:
                # moment is not read, skip
                continue
            first_gate = scan_params['first_gate'][i]
            gate_spacing = scan_params['gate_spacing'][i]
            last_gate = first_gate + gate_spacing * (ngates - 0.5)

            min_first_gate = min(min_first_gate, first_gate)
            min_gate_spacing = min(min_gate_spacing, gate_spacing)
            max_last_gate = max(max_last_gate, last_gate)
    return min_first_gate, min_gate_spacing, max_last_gate


def _find_scans_to_interp(scan_info, first_gate, gate_spacing, filemetadata):
    """ Return a dict indicating what moments/scans need interpolation.  """
    moments = set([m for scan in scan_info for m in scan['moments']])
    interpolate = dict([(moment, []) for moment in moments])
    for scan_num, scan in enumerate(scan_info):
        for moment in moments:
            if moment not in scan['moments']:
                continue
            if filemetadata.get_field_name(moment) is None:
                # moment is not read, skip
                continue
            index = scan['moments'].index(moment)
            first = scan['first_gate'][index]
            spacing = scan['gate_spacing'][index]
            if first != first_gate or spacing != gate_spacing:
                interpolate[moment].append(scan_num)
                # for proper interpolation the gate spacing of the scan to be
                # interpolated should be 1/4th the spacing of the radar
                assert spacing == gate_spacing * 4
                # and the first gate for the scan should be one and half times
                # the radar spacing past the radar first gate
                assert first_gate + 1.5 * gate_spacing == first
    # remove moments with no scans needing interpolation
    interpolate = dict([(k, v) for k, v in interpolate.items() if len(v) != 0])
    return interpolate


def _interpolate_scan(mdata, start, end, moment_ngates, linear_interp=True):
    """ Interpolate a single NEXRAD moment scan from 1000 m to 250 m. """
    fill_value = -9999
    data = mdata.filled(fill_value)
    scratch_ray = np.empty((data.shape[1], ), dtype=data.dtype)
    _fast_interpolate_scan(data, scratch_ray, fill_value,
                           start, end, moment_ngates, linear_interp)
    mdata[:] = np.ma.array(data, mask=(data == fill_value))


class _NEXRADLevel2StagedField(object):
    """
    A class to facilitate on demand loading of field data from a Level 2 file.
    """

    def __init__(self, nfile, moment, max_ngates, scans):
        """ initialize. """
        self.nfile = nfile
        self.moment = moment
        self.max_ngates = max_ngates
        self.scans = scans

    def __call__(self):
        """ Return the array containing the field data. """
        return self.nfile.get_data(
            self.moment, self.max_ngates, scans=self.scans)
