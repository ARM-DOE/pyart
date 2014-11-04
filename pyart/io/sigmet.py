"""
pyart.io.sigmet
===============

Reading and writing of Sigmet (raw format) files

.. autosummary::
    :toctree: generated/

    read_sigmet
    ymds_time_to_datetime
    _time_order_data_and_metadata_full
    _time_order_data_and_metadata_roll

"""

from __future__ import division
import datetime
import warnings

import numpy as np

from ..config import FileMetadata, get_fillvalue
from ..core.radar import Radar
from .common import make_time_unit_str
from ._sigmetfile import SigmetFile, bin4_to_angle, bin2_to_angle
from . import _sigmet_noaa_hh

SPEED_OF_LIGHT = 299793000.0


def read_sigmet(filename, field_names=None, additional_metadata=None,
                file_field_names=False, exclude_fields=None,
                time_ordered='none', full_xhdr=None, noaa_hh_hdr=None,
                debug=False):
    """
    Read a Sigmet (IRIS) product file.

    Parameters
    ----------
    filename : str
        Name of Sigmet (IRIS) product file to read or file-like object
        pointing to the beginning of such a file.
    field_names : dict, optional
        Dictionary mapping Sigmet data type names to radar field names. If a
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
        True to use the Sigmet data type names for the field names. If this
        case the field_names parameter is ignored. The field dictionary will
        likely only have a 'data' key, unless the fields are defined in
        `additional_metadata`.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters.
    time_ordered : 'full', 'none' or 'roll'.
        Parameter controlling the time ordering of the data. The default,
        'none' keep the data ordered in the same manner as it appears in
        the Sigmet file.  'roll' will attempt to time order the data within
        each sweep by rolling the earliest collected ray to be the beginning.
        Sequential ordering of the rays is maintained but strict time
        increasing order is not guaranteed.  'full' will place data within
        each sweep in a strictly time increasing order, but the rays will
        likely become non-sequential.  The 'full' option is not recommended
        unless strict time increasing order is required.
    full_xhdr : bool or None
        Flag to read in all extended headers for possible decoding. None will
        determine if extended headers should be read in automatically by
        examining the extended header type.
    noaa_hh_hdr : bool or None
        Flag indicating if the extended header should be decoded as those
        used by the NOAA Hurricane Hunters aircraft radars.  None will
        determine if the extended header is of this type automatically by
        examining the header. The `full_xhdr` parameter is set to True
        when this parameter is True.
    debug : bool, optional
        Print debug information during read.

    Returns
    -------
    radar : Radar
        Radar object

    """
    # create metadata retrieval object
    filemetadata = FileMetadata('sigmet', field_names, additional_metadata,
                                file_field_names, exclude_fields)

    # open the file
    sigmetfile = SigmetFile(filename, debug=debug)
    ingest_config = sigmetfile.ingest_header['ingest_configuration']
    task_config = sigmetfile.ingest_header['task_configuration']

    # determine if full extended headers should be read
    if noaa_hh_hdr is True:
        full_xhdr = True
    if full_xhdr is None:
        type_mask = task_config['task_dsp_info']['current_data_type_mask']
        if type_mask['extended_header_type'] == 2:
            full_xhdr = True
        else:
            full_xhdr = False

    # read data
    sigmet_data, sigmet_metadata = sigmetfile.read_data(full_xhdr=full_xhdr)
    first_data_type = sigmetfile.data_type_names[0]
    if first_data_type == 'XHDR':   # don't use XHDR as the first data type
        first_data_type = sigmetfile.data_type_names[1]
    sigmetfile.close()
    nsweeps, nrays, nbins = sigmet_data[first_data_type].shape

    if nsweeps == 0:
        raise IOError('File contains no readable sweep data.')

    # parse the extended headers for time
    if full_xhdr and 'XHDR' in sigmet_data:
        # extract the ms timing data and store the full header for later
        # analysis, keep both in the sigmet_data dictionary so time ordering
        # and removal of missing rays is performed on these "fields"
        xhdr = sigmet_data.pop('XHDR')
        sigmet_data['XHDR'] = xhdr[:, :, :2].copy().view('i4')
        sigmet_data['XHDR_FULL'] = xhdr
        xhdr_metadata = {}
        for key in sigmet_metadata['XHDR'].keys():
            xhdr_metadata[key] = sigmet_metadata['XHDR'][key].copy()
        sigmet_metadata['XHDR_FULL'] = xhdr_metadata

    # time order
    if time_ordered == 'full':
        _time_order_data_and_metadata_full(sigmet_data, sigmet_metadata)
    if time_ordered == 'roll':
        _time_order_data_and_metadata_roll(sigmet_data, sigmet_metadata)

    # remove missing rays from the data
    good_rays = (sigmet_metadata[first_data_type]['nbins'] != -1)
    for field_name in sigmet_data.keys():
        sigmet_data[field_name] = sigmet_data[field_name][good_rays]
        field_metadata = sigmet_metadata[field_name]
        for key in field_metadata.keys():
            field_metadata[key] = field_metadata[key][good_rays]

    # sweep_start_ray_index and sweep_end_ray_index
    ray_count = good_rays.sum(axis=1)
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')

    ssri = np.cumsum(np.append([0], ray_count[:-1])).astype('int32')
    sweep_start_ray_index['data'] = ssri
    sweep_end_ray_index['data'] = np.cumsum(ray_count).astype('int32') - 1

    # time
    time = filemetadata('time')

    if 'XHDR' in sigmet_data:   # use time in extended headers
        tdata = sigmet_data.pop('XHDR')
        tdata = (tdata.flatten() / 1000.).astype('float64')
        sigmet_extended_header = True
    else:
        tdata = sigmet_metadata[first_data_type]['time'].astype('float64')
        sigmet_extended_header = False

    # add sweep_start time to all time values in each sweep.
    dts = [ymds_time_to_datetime(d['sweep_start_time'])
           for d in sigmetfile.ingest_data_headers[first_data_type]]
    for i, dt in enumerate(dts):
        start = sweep_start_ray_index['data'][i]
        end = sweep_end_ray_index['data'][i]
        td = (dt - dts[0])
        tdata[start:end + 1] += (td.microseconds + (td.seconds + td.days *
                                 24 * 3600) * 10**6) / 10**6
    time['data'] = tdata
    time['units'] = make_time_unit_str(dts[0])

    # _range
    _range = filemetadata('range')
    range_info = task_config['task_range_info']
    gate_0 = range_info['first_bin_range'] / 100.       # meters
    gate_nbin = range_info['last_bin_range'] / 100.     # meters
    gate_size = round((gate_nbin - gate_0) / (nbins))
    _range['data'] = gate_0 + gate_size * np.arange(nbins, dtype='float32')
    _range['meters_to_center_of_first_gate'] = np.array([gate_0],
                                                        dtype='float32')
    _range['meters_between_gates'] = np.array([gate_size], dtype='float32')

    # fields
    fields = {}
    for data_type_name, fdata in sigmet_data.iteritems():
        if data_type_name == 'XHDR_FULL':
            continue
        field_name = filemetadata.get_field_name(data_type_name)
        if field_name is None:
            continue
        field_dic = filemetadata(field_name)
        field_dic['data'] = fdata.reshape(-1, nbins)
        field_dic['_FillValue'] = get_fillvalue()
        fields[field_name] = field_dic

    # metadata
    metadata = filemetadata('metadata')
    metadata['original_container'] = 'sigmet'
    metadata['instrument_name'] = ingest_config['site_name'].strip()
    metadata['sigmet_task_name'] = (
        sigmetfile.product_hdr['product_configuration']['task_name'])

    if sigmet_extended_header:
        metadata['sigmet_extended_header'] = 'true'
    else:
        metadata['sigmet_extended_header'] = 'false'

    # scan_type
    if task_config['task_scan_info']['antenna_scan_mode'] == 2:
        scan_type = 'rhi'
    else:
        scan_type = 'ppi'

    # latitude
    latitude = filemetadata('latitude')
    lat = bin4_to_angle(ingest_config['latitude_radar'])
    if lat > 180.0:
        lat -= 360.0
    latitude['data'] = np.array([lat], dtype='float64')

    # longitude
    longitude = filemetadata('longitude')
    lon = bin4_to_angle(ingest_config['longitude_radar'])
    if lon > 180.0:
        lon -= 360.0
    longitude['data'] = np.array([lon], dtype='float64')

    # altitude
    altitude = filemetadata('altitude')
    alt = sigmetfile.product_hdr['product_end']['ground_height']
    altitude['data'] = np.array([alt], dtype='float64')

    # sweep_number
    sweep_number = filemetadata('sweep_number')
    sweep_number['data'] = np.arange(nsweeps, dtype='int32')

    # fixed_angle
    fixed_angle = filemetadata('fixed_angle')
    fa = [d['fixed_angle'] for d in
          sigmetfile.ingest_data_headers[first_data_type]]
    fixed_angle['data'] = bin2_to_angle(np.array(fa)).astype('float32')

    # azimuth
    az0 = sigmet_metadata[first_data_type]['azimuth_0']
    az1 = sigmet_metadata[first_data_type]['azimuth_1']
    el0 = sigmet_metadata[first_data_type]['elevation_0']
    el1 = sigmet_metadata[first_data_type]['elevation_1']

    azimuth = filemetadata('azimuth')
    az = (az0 + az1) / 2.
    az[np.where(np.abs(az0 - az1) > 180.0)] += 180.
    az[az > 360.0] -= 360.
    azimuth['data'] = az.astype('float32')

    # elevation
    elevation = filemetadata('elevation')
    elevation['data'] = ((el0 + el1) / 2.).astype('float32')

    # instrument_parameters
    prt = filemetadata('prt')
    prt_mode = filemetadata('prt_mode')
    nyquist_velocity = filemetadata('nyquist_velocity')
    unambiguous_range = filemetadata('unambiguous_range')
    beam_width_h = filemetadata('radar_beam_width_h')
    beam_width_v = filemetadata('radar_beam_width_v')
    pulse_width = filemetadata('pulse_width')

    trays = nsweeps * nrays     # this is correct even with missing rays
    prt_value = 1. / sigmetfile.product_hdr['product_end']['prf']
    prt['data'] = prt_value * np.ones(trays, dtype='float32')

    ur_value = SPEED_OF_LIGHT * prt_value / 2.
    unambiguous_range['data'] = ur_value * np.ones(trays, dtype='float32')

    # TODO Multi PRF mode when
    # task_config['task_dsp_info']['multi_prf_flag'] != 0
    prt_mode['data'] = np.array(nsweeps * ['fixed'])

    wavelength_cm = sigmetfile.product_hdr['product_end']['wavelength']
    nv_value = wavelength_cm / (10000.0 * 4.0 * prt_value)
    nyquist_velocity['data'] = nv_value * np.ones(trays, dtype='float32')
    beam_width_h['data'] = np.array([bin4_to_angle(
        task_config['task_misc_info']['horizontal_beamwidth'])],
        dtype='float32')
    beam_width_v['data'] = np.array([bin4_to_angle(
        task_config['task_misc_info']['vertical_beamwidth'])],
        dtype='float32')
    pulse_width['data'] = np.array(
        [task_config['task_dsp_info']['pulse_width'] * 1e-8] *
        len(time['data']), dtype='float32')

    instrument_parameters = {'unambiguous_range': unambiguous_range,
                             'prt_mode': prt_mode, 'prt': prt,
                             'nyquist_velocity': nyquist_velocity,
                             'radar_beam_width_h': beam_width_h,
                             'radar_beam_width_v': beam_width_v,
                             'pulse_width': pulse_width}

    # decode extended headers
    extended_header_params = {}
    if noaa_hh_hdr is None:
        type_mask = task_config['task_dsp_info']['current_data_type_mask']
        htype = type_mask['extended_header_type']
        if full_xhdr and htype == 2 and sigmet_data['XHDR_FULL'][0, 3] == 136:
            noaa_hh_hdr = True
        else:
            noaa_hh_hdr = False

    if noaa_hh_hdr:
        t = _sigmet_noaa_hh._decode_noaa_hh_hdr(
            sigmet_data['XHDR_FULL'], filemetadata, azimuth, elevation)
        (latitude, longitude, altitude, extended_header_params) = t
        metadata['platform_type'] = 'aircraft'
        # scan_type determined from the antenna_scan_mode parameters
        noaa_hh_scan_modes = {4: 'ppi', 7: 'rhi'}
        scan_mode = task_config['task_scan_info']['antenna_scan_mode']
        if scan_mode not in noaa_hh_scan_modes:
            warnings.warn("Unknown antenna_scan_mode, defaulting to 'rhi'.")
            scan_type = 'rhi'
        else:
            scan_type = noaa_hh_scan_modes[scan_mode]

    # sweep_mode
    sweep_mode = filemetadata('sweep_mode')
    if scan_type == 'ppi':
        sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])
    else:
        sweep_mode['data'] = np.array(nsweeps * ['rhi'])

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters,
        **extended_header_params)


def _time_order_data_and_metadata_roll(data, metadata):
    """
    Put Sigmet data and metadata in time increasing order using a single
    roll.
    """
    # Sigmet data is stored by sweep in azimuth increasing order,
    # to place in time increasing order we must roll each sweep so
    # the earliest collected ray is first.
    # This assuming all the fields have the same timing and the rays
    # were collected in sequentially, which appears to be true.

    if 'XHDR' in data:
        ref_time = data['XHDR'].copy()
        ref_time.shape = ref_time.shape[:-1]
    else:
        ref_time = metadata[metadata.keys()[0]]['time'].astype('int32')
    for i, sweep_time in enumerate(ref_time):

        # determine the number of place by which elements should be shifted.
        sweep_time_diff = np.diff(sweep_time)
        if sweep_time_diff.min() >= 0:
            continue    # already time ordered
        shift = -(sweep_time_diff.argmin() + 1)

        # roll the data and metadata for each field
        for field in data.keys():
            data[field][i] = np.roll(data[field][i], shift, axis=0)

            fmd = metadata[field]
            fmd['time'][i] = np.roll(fmd['time'][i], shift)
            fmd['nbins'][i] = np.roll(fmd['nbins'][i], shift)
            fmd['elevation_1'][i] = np.roll(fmd['elevation_1'][i], shift)
            fmd['elevation_0'][i] = np.roll(fmd['elevation_0'][i], shift)
            fmd['azimuth_0'][i] = np.roll(fmd['azimuth_0'][i], shift)
            fmd['azimuth_1'][i] = np.roll(fmd['azimuth_1'][i], shift)

    return


def _time_order_data_and_metadata_full(data, metadata):
    """
    Put Sigmet data and metadata in time increasing order by sorting the
    time.
    """

    # Sigmet data is stored by sweep in azimuth increasing order,
    # to place in time increasing order we must sort each sweep so
    # that the rays are time ordered.
    # This assuming all the fields have the same timing, rays are not
    # assumed to be collected sequentially.

    if 'XHDR' in data:
        ref_time = data['XHDR'].copy()
        ref_time.shape = ref_time.shape[:-1]
    else:
        ref_time = metadata[metadata.keys()[0]]['time'].astype('int32')
    for i, sweep_time in enumerate(ref_time):

        # determine the indices which sort the sweep time using a stable
        # sorting algorithm to prevent excessive azimuth scrambling.
        sweep_time_diff = np.diff(sweep_time)
        if sweep_time_diff.min() >= 0:
            continue    # already time ordered
        sort_idx = np.argsort(sweep_time, kind='mergesort')

        # sort the data and metadata for each field
        for field in data.keys():
            data[field][i] = data[field][i][sort_idx]

            fmd = metadata[field]
            fmd['time'][i] = fmd['time'][i][sort_idx]
            fmd['nbins'][i] = fmd['nbins'][i][sort_idx]
            fmd['elevation_1'][i] = fmd['elevation_1'][i][sort_idx]
            fmd['elevation_0'][i] = fmd['elevation_0'][i][sort_idx]
            fmd['azimuth_0'][i] = fmd['azimuth_0'][i][sort_idx]
            fmd['azimuth_1'][i] = fmd['azimuth_1'][i][sort_idx]

    return


def ymds_time_to_datetime(ymds):
    """ Return a datetime object from a Sigmet ymds_time dictionary. """
    dt = datetime.datetime(ymds['year'], ymds['month'], ymds['day'])
    delta = datetime.timedelta(seconds=ymds['seconds'])
    return dt + delta
