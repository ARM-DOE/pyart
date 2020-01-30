"""
Reading and writing of Sigmet (raw format) files

"""

import datetime
import warnings

import numpy as np

from ..config import FileMetadata, get_fillvalue
from ..core.radar import Radar
from .common import make_time_unit_str, _test_arguments, prepare_for_read
from ._sigmetfile import SigmetFile, bin4_to_angle, bin2_to_angle
from . import _sigmet_noaa_hh
from ..util import mean_of_two_angles_deg

SPEED_OF_LIGHT = 299793000.0


def read_sigmet(filename, field_names=None, additional_metadata=None,
                file_field_names=False, exclude_fields=None,
                include_fields=None, time_ordered='none', full_xhdr=None,
                noaa_hh_hdr=None, debug=False, ignore_xhdr=False,
                ignore_sweep_start_ms=None, **kwargs):
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
        explicitly included. A value of None, the default, will not
        introduct any addition metadata and the file specific or default
        metadata as specified by the metadata configuration file will be used.
    file_field_names : bool, optional
        True to use the Sigmet data type names for the field names. If this
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
    time_ordered : 'none', 'sequential', 'full', ...,  optional
        Parameter controlling if and how the rays are re-ordered by time.
        The default, 'none' keeps the rays ordered in the same manner as
        they appears in the Sigmet file. 'sequential' will determind and
        apply an operation which maintains a sequential ray order in elevation
        or azimuth yet orders the rays according to time. If no operation can
        be found to accomplish this a warning is issue and the rays are
        returned in their original order. 'roll', 'reverse', and
        'reverse_and_roll' will apply that operation to the rays in order to
        place them in time order, direct use of these is not recommended.
        'full' will order the rays in strictly time increasing order,
        but the rays will likely become non-sequential, thisoption is not
        recommended unless strict time increasing order is required.
    full_xhdr : bool or None
        Flag to read in all extended headers for possible decoding. None will
        determine if extended headers should be read in automatically by
        examining the extended header type.
    noaa_hh_hdr : bool or None
        Flag indicating if the extended header should be decoded as those
        used by the NOAA Hurricane Hunters aircraft radars. None will
        determine if the extended header is of this type automatically by
        examining the header. The `full_xhdr` parameter is set to True
        when this parameter is True.
    ignore_xhdr : bool, optional
        True to ignore all data in the extended headers if they exist.
        False, the default, extracts milliseconds precision times and other
        parameter from the extended headers if they exists in the file.
    ignore_sweep_start_ms : bool or None, optional
        True to ignore the millisecond parameter in the start time for each
        sweep, False will uses this parameter when determining the timing of
        each ray. None, the default, will ignore the millisecond sweep start
        timing only when the file does not contain extended headers or when
        the extended header has been explicity ignored using the `ignore_xhdr`
        parameter. The TRMM RSL library ignores these times so setting this
        parameter to True is required to match the times determined when
        reading Sigmet files with :py:func:`pyart.io.read_rsl`.
        When there are not extended headers ignoring the millisecond sweep
        times provides time data which is always prior to the actual
        collection time with an error from 0 to 2 seconds.
    debug : bool, optional
        Print debug information during read.

    Returns
    -------
    radar : Radar
        Radar object.

    """
    # test for non empty kwargs
    _test_arguments(kwargs)

    # create metadata retrieval object
    filemetadata = FileMetadata('sigmet', field_names, additional_metadata,
                                file_field_names, exclude_fields,
                                include_fields)

    # open the file
    sigmetfile = SigmetFile(prepare_for_read(filename), debug=debug)
    ingest_config = sigmetfile.ingest_header['ingest_configuration']
    task_config = sigmetfile.ingest_header['task_configuration']

    # determine if full extended headers should be read
    if noaa_hh_hdr:
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

    # ignore extended header if user requested
    if ignore_xhdr:
        if 'XHDR' in sigmet_data:
            sigmet_data.pop('XHDR')
            sigmet_metadata.pop('XHDR')

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

    # remove missing rays from the data
    good_rays = (sigmet_metadata[first_data_type]['nbins'] != -1)
    rays_missing = (sigmet_metadata[first_data_type]['nbins'] == -1).sum()
    for field_name in sigmet_data.keys():
        sigmet_data[field_name] = sigmet_data[field_name][good_rays]
        field_metadata = sigmet_metadata[field_name]
        for key in field_metadata.keys():
            field_metadata[key] = field_metadata[key][good_rays]
    rays_per_sweep = good_rays.sum(axis=1)

    # time order
    if time_ordered == 'sequential':
        # Determine appropiate time ordering operation for the sweep and
        # determine if that operation will in fact order the times.
        # If it does not issue a warning and perform no time ordering
        if task_config['task_scan_info']['antenna_scan_mode'] == 2:
            # RHI scan
            if _is_time_ordered_by_reversal(
                    sigmet_data, sigmet_metadata, rays_per_sweep):
                time_ordered = 'reverse'
            else:
                warnings.warn('Rays not collected sequentially in time.')
                time_ordered = 'none'
        else:
            # PPI scan
            if _is_time_ordered_by_roll(
                    sigmet_data, sigmet_metadata, rays_per_sweep):
                time_ordered = 'roll'
            elif _is_time_ordered_by_reverse_roll(
                    sigmet_data, sigmet_metadata, rays_per_sweep):
                time_ordered = 'reverse_and_roll'
            else:
                warnings.warn('Rays not collected sequentially in time.')
                time_ordered = 'none'

    if time_ordered == 'full':
        _time_order_data_and_metadata_full(
            sigmet_data, sigmet_metadata, rays_per_sweep)
    if time_ordered == 'reverse':
        _time_order_data_and_metadata_reverse(
            sigmet_data, sigmet_metadata, rays_per_sweep)
    if time_ordered == 'roll':
        _time_order_data_and_metadata_roll(
            sigmet_data, sigmet_metadata, rays_per_sweep)
    if time_ordered == 'reverse_and_roll':
        # reverse followed by roll
        _time_order_data_and_metadata_reverse(
            sigmet_data, sigmet_metadata, rays_per_sweep)
        _time_order_data_and_metadata_roll(
            sigmet_data, sigmet_metadata, rays_per_sweep)

    # sweep_start_ray_index and sweep_end_ray_index
    ray_count = good_rays.sum(axis=1)
    total_rays = ray_count.sum()
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')

    ssri = np.cumsum(np.append([0], ray_count[:-1])).astype('int32')
    sweep_start_ray_index['data'] = ssri
    sweep_end_ray_index['data'] = np.cumsum(ray_count).astype('int32') - 1

    # time
    # The time of each ray collected in a Sigmet file is stored in the ray
    # header as an offset from the start of of the sweep. If the file contains
    # extended ray headers then the time is truncated to millisecond
    # precision. Without extended headers the time is truncated to second
    # precision. The time of the sweep start is recorded in the ingest data
    # header to millisecond percision.
    # The method below sets the volume starting time to the time of the first
    # sweep truncated to the nearest second so that values time['data'] are
    # always positive. The timing of each ray is then calculated from the
    # offset of the time of the sweep from this volume starting time
    # and the offset from the sweep start to the ray recorded in the ray
    # header. Additional, the user can select to trucate the sweep start
    # times to second precision using the `ignore_sweep_start_ms` parameter.
    # This may be desired to match libraries such as TRMM RSL, which do not
    # read the millisecond sweep start precision. By default this
    # trucation is applied when the file contains no extended header. This
    # ensures that the times listed occur BEFORE the ray was collected with an
    # error from 0 to 2 seconds from the recorded time. Other methods of
    # combining the mixed percision times would result in times which may have
    # a negative error which is undesired.

    # extract times from the ray headers
    if 'XHDR' in sigmet_data:   # use time in extended headers
        tdata = sigmet_data.pop('XHDR')
        tdata = (tdata.flatten() / 1000.).astype('float64')
        sigmet_extended_header = True
    else:   # use times from the standard ray headers
        tdata = sigmet_metadata[first_data_type]['time'].astype('float64')
        sigmet_extended_header = False

    # determine sweep start times and possibly trucate them to sec. precision
    dts = [ymds_time_to_datetime(d['sweep_start_time'])
           for d in sigmetfile.ingest_data_headers[first_data_type]]
    if ignore_sweep_start_ms is None:
        ignore_sweep_start_ms = not sigmet_extended_header
    if ignore_sweep_start_ms:
        dts = [d.replace(microsecond=0) for d in dts]

    # truncate volume start time to second precision of first sweep
    dt_start = dts[0].replace(microsecond=0)

    # add sweep start times to all time values in each sweep.
    for i, dt in enumerate(dts):
        start = sweep_start_ray_index['data'][i]
        end = sweep_end_ray_index['data'][i]
        td = (dt - dt_start)    # time delta from volume start
        tdata[start:end + 1] += (
            td.microseconds + (
                td.seconds + td.days * 24 * 3600) * 10**6) / 10**6

    # stores this data in the time dictionary
    time = filemetadata('time')
    time['data'] = tdata
    time['units'] = make_time_unit_str(dt_start)

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
    for data_type_name, fdata in sigmet_data.items():
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
    metadata['time_ordered'] = time_ordered
    metadata['rays_missing'] = rays_missing

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
    azimuth = filemetadata('azimuth')
    az0 = sigmet_metadata[first_data_type]['azimuth_0']
    az1 = sigmet_metadata[first_data_type]['azimuth_1']
    az_data = mean_of_two_angles_deg(az0, az1).astype('float32')
    az_data[az_data < 0] += 360.0
    azimuth['data'] = az_data

    # elevation
    elevation = filemetadata('elevation')
    el0 = sigmet_metadata[first_data_type]['elevation_0']
    el1 = sigmet_metadata[first_data_type]['elevation_1']
    elevation['data'] = mean_of_two_angles_deg(el0, el1).astype('float32')

    # instrument_parameters
    prt = filemetadata('prt')
    prt_mode = filemetadata('prt_mode')
    nyquist_velocity = filemetadata('nyquist_velocity')
    unambiguous_range = filemetadata('unambiguous_range')
    beam_width_h = filemetadata('radar_beam_width_h')
    beam_width_v = filemetadata('radar_beam_width_v')
    pulse_width = filemetadata('pulse_width')
    prt_ratio = filemetadata('prt_ratio')

    prt_value = 1. / sigmetfile.product_hdr['product_end']['prf']
    prt['data'] = prt_value * np.ones(total_rays, dtype='float32')

    ur_value = SPEED_OF_LIGHT * prt_value / 2.
    unambiguous_range['data'] = ur_value * np.ones(total_rays, dtype='float32')

    multi_prf_flag = task_config['task_dsp_info']['multi_prf_flag']
    prf_multiplier = [1, 2, 3, 4][multi_prf_flag]
    if prf_multiplier != 1:
        prt_mode['data'] = np.array(nsweeps * ['dual'], dtype='S')
        ratio = (prf_multiplier + 1) / (prf_multiplier)  # N+1/N
        prt_ratio['data'] = ratio * np.ones(total_rays, dtype='float32')
    else:
        prt_mode['data'] = np.array(nsweeps * ['fixed'], dtype='S')
        prt_ratio['data'] = np.ones(total_rays, dtype='float32')

    wavelength_cm = sigmetfile.product_hdr['product_end']['wavelength']
    nv_value = wavelength_cm / (10000.0 * 4.0 * prt_value) * prf_multiplier
    nyquist_velocity['data'] = nv_value * np.ones(total_rays, dtype='float32')
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
                             'prt_mode': prt_mode,
                             'prt': prt,
                             'prt_ratio': prt_ratio,
                             'nyquist_velocity': nyquist_velocity,
                             'radar_beam_width_h': beam_width_h,
                             'radar_beam_width_v': beam_width_v,
                             'pulse_width': pulse_width}
    if prf_multiplier != 1:
        prf_flag = filemetadata('prf_flag')
        prf_flag['data'] = sigmet_metadata[first_data_type]['prf_flag']
        instrument_parameters['prf_flag'] = prf_flag

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
        sweep_mode['data'] = np.array(
            nsweeps * ['azimuth_surveillance'], dtype='S')
    else:
        sweep_mode['data'] = np.array(nsweeps * ['rhi'], dtype='S')

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters,
        **extended_header_params)


def _is_time_ordered_by_reversal(data, metadata, rays_per_sweep):
    """
    Returns if volume can be time ordered by reversing some or all sweeps.
    True if the volume can be time ordered, False if not.
    """
    if 'XHDR' in data:
        ref_time = data['XHDR'].copy()
        ref_time.shape = ref_time.shape[:-1]
    else:
        ref_time = metadata[metadata.keys()[0]]['time'].astype('int32')
    start = 0
    for nrays in rays_per_sweep:
        if nrays == 0 or nrays == 1:
            continue    # Do not attempt to order sweeps with no rays
        s = slice(start, start + nrays)     # slice which selects sweep
        start += nrays
        sweep_time_diff = np.diff(ref_time[s])
        if np.all(sweep_time_diff >= 0) or np.all(sweep_time_diff <= 0):
            continue
        else:
            return False
    return True


def _is_time_ordered_by_roll(data, metadata, rays_per_sweep):
    """
    Returns if volume can be time ordered by rolling some or all sweeps.
    True if the volume can be time ordered, False if not.
    """
    if 'XHDR' in data:
        ref_time = data['XHDR'].copy()
        ref_time.shape = ref_time.shape[:-1]
    else:
        ref_time = metadata[metadata.keys()[0]]['time'].astype('int32')
    start = 0
    for nrays in rays_per_sweep:
        if nrays == 0 or nrays == 1:
            continue    # Do not attempt to order sweeps with no rays
        s = slice(start, start + nrays)     # slice which selects sweep
        first = ref_time[start]
        last = ref_time[start + nrays - 1]
        start += nrays
        sweep_time_diff = np.diff(ref_time[s])
        count = np.count_nonzero(sweep_time_diff < 0)
        # compare the first and last times for continuity
        if (first - last) < 0:
            count += 1
        if count != 0 and count != 1:
            return False
    return True


def _is_time_ordered_by_reverse_roll(data, metadata, rays_per_sweep):
    """
    Returns if volume can be time ordered by reversing and rolling some or all
    sweeps. True if the volume can be time ordered, False if not.
    """
    if 'XHDR' in data:
        ref_time = data['XHDR'].copy()
        ref_time.shape = ref_time.shape[:-1]
    else:
        ref_time = metadata[metadata.keys()[0]]['time'].astype('int32')
    start = 0
    for nrays in rays_per_sweep:
        if nrays == 0 or nrays == 1:
            continue    # Do not attempt to order sweeps with no rays
        s = slice(start, start + nrays)     # slice which selects sweep
        first = ref_time[start]
        last = ref_time[start + nrays - 1]
        start += nrays
        sweep_time_diff = np.diff(ref_time[s])
        if sweep_time_diff.min() < 0:   # optional reverse
            sweep_time_diff = np.diff(ref_time[s][::-1])
            first, last = last, first
        count = np.count_nonzero(sweep_time_diff < 0)
        # compare the first and last times for continuity
        if (first - last) < 0:
            count += 1
        if count != 0 and count != 1:
            return False
    return True


def _time_order_data_and_metadata_roll(data, metadata, rays_per_sweep):
    """
    Put Sigmet data and metadata in time increasing order using a roll
    operation.
    """
    # Sigmet data is stored by sweep in azimuth or elevation increasing order.
    # Time ordering PPI scans can typically be achieved by rolling the
    # ray collected first to the beginning of the sweep, which is performed
    # here.  Perfect time ordering is achieved if the rays within the sweep
    # were collected sequentially in a clockwise manner from 0 to 360 degrees
    # regardless of the first azimuth collected.
    if 'XHDR' in data:
        ref_time = data['XHDR'].copy()
        ref_time.shape = ref_time.shape[:-1]
    else:
        ref_time = metadata[metadata.keys()[0]]['time'].astype('int32')

    start = 0
    for nrays in rays_per_sweep:
        if nrays == 0 or nrays == 1:
            continue    # Do not attempt to order sweeps with no rays

        s = slice(start, start + nrays)     # slice which selects sweep
        start += nrays
        # determine the number of place by which elements should be shifted.
        sweep_time = ref_time[s]
        sweep_time_diff = np.diff(sweep_time)
        if sweep_time_diff.min() >= 0:
            continue    # already time ordered
        shift = -(sweep_time_diff.argmin() + 1)

        # roll the data and metadata for each field
        for field in data.keys():
            data[field][s] = np.roll(data[field][s], shift, axis=0)
            field_metadata = metadata[field]
            for key in field_metadata.keys():
                field_metadata[key][s] = np.roll(field_metadata[key][s], shift)
    return


def _time_order_data_and_metadata_reverse(data, metadata, rays_per_sweep):
    """
    Put Sigmet data and metadata in time increasing order by reverse sweep in
    time reversed order.
    """
    # Sigmet data is stored by sweep in azimuth or elevation increasing order.
    # Time ordering RHI scans can typically be achieved by reversing the
    # ray order of sweep collected in time from 180 to 0 degrees.
    # Perfect time ordering is achieved if the rays within all sweeps
    # were collected sequentially from 0 to 180 degree or 180 to 0 degrees.
    if 'XHDR' in data:
        ref_time = data['XHDR'].copy()
        ref_time.shape = ref_time.shape[:-1]
    else:
        ref_time = metadata[metadata.keys()[0]]['time'].astype('int32')

    start = 0
    for nrays in rays_per_sweep:
        if nrays == 0 or nrays == 1:
            continue    # Do not attempt to order sweeps with to few rays

        s = slice(start, start + nrays)     # slice which selects sweep
        start += nrays
        # determine the number of place by which elements should be shifted.
        sweep_time = ref_time[s]
        sweep_time_diff = np.diff(sweep_time)
        if sweep_time_diff.min() >= 0:
            continue    # already time ordered, no reversal needed
        # reverse the data and metadata for each field
        for field in data.keys():
            data[field][s] = data[field][s][::-1]
            field_metadata = metadata[field]
            for key in field_metadata.keys():
                field_metadata[key][s] = field_metadata[key][s][::-1]
    return


def _time_order_data_and_metadata_full(data, metadata, rays_per_sweep):
    """
    Put Sigmet data and metadata in time increasing order by sorting the
    times.
    """
    # Sigmet data is stored by sweep in azimuth or elevation increasing order.
    # When rays within the sweeps are collected non-sequentially or in a
    # complex manner, perfect time ordering can only be achived by sorting
    # the times themselves and reordering the rays according to this sort.
    # This ordering method should only be used as a last resort when perfect
    # time ordering is required in the output and other ordering operations
    # (roll, reverse, reverse-roll) will not order the rays correctly.
    if 'XHDR' in data:
        ref_time = data['XHDR'].copy()
        ref_time.shape = ref_time.shape[:-1]
    else:
        ref_time = metadata[metadata.keys()[0]]['time'].astype('int32')

    start = 0
    for nrays in rays_per_sweep:
        if nrays == 0 or nrays == 1:
            continue    # Do not attempt to order sweeps with no rays

        s = slice(start, start + nrays)     # slice which selects sweep
        start += nrays
        # determine the indices which sort the sweep time using a stable
        # sorting algorithm to prevent excessive azimuth scrambling.
        sweep_time = ref_time[s]
        sweep_time_diff = np.diff(ref_time[s])
        if sweep_time_diff.min() >= 0:
            continue    # already time ordered
        sort_idx = np.argsort(sweep_time, kind='mergesort')

        # sort the data and metadata for each field
        for field in data.keys():
            data[field][s] = data[field][s][sort_idx]
            field_metadata = metadata[field]
            for key in field_metadata.keys():
                field_metadata[key][s] = field_metadata[key][s][sort_idx]
    return


def ymds_time_to_datetime(ymds):
    """ Return a datetime object from a Sigmet ymds_time dictionary. """
    dt = datetime.datetime(ymds['year'], ymds['month'], ymds['day'])
    # lowest 10 bits of millisecond parameter specifies milliseconds
    microsec = 1e3 * (ymds['milliseconds'] & 0b1111111111)
    delta = datetime.timedelta(seconds=ymds['seconds'], microseconds=microsec)
    return dt + delta
