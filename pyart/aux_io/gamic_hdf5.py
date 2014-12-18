"""
pyart.aux_io.read_gamic
=======================

Utilities for reading gamic hdf5 files.

.. autosummary::
    :toctree: generated/

    read_gamic
    _h5_to_dict
    _h5_moments_to_dict

"""

import datetime
import warnings

import h5py
import numpy as np

from ..config import FileMetadata, get_fillvalue
from ..io.common import make_time_unit_str
from ..core.radar import Radar


LIGHT_SPEED = 2.99792458e8  # speed of light in meters per second


def read_gamic(filename, field_names=None, additional_metadata=None,
               file_field_names=False, exclude_fields=None,
               valid_range_from_file=True, units_from_file=True):
    """
    Read a GAMIC hdf5 file.

    Parameters
    ----------
    filename : str
        Name of GAMIC HDF5 file to read data from.
    field_names : dict, optional
        Dictionary mapping field names in the file names to radar field names.
        Unlike other read functions, fields not in this dictionary or having a
        value of None are still included in the radar.fields dictionary, to
        exclude them use the `exclude_fields` parameter. Fields which are
        mapped by this dictionary will be renamed from key to value.
    additional_metadata : dict of dicts, optional
        This parameter is not used, it is included for uniformity.
    file_field_names : bool, optional
        True to force the use of the field names from the file in which
        case the `field_names` parameter is ignored. False will use to
        `field_names` parameter to rename fields.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters.
    valid_range_from_file : bool, optional
        True to extract valid range (valid_min and valid_max) for all
        field from the file when they are present.  False will not extract
        these parameters.
    units_from_file : bool, optional
        True to extract the units for all fields from the file when available.
        False will not extract units using the default units for the fields.

    Returns
    -------
    radar : Radar
        Radar object.

    """
    # TODO
    # * refactor

    # * unit tests
    # * auto-detect file type with pyart.io.read function
    # * move to pyart.io namespace

    # create metadata retrieval object
    filemetadata = FileMetadata('gamic', field_names, additional_metadata,
                                file_field_names, exclude_fields)

    # Open HDF5 file and get handle
    hfile = h5py.File(filename, 'r')

    # metadata, initial empty, filled in throughout the function
    metadata = filemetadata('metadata')

    # /where
    # latitude, longitude and altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')

    h_where = hfile['where'].attrs
    latitude['data'] = np.array([h_where['lat']], dtype='float64')
    longitude['data'] = np.array([h_where['lon']], dtype='float64')
    altitude['data'] = np.array([h_where['height']], dtype='float64')

    # /what
    what_attrs = hfile['what'].attrs
    nsweeps = int(what_attrs['sets'])
    # fill in metadata
    metadata['original_container'] = 'GAMIC-HDF5'
    what_mapping = {
        'version': 'gamic_version', 'sets_scheduled': 'sets_scheduled',
        'object': 'gamic_object', 'date': 'gamic_date'}
    for hkey, mkey in what_mapping.items():
        if hkey in what_attrs:
            metadata[mkey] = what_attrs[hkey]

    # /how
    # fill in metadata
    how_attrs = hfile['how'].attrs
    how_keys = ['software', 'template_name', 'site_name', 'host_name',
                'azimuth_beam', 'elevation_beam', 'sdp_name', 'sw_version',
                'sdp_version', 'simulated']
    for key in how_keys:
        if key in how_attrs:
            metadata[key] = how_attrs[key]

    # /scan0 to /scan[nsweeps-1]
    scans = ['scan%i' % (i) for i in range(nsweeps)]
    # verify that scan0 ... scan[nsweeps-1] are present in file
    for scan in scans:
        assert scan in hfile

    # sweep_start_ray_index, sweep_end_ray_index
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')

    rays_per_sweep = [hfile[s]['how'].attrs['ray_count'] for s in scans]
    total_rays = sum(rays_per_sweep)
    ssri = np.cumsum(np.append([0], rays_per_sweep[:-1])).astype('int32')
    seri = np.cumsum(rays_per_sweep).astype('int32') - 1
    sweep_start_ray_index['data'] = ssri
    sweep_end_ray_index['data'] = seri

    # sweep number
    sweep_number = filemetadata('sweep_number')
    sweep_number['data'] = np.array(
        [hfile[scan]['what'].attrs['set_idx'] for scan in scans],
        dtype='int32')

    # sweep_type
    scan_type = hfile['scan0']['what'].attrs['scan_type'].lower()
    # check that all scans in the volume are the same type
    for scan in scans:
        if hfile[scan]['what'].attrs['scan_type'].lower() != scan_type:
            raise NotImplementedError('Mixed scan_type volume in scan:' + scan)
    if scan_type not in ['ppi', 'rhi']:
        message = "Unknown scan type: %s, reading as RHI scans." % (scan_type)
        warnings.warn(message)
        scan_type = 'rhi'

    # sweep_mode, fixed_angle
    sweep_mode = filemetadata('sweep_mode')
    fixed_angle = filemetadata('fixed_angle')

    if scan_type == 'rhi':
        sweep_mode['data'] = np.array(nsweeps * ['rhi'])
        sweep_az = [hfile[s]['how'].attrs['azimuth'] for s in scans]
        fixed_angle['data'] = np.array(sweep_az, dtype='float32')
    elif scan_type == 'ppi':
        sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])
        sweep_elev = [hfile[s]['how'].attrs['elevation'] for s in scans]
        fixed_angle['data'] = np.array(sweep_elev, dtype='float32')

    # time
    time = filemetadata('time')
    t_data = np.empty((total_rays, ), dtype='int64')
    for scan, start, stop in zip(scans, ssri, seri):
        # retrieve the time for each ray from the header timestamp (usec)
        t_data[start:stop+1] = hfile[scan]['ray_header']['timestamp']
    start_epoch = t_data[0] // 1.e6     # truncate to second resolution
    start_time = datetime.datetime.utcfromtimestamp(start_epoch)
    time['units'] = make_time_unit_str(start_time)
    time['data'] = ((t_data - start_epoch * 1.e6) / 1.e6).astype('float64')

    # range
    _range = filemetadata('range')
    ngates = int(hfile['/scan0/how'].attrs['bin_count'])
    range_start = float(hfile['/scan0/how'].attrs['range_start'])
    range_step = float(hfile['/scan0/how'].attrs['range_step'])
    # range_step may need to be scaled by range_samples
    # XXX This gives distances to start of gates not center, this matches
    # Radx but may be incorrect, add range_step / 2. for center
    _range['data'] = (np.arange(ngates, dtype='float32') * range_step +
                      range_start)
    _range['meters_to_center_of_first_gate'] = range_start
    _range['meters_between_gates'] = range_step

    # elevation
    elevation = filemetadata('elevation')
    el_data = np.empty((total_rays, ), dtype='float32')
    for scan, start, stop in zip(scans, ssri, seri):
        startel = hfile[scan]['ray_header']['elevation_start']
        stopel = hfile[scan]['ray_header']['elevation_stop']
        scan_el = np.angle(
            (np.exp(1.j*np.deg2rad(startel)) +
             np.exp(1.j*np.deg2rad(stopel))) / 2., deg=True)
        el_data[start:stop+1] = scan_el
    elevation['data'] = el_data

    # azimuth
    azimuth = filemetadata('azimuth')
    az_data = np.empty((total_rays, ), dtype='float32')
    for scan, start, stop in zip(scans, ssri, seri):
        startaz = hfile[scan]['ray_header']['azimuth_start']
        stopaz = hfile[scan]['ray_header']['azimuth_stop']
        scan_az = np.angle(
            (np.exp(1.j*np.deg2rad(startaz)) +
             np.exp(1.j*np.deg2rad(stopaz))) / 2., deg=True) % 360.
        az_data[start:stop+1] = scan_az
    azimuth['data'] = az_data

    # fields
    fields = {}
    h_keys = [k for k in hfile['/scan0'] if k.startswith('moment_')]
    h_fields = [hfile['/scan0'][k].attrs['moment'] for k in h_keys]

    for h_field, h_key in zip(h_fields, h_keys):
        field_name = filemetadata.get_field_name(h_field)
        if field_name is None:
            continue
        fdata = np.ma.zeros((total_rays, ngates), dtype='float32')
        # loop over the sweeps, copy data into fdata
        for scan, start, stop in zip(scans, ssri, seri):
            fdata[start:stop+1] = _get_gamic_sweep_data(hfile[scan][h_key])[:]
        field_dic = filemetadata(field_name)
        field_dic['data'] = fdata
        field_dic['_FillValue'] = get_fillvalue()
        if valid_range_from_file:
            try:
                valid_min = hfile['/scan0'][h_key].attrs['dyn_range_min']
                valid_max = hfile['/scan0'][h_key].attrs['dyn_range_max']
                field_dic['valid_min'] = valid_min
                field_dic['valid_max'] = valid_max
            except:
                pass
        if units_from_file:
            try:
                field_dic['units'] = hfile['/scan0'][h_key].attrs['unit']
            except:
                pass
        fields[field_name] = field_dic

    # ray_angle_res
    ray_angle_res = filemetadata('ray_angle_res')
    ray_angle_res['data'] = np.array(
        [hfile[scan]['how'].attrs['angle_step'] for scan in scans],
        dtype='float32')

    # rays_are_indexed
    rays_are_indexed = filemetadata('rays_are_indexed')
    rays_are_indexed['data'] = np.array(
        [['false', 'true'][hfile[scan]['how'].attrs['angle_sync']] for
         scan in scans])

    # target_scan_rate
    target_scan_rate = filemetadata('target_scan_rate')
    target_scan_rate['data'] = np.array(
        [hfile[scan]['how'].attrs['scan_speed'] for scan in scans],
        dtype='float32')

    # scan_rate
    scan_rate = filemetadata('scan_rate')
    if scan_type == 'ppi':
        sweep_rates = [hfile[s]['ray_header']['az_speed'] for s in scans]
        scan_rate['data'] = np.concatenate(sweep_rates).astype('float32')
    elif scan_type == 'rhi':
        sweep_rates = [hfile[s]['ray_header']['el_speed'] for s in scans]
        scan_rate['data'] = np.concatenate(sweep_rates).astype('float32')
    else:
        scan_rate = None

    # instrument_parameters
    frequency = filemetadata('frequency')
    frequency['data'] = np.array(
        [LIGHT_SPEED / hfile['scan0']['how'].attrs['radar_wave_length']],
        dtype='float32')

    beamwidth_h = filemetadata('radar_beam_width_h')
    beamwidth_h['data'] = np.array(
        [hfile['how'].attrs['azimuth_beam']], dtype='float32')

    beamwidth_v = filemetadata('radar_beam_width_v')
    beamwidth_v['data'] = np.array(
        [hfile['how'].attrs['elevation_beam']], dtype='float32')

    pulse_width = filemetadata('pulse_width')
    t = [hfile[s]['how'].attrs['pulse_width_us']*1e-6 for s in scans]
    pulse_width['data'] = np.repeat(t, rays_per_sweep).astype('float32')

    prt = filemetadata('prt')
    t = [1. / hfile[s]['how'].attrs['PRF'] for s in scans]
    prt['data'] = np.repeat(t, rays_per_sweep).astype('float32')

    prt_mode = filemetadata('prt_mode')
    prt_ratio = filemetadata('prt_ratio')
    sweep_unfolding = [hfile[s]['how'].attrs['unfolding'] for s in scans]
    sweep_prt_mode = map(_prt_mode_from_unfolding, sweep_unfolding)
    sweep_prt_ratio = [[1, 2./3., 3./4., 4./5.][i] for i in sweep_unfolding]
    prt_mode['data'] = np.array(sweep_prt_mode)
    prt_ratio['data'] = np.repeat(sweep_prt_ratio, rays_per_sweep)

    unambiguous_range = filemetadata('unambiguous_range')
    t = [hfile[s]['how'].attrs['range'] for s in scans]
    unambiguous_range['data'] = np.repeat(t, rays_per_sweep).astype('float32')

    nyquist_velocity = filemetadata('nyquist_velocity')
    t = [float(hfile[s]['how']['extended'].attrs['nyquist_velocity'])
         for s in scans]
    nyquist_velocity['data'] = np.repeat(t, rays_per_sweep).astype('float32')

    n_samples = filemetadata('n_samples')
    t = [hfile[s]['how'].attrs['range_samples'] *
         hfile[s]['how'].attrs['time_samples'] for s in scans]
    n_samples['data'] = np.repeat(t, rays_per_sweep).astype('int32')

    instrument_parameters = {
        'frequency': frequency,
        'radar_beam_width_h': beamwidth_h,
        'radar_beam_width_v': beamwidth_v,
        'n_samples': n_samples,
        'pulse_width': pulse_width,
        'prt': prt,
        'prt_ratio': prt_ratio,
        'nyquist_velocity': nyquist_velocity,
        'unambiguous_range': unambiguous_range,
        'prt_mode': prt_mode,
    }

    hfile.close()

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters,
        ray_angle_res=ray_angle_res,
        rays_are_indexed=rays_are_indexed,
        scan_rate=scan_rate, target_scan_rate=target_scan_rate)


def _prt_mode_from_unfolding(unfolding):
    """ Return 'fixed' or 'staggered' depending on unfolding flag """
    if unfolding == 0:
        return 'fixed'
    else:
        return 'staggered'


def _get_gamic_sweep_data(group):
    """ Get GAMIC HDF5 sweep data from an HDF5 group. """
    dyn_range_min = group.attrs['dyn_range_min']
    dyn_range_max = group.attrs['dyn_range_max']
    raw_data = group[:]
    fmt = group.attrs['format']
    if fmt == 'UV16':
        # unsigned 16-bit integer data, 0 indicates a masked value
        assert raw_data.dtype == np.uint16
        scale = (dyn_range_max - dyn_range_min) / 65535.
        offset = dyn_range_min
        sweep_data = np.ma.masked_array(
            raw_data * scale + offset, mask=(raw_data == 0), dtype='float32')
    elif fmt == 'UV8':
        # unsigned 8-bit integer data, 0 indicates a masked value
        assert raw_data.dtype == np.uint8
        scale = (dyn_range_max - dyn_range_min) / 255.
        offset = dyn_range_min
        sweep_data = np.ma.masked_array(
            raw_data * scale + offset, mask=(raw_data == 0), dtype='float32')
    else:
        raise NotImplementedError('GAMIC data format: %s', fmt)
    return sweep_data
