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

import h5py
import numpy as np

from ..config import FileMetadata, get_fillvalue
from ..io.common import make_time_unit_str
from ..core.radar import Radar


def read_gamic(filename, field_names=None, additional_metadata=None,
               file_field_names=False, exclude_fields=None):
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

    Returns
    -------
    radar : Radar
        Radar object.

    """
    # TODO
    # * unit tests
    # * add default field mapping, etc to default config
    # * auto-detect file type with pyart.io.read function
    # * move to pyart.io namespace

    # create metadata retrieval object
    filemetadata = FileMetadata('gamic', field_names, additional_metadata,
                                file_field_names, exclude_fields)

    # Open HDF5 file and get handle
    hfile = h5py.File(filename, 'r')

    # metadata, initial empty is filled in throughout
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
    if scan_type not in ['ppi', 'rhi']:
        raise NotImplementedError('Not supported sweep_type: ' + scan_type)
    # check that all scans in the volume are the same type
    for scan in scans:
        if hfile[scan]['what'].attrs['scan_type'].lower() != scan_type:
            raise NotImplementedError('Mixed scan_type volume in scan:' + scan)

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
        fields[field_name] = field_dic

    # instrument_parameters
    instrument_parameters = None

    # calibration
    radar_calibration = None

    hfile.close()

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters,
        radar_calibration=radar_calibration)


def _get_gamic_sweep_data(group):
    """ Get GAMIC HDF5 sweep data from a HDF5 group. """
    return group[:]
