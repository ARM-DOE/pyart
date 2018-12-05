"""
pyart.aux_io.sinarame_h5
========================

Routines for reading sinarame_H5 files.

.. autosummary::
    :toctree: generated/

    read_sinarame_h5
    write_sinarame_cfradial
    _to_str
    _get_SINARAME_h5_sweep_data

"""

from __future__ import print_function

from datetime import datetime
import glob
import os

try:
    from netcdftime import utime
except ImportError:
    from cftime import utime

import numpy as np

try:
    import h5py
    _H5PY_AVAILABLE = True
except ImportError:
    _H5PY_AVAILABLE = False

from ..config import FileMetadata, get_fillvalue
from ..io.common import make_time_unit_str, _test_arguments
from ..io import write_cfradial
from ..core.radar import Radar
from ..exceptions import MissingOptionalDependency


SINARAME_H5_FIELD_NAMES = {
    'TH': 'total_power',        # uncorrected reflectivity, horizontal
    'TV': 'total_power_vertical',        # uncorrected reflectivity, vertical
    'DBZH': 'reflectivity',     # corrected reflectivity, horizontal
    'DBZV': 'reflectivity',     # corrected reflectivity, vertical
    'CM': 'clutter mask',
    'ZDR': 'differential_reflectivity',     # differential reflectivity
    'RHOHV': 'cross_correlation_ratio',
    'LDR': 'linear_polarization_ratio',
    'PHIDP': 'differential_phase',
    'KDP': 'specific_differential_phase',
    'SQI': 'normalized_coherent_power',
    'SNR': 'signal_to_noise_ratio',
    'VRAD': 'velocity',
    'WRAD': 'spectrum_width',
    'QIND': 'quality_index',
}


def read_sinarame_h5(filename, field_names=None, additional_metadata=None,
                     file_field_names=False, exclude_fields=None,
                     include_fields=None, **kwargs):
    """
    Read a SINARAME_H5 file.

    Parameters
    ----------
    filename : str
        Name of the SINARAME_H5 file to read.
    field_names : dict, optional
        Dictionary mapping SINARAME_H5 field names to radar field names. If a
        data type found in the file does not appear in this dictionary or has
        a value of None it will not be placed in the radar.fields dictionary.
        A value of None, the default, will use the mapping defined in the
        Py-ART configuration file.
    additional_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve metadata from during this read.
        This metadata is not used during any successive file reads unless
        explicitly included.  A value of None, the default, will not
        introduct any addition metadata and the file specific or default
        metadata as specified by the Py-ART configuration file will be used.
    file_field_names : bool, optional
        True to use the MDV data type names for the field names. If this
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


    Returns
    -------
    radar : Radar
        Radar object containing data from SINARAME_H5 file.

    """
    # TODO before moving to pyart.io
    # * unit test
    # * add default field mapping, etc to default config
    # * auto-detect file type with pyart.io.read function
    # * instrument parameters
    # * add additional checks for HOW attributes
    # * support for other objects (SCAN, XSEC)

    # check that h5py is available
    if not _H5PY_AVAILABLE:
        raise MissingOptionalDependency(
            "h5py is required to use read_sinarame_h5 but is not installed")

    # test for non empty kwargs
    _test_arguments(kwargs)

    # create metadata retrieval object
    if field_names is None:
        field_names = SINARAME_H5_FIELD_NAMES
    filemetadata = FileMetadata('SINARAME_h5', field_names,
                                additional_metadata,
                                file_field_names, exclude_fields,
                                include_fields)

    # open the file
    hfile = h5py.File(filename, 'r')
    SINARAME_object = _to_str(hfile['what'].attrs['object'])
    if SINARAME_object not in ['PVOL', 'SCAN', 'ELEV', 'AZIM']:
        raise NotImplementedError(
            'object: %s not implemented.' % (SINARAME_object))

    # determine the number of sweeps by the number of groups which
    # begin with dataset
    datasets = [k for k in hfile if k.startswith('dataset')]
    datasets.sort(key=lambda x: int(x[7:]))
    nsweeps = len(datasets)

    # latitude, longitude and altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')

    h_where = hfile['where'].attrs
    latitude['data'] = np.array([h_where['lat']], dtype='float64')
    longitude['data'] = np.array([h_where['lon']], dtype='float64')
    altitude['data'] = np.array([h_where['height']], dtype='float64')

    # metadata
    metadata = filemetadata('metadata')
    metadata['source'] = _to_str(hfile['what'].attrs['source'])
    metadata['original_container'] = 'SINARAME_h5'
    metadata['SINARAME_conventions'] = _to_str(hfile.attrs['Conventions'])

    h_what = hfile['what'].attrs
    metadata['version'] = _to_str(h_what['version'])
    metadata['source'] = _to_str(h_what['source'])

    try:
        ds1_how = hfile[datasets[0]]['how'].attrs
    except KeyError:
        # if no how group exists mock it with an empty dictionary
        ds1_how = {}
    if 'system' in ds1_how:
        metadata['system'] = ds1_how['system']
    if 'software' in ds1_how:
        metadata['software'] = ds1_how['software']
    if 'sw_version' in ds1_how:
        metadata['sw_version'] = ds1_how['sw_version']

    # sweep_start_ray_index, sweep_end_ray_index
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')

    if SINARAME_object in ['AZIM', 'SCAN', 'PVOL']:
        rays_per_sweep = [
            int(hfile[d]['where'].attrs['nrays']) for d in datasets]
    elif SINARAME_object == 'ELEV':
        rays_per_sweep = [
            int(hfile[d]['where'].attrs['angles'].size) for d in datasets]
    total_rays = sum(rays_per_sweep)
    ssri = np.cumsum(np.append([0], rays_per_sweep[:-1])).astype('int32')
    seri = np.cumsum(rays_per_sweep).astype('int32') - 1
    sweep_start_ray_index['data'] = ssri
    sweep_end_ray_index['data'] = seri

    # sweep_number
    sweep_number = filemetadata('sweep_number')
    sweep_number['data'] = np.arange(nsweeps, dtype='int32')

    # sweep_mode
    sweep_mode = filemetadata('sweep_mode')
    sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])

    # scan_type
    if SINARAME_object == 'ELEV':
        scan_type = 'rhi'
    else:
        scan_type = 'ppi'

    # fixed_angle
    fixed_angle = filemetadata('fixed_angle')
    if SINARAME_object == 'ELEV':
        sweep_el = [hfile[d]['where'].attrs['az_angle'] for d in datasets]
    else:
        sweep_el = [hfile[d]['where'].attrs['elangle'] for d in datasets]
    fixed_angle['data'] = np.array(sweep_el, dtype='float32')

    # elevation
    elevation = filemetadata('elevation')
    if 'elangles' in ds1_how:
        edata = np.empty(total_rays, dtype='float32')
        for d, start, stop in zip(datasets, ssri, seri):
            edata[start:stop+1] = hfile[d]['how'].attrs['elangles'][:]
        elevation['data'] = edata
    elif SINARAME_object == 'ELEV':
        edata = np.empty(total_rays, dtype='float32')
        for d, start, stop in zip(datasets, ssri, seri):
            edata[start:stop+1] = hfile[d]['where'].attrs['angles'][:]
        elevation['data'] = edata
    else:
        elevation['data'] = np.repeat(sweep_el, rays_per_sweep)

    # range
    _range = filemetadata('range')
    if 'rstart' in hfile['dataset1/where'].attrs:
        # derive range from rstart and rscale attributes if available

        # check that the gate spacing is constant between sweeps
        rstart = [hfile[d]['where'].attrs['rstart'] for d in datasets]
        if any(rstart != rstart[0]):
            raise ValueError('range start changes between sweeps')
        rscale = [hfile[d]['where'].attrs['rscale'] for d in datasets]
        if any(rscale != rscale[0]):
            raise ValueError('range scale changes between sweeps')
        nbins = int(hfile['dataset1']['where'].attrs['nbins'])
        _range['data'] = (np.arange(nbins, dtype='float32') * rscale[0] +
                          rstart[0])
        _range['meters_to_center_of_first_gate'] = rstart[0]
        _range['meters_between_gates'] = float(rscale[0])
    else:
        # if not defined use range attribute which defines the maximum range in
        # km.  There is no information on the starting location of the range
        # bins so we assume this to be 0.
        # This most often occurs in RHI files, which technically do not meet
        # the ODIM 2.2 specs. Section 7.4 requires that these files include
        # the where/rstart, where/rscale and where/nbins attributes.
        max_range = [hfile[d]['where'].attrs['range'] for d in datasets]
        if any(max_range != max_range[0]):
            raise ValueError('maximum range changes between sweeps')
        # nbins is required
        nbins = hfile['dataset1/data1/data'].shape[1]
        _range['data'] = np.linspace(
            0, max_range[0] * 1000., nbins, dtype='float32')
        _range['meters_to_center_of_first_gate'] = 0
        _range['meters_between_gates'] = max_range[0] * 1000. / nbins

    # azimuth
    azimuth = filemetadata('azimuth')
    az_data = np.ones((total_rays, ), dtype='float32')
    for dset, start, stop in zip(datasets, ssri, seri):
        if SINARAME_object == 'ELEV':
            # all azimuth angles are the sweep azimuth angle
            sweep_az = hfile[dset]['where'].attrs['az_angle']
        elif SINARAME_object == 'AZIM':
            # Sector azimuths are specified in the startaz and stopaz
            # attribute of dataset/where.
            # Assume that the azimuth angles do not pass through 0/360 degrees
            startaz = hfile[dset]['where'].attrs['startaz']
            stopaz = hfile[dset]['where'].attrs['stopaz']
            nrays = stop - start + 1
            sweep_az = np.linspace(startaz, stopaz, nrays, endpoint=True)
        elif ('startazA' in ds1_how) and ('stopazA' in ds1_how):
            # average between start and stop azimuth angles
            startaz = hfile[dset]['how'].attrs['startazA']
            stopaz = hfile[dset]['how'].attrs['stopazA']
            sweep_az = np.angle(
                (np.exp(1.j*np.deg2rad(startaz)) +
                 np.exp(1.j*np.deg2rad(stopaz))) / 2., deg=True)
        else:
            # according to section 5.1 the first ray points north (0 degrees)
            # and proceeds clockwise for a complete 360 rotation.
            nrays = stop - start + 1
            sweep_az = np.linspace(0, 360, nrays, endpoint=False)
        az_data[start:stop+1] = sweep_az
    azimuth['data'] = az_data

    # time
    _time = filemetadata('time')
    if ('startazT' in ds1_how) and ('stopazT' in ds1_how):
        # average between startazT and stopazT
        t_data = np.empty((total_rays, ), dtype='float32')
        for dset, start, stop in zip(datasets, ssri, seri):
            t_start = hfile[dset]['how'].attrs['startazT']
            t_stop = hfile[dset]['how'].attrs['stopazT']
            t_data[start:stop+1] = (t_start + t_stop) / 2
        start_epoch = t_data.min()
        start_time = datetime.utcfromtimestamp(start_epoch)
        _time['units'] = make_time_unit_str(start_time)
        _time['data'] = t_data - start_epoch
    else:
        t_data = np.empty((total_rays, ), dtype='int32')
        # interpolate between each sweep starting and ending time
        for dset, start, stop in zip(datasets, ssri, seri):
            dset_what = hfile[dset]['what'].attrs
            start_str = _to_str(
                dset_what['startdate'] + dset_what['starttime'])
            end_str = _to_str(dset_what['enddate'] + dset_what['endtime'])
            start_dt = datetime.strptime(start_str, '%Y%m%d%H%M%S')
            end_dt = datetime.strptime(end_str, '%Y%m%d%H%M%S')

            time_delta = end_dt - start_dt
            delta_seconds = time_delta.seconds + time_delta.days * 3600 * 24
            rays = stop - start + 1
            sweep_start_epoch = (
                start_dt - datetime(1970, 1, 1)).total_seconds()
            t_data[start:stop+1] = (sweep_start_epoch +
                                    np.linspace(0, delta_seconds, rays))
        start_epoch = t_data.min()
        start_time = datetime.utcfromtimestamp(start_epoch)
        _time['units'] = make_time_unit_str(start_time)
        _time['data'] = (t_data - start_epoch).astype('float32')

    # fields
    fields = {}
    h_field_keys = [k for k in hfile['dataset1'] if k.startswith('data')]
    SINARAME_fields = [hfile['dataset1'][d]['what'].attrs['quantity'] for d in
                       h_field_keys]
    for SINARAME_field, h_field_key in zip(SINARAME_fields, h_field_keys):
        field_name = filemetadata.get_field_name(_to_str(SINARAME_field))
        if field_name is None:
            continue
        fdata = np.ma.zeros((total_rays, nbins), dtype='float32')
        start = 0
        # loop over the sweeps, copy data into correct location in data array
        for dset, rays_in_sweep in zip(datasets, rays_per_sweep):
            sweep_data = _get_SINARAME_h5_sweep_data(hfile[dset][h_field_key])
            sweep_nbins = sweep_data.shape[1]
            fdata[start:start + rays_in_sweep, :sweep_nbins] = sweep_data[:]
            start += rays_in_sweep
        # create field dictionary
        field_dic = filemetadata(field_name)
        field_dic['data'] = fdata
        field_dic['_FillValue'] = get_fillvalue()
        fields[field_name] = field_dic
        # Add missing metadata
        if file_field_names:
            fields[field_name].update(
                filemetadata.get_metadata(field_names[field_name]))

    # instrument_parameters
    instrument_parameters = None

    return Radar(
        _time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)


def write_sinarame_cfradial(path):
    """
    This function takes SINARAME_H5 files (where every file has only one field
    and one volume) from a folder and writes a CfRadial file for each volume
    including all fields.

    Parameters
    ----------
    path : str
        Where the SINARAME_H5 files are.

    """

    path_user = os.path.expanduser(path)

    TH_list = glob.glob(path_user + '/*_TH_*.H5')

    file_date = [i.split('_')[-1][9:-4] for i in TH_list]
    file_date.sort()

    for i in file_date:
        files = glob.glob(path_user + '/*' + i + 'Z.H5')

        # I want it to start with TV or TH because of some range issue
        files.sort(reverse=True)
        files.sort(key=lambda x: len(x.split('_')[-2]))

        for j in np.arange(len(files)):
            basename = os.path.basename(files[j])
            bs = basename.split('_')
            base1 = '{b1}_{b2}_{b3}_{fn}_{b4}'.format(
                b1=bs[0], b2=bs[1], b3=bs[2], fn=bs[3], b4=bs[4])
            file = '{path}/{base1}'.format(path=path_user, base1=base1)

            if j == 0:
                try:
                    radar = read_sinarame_h5(file, file_field_names=True)
                    azi_shape, range_shape = radar.fields['TV']['data'].shape

                except ValueError:
                    print("x - this Radar wasn't created",
                          base1, sep='\t')

                # In case the H5 file was badly created with bbufr it
                # throws a KeyError
                except KeyError:
                    print("x - Wrong BUFR conversion", base1, sep='\t')

            else:
                try:
                    radar_tmp = read_sinarame_h5(file, file_field_names=True)

                    field = radar_tmp.fields.keys()[0]

                    # Add missing gates
                    if radar_tmp.fields[field]['data'].shape[1] != range_shape:
                        n_missing_gates = \
                            (range_shape
                             - radar_tmp.fields[field]['data'].shape[1])
                        fill = np.ma.masked_all((azi_shape, n_missing_gates))
                        data = np.ma.concatenate(
                            [radar_tmp.fields[field]['data'], fill], 1)
                        radar_tmp.fields[field]['data'] = data

                    radar.fields.update(radar_tmp.fields)

                # Same as above, bbufr conversion errors
                except KeyError:
                    print("x - Wrong BUFR conversion", base1, sep='\t')

                except ValueError:
                    print("x - this Radar wasn't created",
                          base1, sep='\t')

                # In case the first file didn't create a Radar object
                except NameError:
                    print("Radar didn't exist, creating")
                    radar = read_sinarame_h5(file, file_field_names=True)

        cal_temps = u"gregorian"
        cdftime = utime(radar.time['units'])

        time1 = cdftime.num2date(radar.time['data'][0]).strftime(
            '%Y%m%d_%H%M%S')
        time2 = cdftime.num2date(radar.time['data'][-1]).strftime(
            '%Y%m%d_%H%M%S')

        radar._DeflateLevel = 5

        # cfrad.TIME1_to_TIME2_NAME_VCP_RANGE.nc
        cffile = 'cfrad.{time1}.0000_to_{time2}.0000'\
                 '_{b1}_{est}_{ran}'.format(time1=time1, time2=time2,
                                            b1=bs[0], est=bs[1], ran=bs[2])

        print('Writing to {path}{filename}.nc'.format(path=path_user,
                                                      filename=cffile))
        write_cfradial(path_user + '/' + cffile + '.nc',
                       radar, format='NETCDF4_CLASSIC')


def _to_str(text):
    """ Convert bytes to str if necessary. """
    if hasattr(text, 'decode'):
        return text.decode('utf-8')
    else:
        return text


def _get_SINARAME_h5_sweep_data(group):
    """ Get SINARAME_H5 sweet data from an HDF5 group. """

    # mask raw data
    what = group['what']
    raw_data = group['data'][:]

    if 'nodata' in what.attrs:
        nodata = what.attrs.get('nodata')
        data = np.ma.masked_equal(raw_data, nodata)
    else:
        data = np.ma.masked_array(raw_data)
    if 'undetect' in what.attrs:
        undetect = what.attrs.get('undetect')
        data[data == undetect] = np.ma.masked

    offset = 0.0
    gain = 1.0
    if 'offset' in what.attrs:
        offset = what.attrs.get('offset')
    if 'gain' in what.attrs:
        gain = what.attrs.get('gain')
    return data * gain + offset
