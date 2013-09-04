"""
pyart.io.sigmet
===============

Reading and writing of Sigmet (raw format) files

.. autosummary::
    :toctree: generated/

    read_sigmet
    ymds_time_to_datetime
    _time_order_data_and_metadata

"""

import datetime

import numpy as np

from .common import get_metadata, make_time_unit_str
from .radar import Radar
from ._sigmetfile import SigmetFile, bin4_to_angle, bin2_to_angle

SPEED_OF_LIGHT = 299793000.0


def read_sigmet(filename, field_names=None, field_metadata=None,
                sigmet_field_names=False, time_ordered=True, debug=False):
    """
    Read a Sigmet (IRIS) product file.

    Parameters
    ----------
    filename : str
        Name of Sigmet (IRIS) product file to read.
    field_names : dict, optional
        Dictionary mapping Sigmet data type names to radar field names. If a
        data type found in the file does not appear in this dictionary it will
        not be placed in the radar.fields dictionary.  If None (default) a
        standard mapping will be used (SIGMET_TO_STANDARD)
    field_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve field metadata from, if a mapped
        field does not appear in this list standard metadata will be used, or
        none if these is no metadata.  None (default) uses a blank dictionary.
    sigmet_field_names : bool, optional
        True to use the Sigmet data type names for the field names. If this
        case the field_metadata and field_names parameters are ignored and the
        returned radar object has a fields attribute filled with the sigmet
        data type names with no metadata.
    time_ordered : bool, optional
        True (default) to return data in Radar object in time increasing
        order.  False will return the data ordered as present in the file.
    debug : bool, optional

    Returns
    -------
    radar : Radar
        Radar object

    """
    # parse parameters
    if field_names is None:
        field_names = SIGMET_TO_STANDARD
    if field_metadata is None:
        field_metadata = {}

    # open the file, read data
    sigmetfile = SigmetFile(filename, debug=debug)
    sigmet_data, sigmet_metadata = sigmetfile.read_data()
    ingest_config = sigmetfile.ingest_header['ingest_configuration']
    task_config = sigmetfile.ingest_header['task_configuration']
    first_data_type = sigmetfile.data_type_names[0]
    sigmetfile.close()
    nsweeps, nrays, nbins = sigmet_data[first_data_type].shape

    # mark times as missing on rays not collected
    bad_rays = (sigmet_metadata[first_data_type]['nbins'] == -1)
    for k, v in sigmet_metadata.iteritems():
        sigmet_metadata[k]['time'] = np.ma.masked_array(v['time'], bad_rays)

    # time order
    if time_ordered:
        _time_order_data_and_metadata(sigmet_data, sigmet_metadata)

    # remove missing rays from the data
    good_rays = (sigmet_metadata[first_data_type]['nbins'] != -1)
    for field_name in sigmet_data.keys():
        sigmet_data[field_name] = sigmet_data[field_name][good_rays]

        field_metadata = sigmet_metadata[field_name]
        for key in field_metadata.keys():
            field_metadata[key] = field_metadata[key][good_rays]

    # sweep_start_ray_index and sweep_end_ray_index
    ray_count = good_rays.sum(axis=1)
    sweep_start_ray_index = get_metadata('sweep_start_ray_index')
    sweep_end_ray_index = get_metadata('sweep_end_ray_index')

    ssri = np.cumsum(np.append([0], ray_count[:-1])).astype('int32')
    sweep_start_ray_index['data'] = ssri
    sweep_end_ray_index['data'] = np.cumsum(ray_count).astype('int32') - 1

    # time
    time = get_metadata('time')

    # add sweep_start time to all time values in each sweep.
    dts = [ymds_time_to_datetime(d['sweep_start_time'])
           for d in sigmetfile.ingest_data_headers[first_data_type]]

    tdata = sigmet_metadata[first_data_type]['time'].astype('float64')
    for i, dt in enumerate(dts):
        start = sweep_start_ray_index['data'][i]
        end = sweep_end_ray_index['data'][i]
        tdata[start: end + 1] += (dt - dts[0]).total_seconds()
    time['data'] = tdata.filled()
    time['units'] = make_time_unit_str(dts[0])

    # _range
    _range = get_metadata('range')
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
        if sigmet_field_names:
            fields[data_type_name] = {'data': fdata.reshape(-1, nbins)}
        elif data_type_name in field_names:

            field_name = field_names[data_type_name]

            if field_name in field_metadata:
                field_dic = field_metadata[field_name].copy()
            else:
                field_dic = get_metadata(field_name)

            field_dic['data'] = fdata
            field_dic['_FillValue'] = -9999.0
            fields[field_name] = field_dic

    # metadata
    metadata = {'title': '', 'institution': '', 'references': '',
                'source': '', 'history': '', 'comment': ''}
    metadata['original_container'] = 'sigmet'
    metadata['instrument_name'] = ingest_config['site_name'].strip()

    # scan_type
    if task_config['task_scan_info']['antenna_scan_mode'] == 2:
        scan_type = 'rhi'
    else:
        scan_type = 'ppi'

    # latitude
    latitude = get_metadata('latitude')
    lat = bin4_to_angle(ingest_config['latitude_radar'])
    if lat > 180.0:
        lat -= 360.0
    latitude['data'] = np.array([lat], dtype='float64')

    # longitude
    longitude = get_metadata('longitude')
    lon = bin4_to_angle(ingest_config['longitude_radar'])
    if lon > 180.0:
        lon -= 360.0
    longitude['data'] = np.array([lon], dtype='float64')

    # altitude
    altitude = get_metadata('altitude')
    alt = sigmetfile.product_hdr['product_end']['ground_height']
    altitude['data'] = np.array([alt], dtype='float64')

    # sweep_number
    sweep_number = get_metadata('sweep_number')
    sweep_number['data'] = np.arange(nsweeps, dtype='int32')

    # sweep_mode
    sweep_mode = get_metadata('sweep_mode')
    if scan_type == 'ppi':
        sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])
    else:
        sweep_mode['data'] = np.array(nsweeps * ['rhi'])

    # fixed_angle
    fixed_angle = get_metadata('fixed_angle')
    fa = [d['fixed_angle'] for d in
          sigmetfile.ingest_data_headers[first_data_type]]
    fixed_angle['data'] = bin2_to_angle(np.array(fa, dtype='float32'))

    # azimuth
    az0 = sigmet_metadata[first_data_type]['azimuth_0']
    az1 = sigmet_metadata[first_data_type]['azimuth_1']
    el0 = sigmet_metadata[first_data_type]['elevation_0']
    el1 = sigmet_metadata[first_data_type]['elevation_1']

    azimuth = get_metadata('azimuth')
    az = (az0 + az1) / 2.
    az[np.where(np.abs(az0 - az1) > 180.0)] += 180.
    az[az > 360.0] -= 360.
    azimuth['data'] = az.astype('float32')

    # elevation
    elevation = get_metadata('elevation')
    elevation['data'] = ((el0 + el1) / 2.).astype('float32')

    # instrument_parameters
    prt = get_metadata('prt')
    prt_mode = get_metadata('prt_mode')
    nyquist_velocity = get_metadata('nyquist_velocity')
    unambiguous_range = get_metadata('unambiguous_range')

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

    instrument_parameters = {'unambiguous_range': unambiguous_range,
                             'prt_mode': prt_mode, 'prt': prt,
                             'nyquist_velocity': nyquist_velocity}

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)


def _time_order_data_and_metadata(data, metadata):
    """ Put Sigmet data and metadata in time increasing order. """

    # Sigmet data is stored by sweep in azimuth increasing order,
    # to place in time increasing order we must roll each sweep so
    # the earliest collected ray is first.
    # This assuming all the fields have the same timing and the rays
    # were collected in sequentially, which appears to be true.

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


def ymds_time_to_datetime(ymds):
    """ Return a datetime object from a Sigmet ymds_time dictionary. """
    dt = datetime.datetime(ymds['year'], ymds['month'], ymds['day'])
    delta = datetime.timedelta(seconds=ymds['seconds'])
    return dt + delta


# This dictionary maps sigmet data types -> radar field names
# Users can pass their own version of this dictionary to the read_sigmet
# function as the field_names parameter.
SIGMET_TO_STANDARD = {
    #'DBT': 'DBT',                       # (1) Total Power
    #'DBZ': 'DBZ',                       # (2) Reflectivity
    #'VEL': 'VEL',                       # (3) Velocity
    #'WIDTH': 'WIDTH',                   # (4) Width
    #'ZDR': 'ZDR',                       # (5) Differential reflectivity
    #'DBZC': 'DBZC',                     # (7) Corrected reflectivity
    #'DBT2': 'DBT2',                     # (8) Total Power
    'DBZ2': 'reflectivity_horizontal',  # (9) Reflectivity
    'VEL2': 'mean_doppler_velocity',    # (10) Velocity
    #'WIDTH2': 'WIDTH2',                 # (11) Width
    'ZDR2': 'diff_reflectivity',        # (12) Differential reflectivity
    #'RAINRATE2': 'RAINRATE2',           # (13) Rainfall rate
    #'KDP': 'KDP',                       # (14) KDP (differential phase)
    'KDP2': 'diff_phase',               # (15) KDP (differential phase)
    #'PHIDP': 'PHIDP',                   # (16) PhiDP (differential phase)
    #'VELC': 'VELC',                     # (17) Corrected velocity
    #'SQI': 'SQI',                       # (18) SQI
    #'RHOHV': 'RHOHV',                   # (19) RhoHV
    'RHOHV2': 'copol_coeff',            # (20) RhoHV
    'DBZC2': 'reflectivity_horizontal_filtered',    # (21) Corrected Reflec.
    #'VELC2': 'VELC2',                   # (21) Corrected Velocity
    'SQI2': 'norm_coherent_power',      # (23) SQI
    'PHIDP2': 'dp_phase_shift',         # (24) PhiDP (differential phase)
    #'LDRH': 'LDRH',                     # (25) LDR xmt H, rcv V
    #'LDRH2': 'LDRH2',                   # (26) LDR xmt H, rcv V
    #'LDRV': 'LDRV',                     # (27) LDR xmt V, rcv H
    #'LDRV2': 'LDRV2',                   # (28) LDR xmt V, rcv H
    #'HEIGHT': 'HEIGHT',                 # (32) Height (1/10 km)
    #'VIL2': 'VIL2',                     # (33) Linear Liquid
    #'RAW': 'RAW',                       # (34) Raw Data
    #'SHEAR': 'SHEAR',                   # (35) Wind Shear
    #'DIVERGE2': 'DIVERGE2',             # (36) Divergence
    #'FLIQUID2': 'FLIQUID2',             # (37) Floated liquid
    #'USER': 'USER',                     # (38) User type
    #'OTHER': 'OTHER',                   # (39) Unspecified
    #'DEFORM2': 'DEFORM2',               # (40) Deformation
    #'VVEL2': 'VVEL2',                   # (41) Vertical velocity
    #'HVEL2': 'HVEL2',                   # (42) Horizontal velocity
    #'HDIR2': 'HDIR2',                   # (43) Horizontal wind direction
    #'AXDIL2': 'AXDIL2',                 # (44) Axis of dilation
    #'TIME2': 'TIME2',                   # (45) Time in seconds
    #'RHOH': 'RHOH',                     # (46) Rho, xmt H, rcv V
    #'RHOH2': 'RHOH2',                   # (47) Rho, xmt H, rcv V
    #'RHOV': 'RHOV',                     # (48) Rho, xmt V, rcv H
    #'RHOV2': 'RHOV2',                   # (49) Rho, xmt V, rcv H
    #'PHIH': 'PHIH',                     # (50) Phi, xmt H, rcv V
    #'PHIH2': 'PHIH2',                   # (51) Phi, xmt H, rcv V
    #'PHIV': 'PHIV',                     # (52) Phi, xmt V, rcv H
    #'PHIV2': 'PHIV2',                   # (53) Phi, xmt V, rcv H
    #'USER2': 'USER2',                   # (54) User type
    #'HCLASS': 'HCLASS',                 # (55) Hydrometeor class
    #'HCLASS2': 'HCLASS2',               # (56) Hydrometeor class
    #'ZDRC': 'ZDRC',                     # (57) Corrected diff. refl.
    #'ZDRC2': 'ZDRC2'                    # (58) Corrected diff. refl.
}
