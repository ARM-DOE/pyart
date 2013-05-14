"""
pyart.io.radar
==============

A general central radial scanning (or dwelling) instrument class.

.. autosummary::
    :toctree: generated/

    Radar
    join_radar


"""

import copy

import numpy as np


class Radar:
    """
    A class for storing antenna coordinate radar data.

    The structure of the Radar class is based on the CF/Radial Data file
    format.  Global attributes and variables (section 4.1 and 4.3) are
    represented as a dictionary in the metadata attribute.  Other required and
    optional variables are represented as dictionaries in a attribute with the
    same name as the variable in the CF/Radial standard.  When a optional
    attribute not present the attribute has a value of None.  The data for a
    given variable is stored in the dictionary under the 'data' key.  Moment
    field data is stored as a dictionary of dictionaries in the fields
    attribute.  Sub-convention variables are stored as a dictionary of
    dictionaries under the meta_group attribute.

    Refer to the attribute section for information on the parameters.

    Attributes
    ----------
    time : dict
        Time at the center of each ray.
    range : dict
        Range to the center of each gate (bin).
    fields : dict of dicts
        Moment fields.
    metadata : dict
        Metadata describing the instrument and data.
    scan_type : str
        Type of scan, one of 'ppi', 'rhi', 'sector' or 'other'.  If the scan
        volume contains multiple sweep modes this should be 'other'.
    latitude : dict
        Latitude of the instrument.
    longitude : dict
        Longitude of the instrument.
    altitude :
        Altitude of the instrument, above sea level.
    altitude_agl : dict or None
        Altitude of the instrument above ground level.  If not provided this
        attribute is set to None, indicating this parameter not available.
    sweep_number : dict
        The number of the sweep in the volume scan, 0-based.
    sweep_mode : dict
        Sweep mode for each mode in the volume scan.
    fixed_angle : dict
        Target angle for thr sweep.  Azimuth angle in RHI modes, elevation
        angle in all other modes.
    sweep_start_ray_index : dict
        Index of the first ray in each sweep relative to the start of the
        volume, 0-based.
    sweep_end_ray_index : dict
        Index of the last ray in each sweep relative to the start of the
        volume, 0-based.
    target_scan_rate : dict or None
        Intended scan rate for each sweep.  If not provided this attribute is
        set to None, indicating this parameter is not available.
    azimuth : dict
        Azimuth of antenna, relative to true North.
    elevation : dict
        Elevation of antenna, relative to the horizontal plane.
    scan_rate : dict or None
        Actual antenna scan rate.  If not provided this attribute is set to
        None, indicating this parameter is not available.
    antenna_transition : dict or None
        Flag indicating if the antenna is in transition, 1 = tes, 0 = no.
        If not provided this attribute is set to None, indicating this
        parameter is not available.
    instruments_parameters : dict of dicts or None
        Instrument parameters, if not provided this attribute is set to None,
        indicating these parameters are not avaiable.  This dictionary also
        includes variables in the radar_parameters CF/Radial subconvention.
    radar_calibration : dict of dicts or None
        Instrument calibration parameters.  If not provided this attribute is
        set to None, indicating these parameters are not available
    ngates : int
        Number of gates (bins) in the volume.
    nrays : int
        Number of rays in the volume.
    nsweeps : int
        Number of sweep in the volume.

    """

    def __init__(self, time, _range, fields, metadata, scan_type,
                 latitude, longitude, altitude,

                 sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
                 sweep_end_ray_index,

                 azimuth, elevation,

                 altitude_agl=None,
                 target_scan_rate=None,

                 scan_rate=None, antenna_transition=None,

                 instrument_parameters=None,
                 radar_calibration=None,

                 ):

        if 'calendar' not in time:
            time['calendar'] = 'gregorian'
        self.time = time
        self.range = _range

        self.fields = fields
        self.metadata = metadata
        self.scan_type = scan_type

        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude
        self.altitude_agl = altitude_agl  # optional

        self.sweep_number = sweep_number
        self.sweep_mode = sweep_mode
        self.fixed_angle = fixed_angle
        self.sweep_start_ray_index = sweep_start_ray_index
        self.sweep_end_ray_index = sweep_end_ray_index
        self.target_scan_rate = target_scan_rate  # optional

        self.azimuth = azimuth
        self.elevation = elevation
        self.scan_rate = scan_rate  # optional
        self.antenna_transition = antenna_transition  # optional

        self.instrument_parameters = instrument_parameters  # optional
        self.radar_calibration = radar_calibration  # optional

        self.ngates = len(_range['data'])
        self.nrays = len(time['data'])
        self.nsweeps = len(sweep_number['data'])


def join_radar(radar1, radar2):

    #must have same gate spacing
    new_radar = copy.deepcopy(radar1)
    new_radar.azimuth['data'] = np.append(radar1.azimuth['data'],
                                          radar2.azimuth['data'])
    new_radar.elevation['data'] = np.append(radar1.elevation['data'],
                                            radar2.elevation['data'])

    if len(radar1.range['data']) >= len(radar2.range['data']):
        new_radar.range['data'] = radar1.range['data']
    else:
        new_radar.range['data'] = radar2.range['data']
    new_radar.time['data'] = np.append(radar1.time['data'],
                                       radar2.time['data'])

    for var in new_radar.fields.keys():
        sh1 = radar1.fields[var]['data'].shape
        sh2 = radar2.fields[var]['data'].shape
        print sh1, sh2
        new_field = np.ma.zeros([sh1[0] + sh2[0],
                                max([sh1[1], sh2[1]])]) - 9999.0
        new_field[0:sh1[0], 0:sh1[1]] = radar1.fields[var]['data']
        new_field[sh1[0]:, 0:sh2[1]] = radar2.fields[var]['data']
        new_radar.fields[var]['data'] = new_field

    # This will not work for two already moving platforms..
    # need to enhance later
    if radar1.location['latitude']['data'] != radar2.location['latitude']['data'] or radar1.location['longitude']['data'] != radar2.location['longitude']['data'] or radar1.location['altitude']['data'] != radar2.location['altitude']['data']:
        for key in radar1.location.keys():
            new_radar.location[key]['data'] = np.append(
                np.zeros(len(radar1.time['data'])) +
                radar1.location[key]['data'],
                np.zeros(len(radar2.time['data'])) +
                radar2.location[key]['data'])
    return new_radar
