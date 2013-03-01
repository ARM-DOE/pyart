""" A general central radial scanning (or dwelling) instrument class. """

import copy

import numpy as np


class Radar:
    """
    A class for storing antenna coordinate radar data which will interact
    nicely with CF-Radial files and other pyart code

    Dictionary attributes
    ---------------------
    * azimuth
    * elevation
    * range (with additional)
    * time (with additional)
    * sweep_mode
        same as sweep_info['sweep_mode'] ...sorta
    * sweep_number
        same as sweep_info[sweep_number'] ... sorta

    Value attributes
    ----------------
    * tu - string like "seconds since ...."
        same as time['tu']
    * cal : str
        Calender (gregorian), same as time['cal']

    * nele : int
    * naz : int
    * ngates : int
    * nrays : int
    * nsweeps : int
    * scan_type : str

    Special attributes
    ------------------
    * fields
    * sweep_info : dict of dicts
    * location : dict of dicts

    #   metadata
    #   inst_params

    """

    def __init__(self, nsweeps, nrays, ngates, scan_type, naz, nele, _range,
                 azimuth, elevation, tu, cal, time, fields, sweep_info,
                 sweep_mode, sweep_number, location, inst_params, metadata):

        self.azimuth = azimuth
        self.elevation = elevation
        self.fields = fields
        self.location = location
        self.metadata = metadata
        self.naz = naz
        self.nele = nele
        self.ngates = ngates
        self.nsweeps = nsweeps
        self.range = _range
        self.scan_type = scan_type
        self.sweep_info = sweep_info
        self.sweep_mode = sweep_mode
        self.sweep_number = sweep_number
        self.time = time
        self.inst_params = inst_params

        self.cal = cal
        self.nrays = nrays
        self.tu = tu


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
