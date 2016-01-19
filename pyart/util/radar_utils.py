""" Functions for working radar instances. """

from __future__ import print_function

import copy

import numpy as np
from netCDF4 import num2date, date2num


def join_radar(radar1, radar2):

    # must have same gate spacing
    new_radar = copy.deepcopy(radar1)
    new_radar.azimuth['data'] = np.append(radar1.azimuth['data'],
                                          radar2.azimuth['data'])
    new_radar.elevation['data'] = np.append(radar1.elevation['data'],
                                            radar2.elevation['data'])

    if len(radar1.range['data']) >= len(radar2.range['data']):
        new_radar.range['data'] = radar1.range['data']
    else:
        new_radar.range['data'] = radar2.range['data']

    # to combine times we need to reference them to a standard
    # for this we'll use epoch time
    estring = "seconds since 1970-01-01T00:00:00Z"
    r1dt = num2date(radar1.time['data'], radar1.time['units'])
    r2dt = num2date(radar2.time['data'], radar2.time['units'])
    r1num = date2num(r1dt, estring)
    r2num = date2num(r2dt, estring)
    new_radar.time['data'] = np.append(r1num, r2num)
    new_radar.time['units'] = estring
    ### TODO Use new updated datetime_utils if accepted

    for var in new_radar.fields.keys():
        sh1 = radar1.fields[var]['data'].shape
        sh2 = radar2.fields[var]['data'].shape
#        print(sh1, sh2)
        new_field = np.ma.zeros([sh1[0] + sh2[0],
                                max([sh1[1], sh2[1]])]) - 9999.0
        new_field[0:sh1[0], 0:sh1[1]] = radar1.fields[var]['data']
        new_field[sh1[0]:, 0:sh2[1]] = radar2.fields[var]['data']
        new_radar.fields[var]['data'] = new_field

    # radar locations
    # TODO moving platforms - any more?
    if (len(radar1.latitude['data']) == 1 &
        len(radar2.latitude['data']) == 1 &
        len(radar1.longitude['data']) == 1 &
        len(radar2.longitude['data']) == 1 &
        len(radar1.altitude['data']) == 1 &
        len(radar2.altitude['data']) == 1):

        lat1 = float(radar1.latitude['data'])
        lon1 = float(radar1.longitude['data'])
        alt1 = float(radar1.altitude['data'])
        lat2 = float(radar2.latitude['data'])
        lon2 = float(radar2.longitude['data'])
        alt2 = float(radar2.altitude['data'])

        if (lat1 != lat2) or (lon1 != lon2) or (alt1 != alt2):
            ones1 = np.ones(len(radar1.time['data']), dtype='float32')
            ones2 = np.ones(len(radar2.time['data']), dtype='float32')
            new_radar.latitude['data'] = np.append(ones1 * lat1, ones2 * lat2)
            new_radar.longitude['data'] = np.append(ones1 * lon1, ones2 * lon2)
            new_radar.latitude['data'] = np.append(ones1 * alt1, ones2 * alt2)
        else:
            new_radar.latitude['data'] = radar1.latitude['data']
            new_radar.longitude['data'] = radar1.longitude['data']
            new_radar.altitude['data'] = radar1.altitude['data']

    else:
        new_radar.latitude['data'] = np.append(radar1.latitude['data'],
                                               radar2.latitude['data'])
        new_radar.longitude['data'] = np.append(radar1.longitude['data'],
                                               radar2.longitude['data'])
        new_radar.altitude['data'] = np.append(radar1.altitude['data'],
                                               radar2.altitude['data'])
    return new_radar
