"""
Retrieval of QVPs from a radar object

"""

import numpy as np

from ..core.transforms import antenna_to_cartesian


def quasi_vertical_profile(radar, desired_angle=None, fields=None, gatefilter=None):

    """
    Quasi Vertical Profile.

    Creates a QVP object containing fields from a radar object that can
    be used to plot and produce the quasi vertical profile

    Parameters
    ----------
    radar : Radar
        Radar object used.
    field : string
        Radar field to use for QVP calculation.
    desired_angle : float
        Radar tilt angle to use for indexing radar field data.
        None will result in wanted_angle = 20.0

    Other Parameters
    ----------------
    gatefilter : GateFilter
        A GateFilter indicating radar gates that should be excluded
        from the import qvp calculation

    Returns
    -------
    qvp : Dictonary
        A quasai vertical profile object containing fields
        from a radar object

    References
    ----------
    Troemel, S., M. Kumjian, A. Ryzhkov, and C. Simmer, 2013: Backscatter
    differential phase - estimation and variability. J Appl. Meteor. Clim..
    52, 2529 - 2548.

    Troemel, S., A. Ryzhkov, P. Zhang, and C. Simmer, 2014: Investigations
    of backscatter differential phase in the melting layer. J. Appl. Meteorol.
    Clim. 54, 2344 - 2359.

    Ryzhkov, A., P. Zhang, H. Reeves, M. Kumjian, T. Tschallener, S. Tromel,
    C. Simmer, 2015: Quasi-vertical profiles - a new way to look at polarimetric
    radar data. Submitted to J. Atmos. Oceanic Technol.

    """

    # Creating an empty dictonary
    qvp = {}

    # Setting the desired radar angle and getting index value for desired radar angle
    if desired_angle is None:
        desired_angle = 20.0
    index = abs(radar.fixed_angle['data'] - desired_angle).argmin()
    radar_slice = radar.get_slice(index)

    # Printing radar tilt angles and radar elevation
    print(radar.fixed_angle['data'])
    print(radar.elevation['data'][-1])

    # Setting field parameters
    # If fields is None then all radar fields pulled else defined field is used
    if fields is None:
        fields = radar.fields

        for field in fields:

            # Filtering data based on defined gatefilter
            # If none is defined goes to else statement
            if gatefilter is not None:
                get_fields = radar.get_field(index, field)
                mask_fields = np.ma.masked_where(gatefilter.gate_excluded[radar_slice],
                                                 get_fields)
                radar_fields = np.ma.mean(mask_fields, axis=0)
            else:
                radar_fields = radar.get_field(index, field).mean(axis=0)

            qvp.update({field:radar_fields})

    else:
        # Filtereing data based on defined gatefilter
        # If none is defined goes to else statement
        if gatefilter is not None:
            get_field = radar.get_field(index, fields)
            mask_field = np.ma.masked_where(gatefilter.gate_excluded[radar_slice],
                                            get_field)
            radar_field = np.ma.mean(mask_field, axis=0)
        else:
            radar_field = radar.get_field(index, fields).mean(axis=0)

        qvp.update({fields:radar_field})

    # Adding range, time, and height fields
    qvp.update({'range': radar.range['data'], 'time': radar.time})
    _, _, z = antenna_to_cartesian(qvp['range']/1000.0, 0.0,
                                   radar.fixed_angle['data'][index])
    qvp.update({'height': z})
    return qvp
