"""
Functions for image muting radar objects before plotting
"""
import numpy as np


def image_mute_radar(
        radar, field, mute_field, mute_threshold, field_threshold=None):
    """
    This function will split a field based on thresholds from another field

    Specifically, it was designed to separate areas of reflectivity where
    the correlation coefficient is less than a certain threshold to discern
    melting precipitation

    Parameters
    ----------
    radar : Radar
        Radar instance which provides the fields for muting
    field : str
        Name of field to image mute
    mute_field : str
        Name of field to image mute by
    mute_threshold : float
        Threshold value to mute by
    field_threshold : float
        Additional threshold to mask

    Returns
    -------
    radar : Radar
        Radar object with 2 new fields from input field, one muted and one not muted
    """

    # add checks for field availability
    if field not in radar.fields.keys():
        raise KeyError('Failed - ', field, ' field to mute not found in Radar object.')

    if mute_field not in radar.fields.keys():
        raise KeyError('Failed - ', mute_field, ' field to mute by not found in Radar object.')

    # get data from fields
    data_to_mute = radar.fields[field]['data']
    data_mute_by = radar.fields[mute_field]['data']

    # create filters
    # field_filter is used if user wants to use additional criteria in the original field
    if field_threshold is not None:
        field_filter = data_to_mute >= field_threshold
    else:
        field_filter = None

    # mute_filter will be the primary filter for determining muted regions
    mute_filter = data_mute_by <= mute_threshold

    # mute_mask is the combined filter
    if field_filter is None:
        mute_mask = mute_filter
    else:
        mute_mask = mute_filter & field_filter

    # break up the field into muted regions and non muted regions
    non_muted_field = np.ma.masked_where(mute_mask, data_to_mute)
    non_muted_field = np.ma.masked_invalid(non_muted_field)

    muted_field = np.ma.masked_where(~mute_mask, data_to_mute)
    muted_field = np.ma.masked_invalid(muted_field)

    # add fields to a dictionary and save to radar object
    non_muted_dict = radar.fields[field].copy()
    non_muted_dict['data'] = non_muted_field
    non_muted_dict['long_name'] = 'Non-muted ' + field
    non_muted_dict['standard_name'] = 'Non-muted ' + field
    radar.add_field('nonmuted_'+field, non_muted_dict)

    muted_dict = radar.fields[field].copy()
    muted_dict['data'] = muted_field
    muted_dict['long_name'] = 'Muted ' + field
    muted_dict['standard_name'] = 'Muted ' + field
    radar.add_field('muted_'+field, muted_dict)

    return radar

