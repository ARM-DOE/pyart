"""
Find contiguous objects in scans and despeckle away ones that are too small.

"""

import numpy as np
from scipy.ndimage import label
from scipy.signal import convolve2d

from ..config import get_fillvalue
from ..filters.gatefilter import GateFilter

BAD = get_fillvalue() # Get default fill value.
DELTA = 5.0  # deg, allowable gap between PPI edges to be considered full 360
# To do:
# Testing


def find_objects(radar, field, threshold, sweeps=None, smooth=None,
                 gatefilter=None, delta=DELTA):
    """
    Find objects (i.e., contiguous gates) in one or more sweeps that match
    thresholds. Filtering & smoothing are available prior to labeling objects.
    In addition, periodic boundaries are accounted for if they exist
    (e.g., 360-deg PPIs). Requires scipy to be installed.

    Parameters
    ----------
    radar : pyart.core.Radar object
        Radar object to query.
    field : str
        Name of field to investigate for objects.
    threshold : int or float, or 2-element tuple of ints or floats
        Threshold values above (if single value) or between (if tuple)
        for objects to be identified.

    Other Parameters
    ----------------
    sweeps : int or array of ints or None, optional
        Sweep numbers to examine. If None, all sweeps are examined.
    smooth : int or None, optional
        Number of gates included in a smoothing box filter along a ray.
        If None, no smoothing is done prior to labeling objects.
    gatefilter : None or pyart.filters.GateFilter object, optional
        Py-ART GateFilter object to apply before labeling objects.
        If None, no filtering will be performed. Note: Filtering always occurs
        before smoothing.
    delta : int or float, optional
        Size of allowable gap near PPI edges, in deg, to consider it full 360.
        If gap is small, then PPI edges will be checked for matching objects
        along the periodic boundary.

    Returns
    -------
    label_dict : dict
        Dictionary that contains all the labeled objects. If this function is
        performed on the full Radar object, then the dict is ready to be added
        as a field.

    """
    if field not in radar.fields.keys():
        raise KeyError('Failed -', field, 'field not found in Radar object.')
    sweeps = _check_sweeps(sweeps, radar)
    tlo, thi = _check_threshold(threshold)
    objcnt = 0
    label_storage = []
    for iswp in sweeps:
        data = _get_data(radar, iswp, field, tlo, thi, smooth,
                         gatefilter=gatefilter)
        az = radar.get_azimuth(iswp, copy=False)
        if _check_for_360(az, delta):
            # If 360 or close, account for the periodic boundary
            labels, nobj = _adjust_for_periodic_boundary(data)
        else:
            labels, nobj = _get_labels(data)
        labels[labels != 0] += objcnt
        objcnt += nobj
        label_storage = _append_labels(labels, label_storage)
    label_storage = np.ma.masked_where(
        label_storage == 0, label_storage)
    return _generate_dict(label_storage)


def despeckle_field(radar, field, label_dict=None, threshold=-100,
                    size=10, gatefilter=None, delta=DELTA):
    """
    Despeckle a radar volume by identifying small objects in each scan and
    masking them out. User can define which field to investigate, as well as
    various thresholds to use on that field and any objects found within.
    Requires scipy to be installed, and returns a GateFilter object.

    Parameters
    ----------
    radar : pyart.core.Radar object
        Radar object to query.
    field : str
        Name of field to investigate for speckles.

    Other Parameters
    ----------------
    label_dict : dict or None, optional
        Dictionary that is produced by find_objects.
        If None, find_objects will be called to produce it.
    threshold : int or float, or 2-element tuple of ints or floats
        Threshold values above (if single value) or between (if tuple)
        for objects to be identified. Default value assumes reflectivity.
    size : int, optional
        Number of contiguous gates in an object, below which it is a speckle.
    gatefilter : None or pyart.filters.GateFilter object
        Py-ART GateFilter object to which to add the despeckling mask. The
        GateFilter object will be permanently modified with the new filtering.
        If None, creates a new GateFilter.
    delta : int or float, optional
        Size of allowable gap near PPI edges, in deg, to consider it full 360.
        If gap is small, then PPI edges will be checked for matching objects.

    Returns
    -------
    gatefilter : pyart.filters.GateFilter object
        Py-ART GateFilter object that includes the despeckling mask

    """
    if field not in radar.fields.keys():
        raise KeyError('Failed -', field, 'field not found in Radar object.')
    if label_dict is None:
        # Label everything in the radar object's field
        label_dict = find_objects(radar, field, threshold,
                                  gatefilter=gatefilter, delta=delta)
    if gatefilter is None:
        gatefilter = GateFilter(radar)
    labels = label_dict['data']

    # Get a copy of the field in the volume
    data = 1.0 * radar.fields[field]['data']
    mask_filter = gatefilter.gate_excluded
    data = np.ma.masked_array(data, mask_filter)
    data = data.filled(fill_value=BAD)
    labf = labels.filled(fill_value=0)

    # First reduce array size to speed up processing
    cond1 = np.logical_and(data != BAD, labf > 0)
    labr = labf[cond1]
    data_r = data[cond1]

    # Now loop thru all objects in volume, mask ones that are too small
    # These are the speckles
    iterarray = np.unique(labr)
    for i, lab in enumerate(iterarray):
        cond = labr == lab
        if np.size(labr[cond]) < size:
            data_r[cond] = BAD
    data[cond1] = data_r
    data = np.ma.masked_where(data == BAD, data)
    gatefilter.exclude_gates(data.mask)
    return gatefilter


def _adjust_for_periodic_boundary(data):
    """
    Identify all the contiguous objects in a sweep, accounting for the
    periodic boundary in a 360-deg PPI. Contiguous means corners or sides
    of gates touch. The algorithm appends the sweep to itself, then looks
    for contiguous objects near the original PPI edges and relabels them.
    Then, the extra sweep is discarded before returning all the labels.

    Parameters
    ----------
    data : 2D array of ints
        Sweep that will be checked for objects. Sweep has already been
        converted to binary 0s/1s based on user-supplied thresholds.

    Returns
    -------
    labels : 2D array of ints
        Numeric object labels, corrected for the periodic boundary.
        Zero values mean no object at that location.
    nobj : int
        Number of distinct objects identified in sweep.

    """
    data = np.append(data, data, axis=0)
    labels, nobj = _get_labels(data)
    i1 = 0
    # i2 = int(np.shape(labels)[0] / 2)
    i2 = labels.shape[0] // 2
    old_labs = np.unique(labels[i2][labels[i2] > 0])
    for i, lab in enumerate(old_labs):
        cond = labels == lab
        indices = np.where(labels[i2] == lab)
        new_lab = np.unique(labels[i1][indices[0]])[0]
        labels[labels == lab] = new_lab
    labels = labels[0:i2]
    nobj = len(np.unique(labels)) - 1  # labels == 0 doesnt count
    return labels, nobj


def _append_labels(labels, label_storage):
    """
    Appends consecutive sweeps of labels, creating a multi-sweep 2D array.
    Typically called iteratively.

    Parameters
    ----------
    labels : 2D array of ints
        Sweep containing object labels.
    label_storage : Empty list or 2D array of ints
        Array to append new sweep of labels to.

    Returns
    -------
    label_storage : 2D array of ints
        Updated array of object labels.

    """
    if len(label_storage) == 0:
        label_storage.append(labels)
        label_storage = np.array(label_storage[0])
    else:
        label_storage = np.append(label_storage, labels, axis=0)
    return label_storage


def _check_for_360(az, delta):
    """
    Check if an array of azimuths indicates the sweep is a full 360 PPI.
    This should also spot RHIs (effectively, a narrow azimuth sector sweep).

    Parameters
    ----------
    az : array of int or float
        Azimuths in the sweep
    delta : int or float
        Size of allowable gap near PPI edges, in deg, to consider it full 360.

    Returns
    -------
    Flag : bool
        True - Sweep is a 360 PPI.

        False - Sweep is not a 360 PPI.

    """
    # Check for small gap in azimuths
    if np.abs(az[0]-az[-1]) < delta or np.abs(az[0]-az[-1]) > 360 - delta:
        # Confirm small gap and not just narrow sector
        if np.max(az) - np.min(az) > 360 - delta:
            # Confirm not narrow sector near true north
            if True not in (np.sin(np.deg2rad(az)) <
                            np.sin(np.deg2rad(360-delta))) or \
                    True not in (np.sin(np.deg2rad(az)) >
                                 np.sin(np.deg2rad(delta))):
                return False
            else:
                return True
        else:
            return False
    else:
        return False


def _check_sweeps(sweeps, radar):
    """
    Parse the sweeps keyword and convert it to a list of ints.
    The output will be iterated over.

    Parameters
    ----------
    sweeps : int or list of ints or None
        Sweep numbers to put into an iterable list. If None, all sweeps in the
        radar object will be examined.
    radar : pyart.core.Radar object
        Radar object to query.

    Returns
    -------
    sweeps : list of ints
        Sweep numbers as an iterable list.

    """
    if sweeps is None:
        sweeps = np.arange(len(radar.sweep_number['data']))
    else:
        if hasattr(sweeps, '__len__'):
            sweeps = np.asarray(sweeps)
        else:
            sweeps = np.asarray([sweeps])
    return sweeps


def _check_threshold(threshold):
    """
    Parse the threshold keyword and return the lower and upper boundaries for
    the object search.

    Parameters
    ----------
    threshold : int or float, or 2-element tuple of ints or floats
        Threshold values above (if single value) or between (if tuple)
        for objects to be identified.

    Returns
    -------
    tlo : int or float
        Lower bound for the threshold. Values below this will not be included
        in the hunt for objects.
    thi : int or float or None
        Upper bound for the threshold. Values above this will not be included
        in the hunt for objects. None means no upper bound.

    """
    if not hasattr(threshold, '__len__'):
        threshold = np.asarray([threshold])
    if len(threshold) == 2:
        tlo = threshold[0]
        thi = threshold[1]
    elif len(threshold) > 2 or np.ndim(threshold) > 1:
        raise IndexError('Fix threshold argument! Must be single scalar ' +
                         'or 2-element tuple')
    else:
        tlo = threshold[0]
        thi = None
    return tlo, thi


def _generate_dict(label_storage):
    """
    Build the dictionary that includes all the object label information.
    If the entire Radar object was searched, the dictionary is ready to
    be added as a new field.

    Parameters
    ----------
    label_storage : 2D array of ints
        Object labels as a 2D array.

    Returns
    -------
    label_dict : dict
        Dictionary containing object labels and associated metadata.

    """
    label_dict = {}
    label_dict['data'] = label_storage
    label_dict['units'] = 'None'
    label_dict['long_name'] = 'Objects in Scan'
    label_dict['standard_name'] = 'objects_in_scan'
    label_dict['coordinates'] = 'elevation azimuth range'
    label_dict['valid_max'] = np.max(label_storage)
    label_dict['valid_min'] = 1
    return label_dict


def _get_data(radar, iswp, field, tlo, thi, window, gatefilter=None):
    """
    Get data for a field from a given sweep in a Radar object.
    Data are smoothed if desired, then converted to binary 0s/1s based
    on whether valid values are present.

    Parameters
    ----------
    radar : pyart.core.Radar object
        Radar object to query.
    iswp : int
        Sweep number to query.
    field : str
        Name of field to investigate for speckles.
    tlo : int or float
        Lower bound for the threshold. Values below this will not be included
        in the hunt for objects.
    thi : int or float or None
        Upper bound for the threshold. Values above this will not be included
        in the hunt for objects. None means no upper bound.
    window : int or None
        Number of gates included in a smoothing box filter along a ray.
        If None, no smoothing is done.

    Other Parameters
    ----------------
    gatefilter : None or pyart.filters.GateFilter object, optional
        Py-ART GateFilter object to apply before labeling objects.
        If None, no filtering will be performed.

    Returns
    -------
    data : 2D array of ints
        Sweep as array of binary 0s/1s based on whether valid values exist.

    """
    data = radar.get_field(iswp, field, copy=True)
    if gatefilter is not None:
        start, end = radar.get_start_end(iswp)
        mask_filter = gatefilter.gate_excluded[start:end+1]
        data = np.ma.masked_array(data, mask_filter)
    else:
        data = np.ma.masked_array(data)
    data = _smooth_data(data, window)
    data = data.filled(fill_value=BAD)
    if thi is None:
        cond = np.logical_or(data < tlo, data == BAD)
    else:
        cond = np.logical_or(
            data == BAD, np.logical_or(data < tlo, data > thi))
    data[cond] = 0
    data[~cond] = 1
    return data


def _get_labels(data):
    """
    Identify all the contiguous objects in a sweep. Contiguous means corners
    or sides of gates touch. Uses scipy.ndimage.label.

    Parameters
    ----------
    data : 2D array of ints
        Sweep that will be checked for objects. Sweep has already been
        converted to binary 0s/1s based on user-supplied thresholds.

    Returns
    -------
    labels : 2D array of ints
        Numeric object labels.
        Zero values mean no object at that location.
    nobj : int
        Number of distinct objects identified in sweep.

    """
    matrix = np.ones((3, 3), dtype='int16')
    labels, nobj = label(data, structure=matrix)
    return labels, nobj


def _smooth_data(data, window):
    """
    Perform box filtering along each ray of a sweep, and return the
    smoothed field. Uses scipy.signal.convolve2d which provides excellent
    performance.

    Parameters
    ----------
    data : 2D array of ints or floats
        Sweep of data for a specific field. Will be masked.
    window : int or None
        Number of gates included in a smoothing box filter along a ray.
        If None, no smoothing is done.

    Returns
    -------
    data : 2D array of ints or floats
        Smoothed sweep of data.

    """
    if window is not None:
        return np.ma.masked_array(convolve2d(
            data, np.ones((1, window))/np.float(window),
            mode='same', boundary='symm'))
    else:
        return data
