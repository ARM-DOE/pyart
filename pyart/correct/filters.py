"""
pyart.correct.filters
=====================

Functions for creating filter fields (masks) which can be used it various
corrections routines in Py-ART.

.. autosummary::
    :toctree: generated/

"""

import numpy as np


class GateFilter:
    """
    A class for building a boolean arrays for filtering gates based on
    a set of condition typically based on the values in the radar fields.
    This filter can be used in various algorithms and calculations within
    Py-ART.

    See :py:func:`pyart.correct.GateFilter.exclude_below`: for method
    parameter details.

    Parameters
    ----------
    radar : Radar

    Attributes
    ----------
    gate_excluded : array, dtype=bool
        Boolean array indicating if a gate should be excluded from a
        calculation. Elements marked True indicate the corresponding gate
        should be excluded.  Those marked False should be included.

    Examples
    --------
    >>> import pyart
    >>> radar = pyart.io.read('radar_file.nc')
    >>> gate_filter = pyart.correct.GateFilter(radar)
    >>> gate_filter.exclude_below('reflectivity', 10)
    >>> gate_filter.exclude_below('normalized_coherent_power', 0.75)

    """

    def __init__(self, radar):
        """ initialize """
        self._radar = radar
        shape = (radar.nrays, radar.ngates)
        self.gate_excluded = np.zeros(shape, dtype=np.bool)

    def _get_fdata(self, field):
        """ Check that the field exists and retrieve field data. """
        self._radar.check_field_exists(field)
        return self._radar.fields[field]['data']

    def _merge(self, marked, op, exclude_masked):
        """ Merge an array of marked gates with the exclude array. """
        # exclude masked elements in marked by replacing them with the value
        # of the exclude_masked flag.  This does nothing if marked is a
        # non-masked array.
        if exclude_masked not in [True, False]:
            raise ValueError("exclude_masked must be 'True' or 'False'")
        marked = np.ma.filled(marked, exclude_masked)

        # merge array of marked gates with existing excluded gates
        # using the specified operation.
        if op == 'or':
            self.gate_excluded = np.logical_or(self.gate_excluded, marked)
        elif op == 'and':
            self.gate_excluded = np.logical_and(self.gate_excluded, marked)
        elif op == 'new':
            self.gate_excluded = marked
        else:
            raise ValueError("invalid 'op' parameter: ", op)
        return

    def exclude_below(self, field, value, exclude_masked=True, op='or'):
        """
        Exclude gates where a given field is below a given value.

        Parameters
        ----------
        field : str
            Name of field compared against the value.
        value : float
            Gates with a value below this value in the specified field will
            be maked for exclusion in the filter.
        filter_masked : bool, optional
            True to filter masked values in the specified field if the data is
            a masked array, False to include any masked values.
        op : {'and', 'or', 'new'}
            Operation to perform when merging a existing filter with the newly
            excluded gates. 'and' will perform a logical AND
            operation, 'or' a logical OR, and 'new' will replace the existing
            excluded gates with the one generated here. 'or', the default, is
            typically desired when building a list of conditions for excluding
            gates where the desired effect is to exclude gates which meet any
            of the conditions.

        """
        marked = np.ma.less(self._get_fdata(field), value)
        return self._merge(marked, op, exclude_masked)

    def exclude_above(self, field, value, exclude_masked=True, op='or'):
        """ Exclude gates where a given field is above a given value. """
        marked = self._get_fdata(field) > value
        return self._merge(marked, op, exclude_masked)

    def exclude_inside(self, field, v1, v2, exclude_masked=True, op='or'):
        """ Exclude gates where a given field is inside a given interval. """
        if v2 < v1:
            (v1, v2) = (v2, v1)
        fdata = self._get_fdata(field)
        marked = (fdata >= v1) & (fdata <= v2)
        return self._merge(marked, op, exclude_masked)

    def exclude_outside(self, field, v1, v2, exclude_masked=True, op='or'):
        """ Exclude gates where a given field is outside a given interval. """
        if v2 < v1:
            (v1, v2) = (v2, v1)
        fdata = self._get_fdata(field)
        marked = (fdata < v1) | (fdata > v2)
        return self._merge(marked, op, exclude_masked)

    def exclude_equal(self, field, value, exclude_masked=True, op='or'):
        """ Exclude gates where a given field is equal to a value. """
        marked = (self._get_fdata(field) == value)
        return self._merge(marked, op, exclude_masked)

    def exclude_not_equal(self, field, value, exclude_masked=True, op='or'):
        """ Exclude gates where a given field is not equal to a value. """
        marked = (self._get_fdata(field) != value)
        return self._merge(marked, op, exclude_masked)

    def exclude_all(self):
        """ Exclude all gates. """
        self.gate_excluded = np.ones_like(self.gate_excluded)
        return

    def exclude_none(self):
        """ Exclude no gates, include all gates. """
        self.gate_excluded = np.zeros_like(self.gate_excluded)
        return

    def exclude_masked(self, field, exclude_masked=True, op='or'):
        """ Exclude gates where a given field is masked. """
        marked = np.ma.getmaskarray(self._get_fdata(field))
        return self._merge(marked, op, exclude_masked)

    def exclude_invalid(self, field, exclude_masked=True, op='or'):
        """
        Exclude gates where an invalid value occurs in a field (NaNs or infs).
        """
        marked = ~np.isfinite(self._get_fdata(field))
        return self._merge(marked, op, exclude_masked)
