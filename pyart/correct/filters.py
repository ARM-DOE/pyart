"""
pyart.correct.filters
=====================

Functions for creating gate filters (masks) which can be used it various
corrections routines in Py-ART.

.. autosummary::
    :toctree: generated/
    :tempate: dev_template.rst

    GateFilter

"""

import numpy as np


class GateFilter(object):
    """
    A class for building a boolean arrays for filtering gates based on
    a set of condition typically based on the values in the radar fields.
    These filter can be used in various algorithms and calculations within
    Py-ART.

    See :py:func:`pyart.correct.GateFilter.exclude_below`: for method
    parameter details.

    Parameters
    ----------
    radar : Radar
        Radar object from which gate filter will be build.
    exclude_based : bool, optional
        True, the default and suggested method, will begin with all gates
        included and then use the exclude_ method to exclude gates based on
        conditions.  False will begin with all gates excluded from which
        a set of gates to include should be set using the include_ methods.

    Attributes
    ----------
    gate_excluded : array, dtype=bool
        Boolean array indicating if a gate should be excluded from a
        calculation. Elements marked True indicate the corresponding gate
        should be excluded.  Those marked False should be included.
        This is read-only attribute, any changes to the array will NOT
        be reflected in gate_included and will be lost when the attribute is
        accessed again.
    gate_included : array, dtype=bool
        Boolean array indicating if a gate should be included in a
        calculation. Elements marked True indicate the corresponding gate
        should be include.  Those marked False should be excluded.
        This is read-only attribute, any changes to the array will NOT
        be reflected in gate_excluded and will be lost when the attribute is
        accessed again.

    Examples
    --------
    >>> import pyart
    >>> radar = pyart.io.read('radar_file.nc')
    >>> gate_filter = pyart.correct.GateFilter(radar)
    >>> gate_filter.exclude_below('reflectivity', 10)
    >>> gate_filter.exclude_below('normalized_coherent_power', 0.75)

    """

    def __init__(self, radar, exclude_based=True):
        """ initialize """
        self._radar = radar
        shape = (radar.nrays, radar.ngates)
        if exclude_based:
            # start with all gates included, exclude gates based on a set
            # of rules using the exclude_ methods.
            self._gate_excluded = np.zeros(shape, dtype=np.bool)
        else:
            # start with all gates excluded, include gates based on a set
            # of rules using the include_ methods.
            self._gate_excluded = np.ones(shape, dtype=np.bool)

    # Implemetation is based on marking excluded gates stored in the private
    # _gate_excluded attribute. The gate_included attribute can be found
    # by taking the ones complement of gates_included.
    @property
    def gate_included(self):
        return ~self._gate_excluded.copy()

    @property
    def gate_excluded(self):
        return self._gate_excluded.copy()

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
            self._gate_excluded = np.logical_or(self._gate_excluded, marked)
        elif op == 'and':
            self._gate_excluded = np.logical_and(self._gate_excluded, marked)
        elif op == 'new':
            self._gate_excluded = marked
        else:
            raise ValueError("invalid 'op' parameter: ", op)
        return

    ###################
    # exclude methods #
    ###################

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
            Operation to perform when merging the existing set of excluded
            gates with the excluded gates from the current operation.
            'and' will perform a logical AND operation, 'or' a logical OR,
            and 'new' will replace the existing excluded gates with the one
            generated here. 'or', the default for exclude_ methods, is
            typically desired when building up a set of conditions for
            excluding gates where the desired effect is to exclude gates which
            meet any of the conditions. 'and', the default for include_
            methods, is typically desired when building up a set of conditions
            where the desired effect is to include gates which meet any of the
            conditions.  Note that the 'and' method MAY results in including
            gates which have peviously been excluded because they were masked
            or invalid.

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
        self._gate_excluded = np.ones_like(self._gate_excluded)
        return

    def exclude_none(self):
        """ Exclude no gates, include all gates. """
        self._gate_excluded = np.zeros_like(self._gate_excluded)
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

    ####################
    # include_ methods #
    ####################

    def include_below(self, field, value, exclude_masked=True, op='and'):
        """ Include gates where a given field is below a given value. """
        marked = self._get_fdata(field) < value
        self._merge(~marked, op, exclude_masked)

    def include_above(self, field, value, exclude_masked=True, op='and'):
        """ Include gates where a given field is above a given value. """
        marked = self._get_fdata(field) > value
        self._merge(~marked, op, exclude_masked)

    def include_inside(self, field, v1, v2, exclude_masked=True, op='and'):
        """ Include gates where a given field is inside a given interval. """
        if v2 < v1:
            (v1, v2) = (v2, v1)
        fdata = self._get_fdata(field)
        marked = (fdata >= v1) & (fdata <= v2)
        return self._merge(~marked, op, exclude_masked)

    def include_outside(self, field, v1, v2, exclude_masked=True, op='and'):
        """ Include gates where a given field is outside a given interval. """
        if v2 < v1:
            (v1, v2) = (v2, v1)
        fdata = self._get_fdata(field)
        marked = (fdata < v1) | (fdata > v2)
        return self._merge(~marked, op, exclude_masked)

    def include_equal(self, field, value, exclude_masked=True, op='and'):
        """ Include gates where a given field is equal to a value. """
        marked = (self._get_fdata(field) == value)
        return self._merge(~marked, op, exclude_masked)

    def include_not_equal(self, field, value, exclude_masked=True, op='and'):
        """ Include gates where a given field is not equal to a value. """
        marked = (self._get_fdata(field) != value)
        return self._merge(~marked, op, exclude_masked)

    def include_all(self):
        """ Include all gates. """
        self._gate_excluded = np.zeros_like(self._gate_excluded)

    def include_none(self):
        """ Include no gates, exclude all gates. """
        self._gate_excluded = np.ones_like(self._gate_excluded)

    def include_not_masked(self, field, exclude_masked=True, op='and'):
        """ Include gates where a given field in not masked. """
        marked = np.ma.getmaskarray(self._get_fdata(field))
        return self._merge(marked, op, exclude_masked)

    def include_valid(self, field, exclude_masked=True, op='and'):
        """
        Include gates where a valid value occurs in a field (not NaN or inf).
        """
        marked = np.isfinite(self._get_fdata(field))
        return self._merge(~marked, op, exclude_masked)
