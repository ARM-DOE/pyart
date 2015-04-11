"""
pyart.correct.filters
=====================

Functions for creating gate filters (masks) which can be used it various
corrections routines in Py-ART.

.. autosummary::
    :toctree: generated/

    moment_based_gate_filter

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    GateFilter

"""

import numpy as np

from ..config import get_field_name


def moment_based_gate_filter(
        radar, ncp_field=None, rhv_field=None, refl_field=None,
        min_ncp=0.5, min_rhv=None, min_refl=-20., max_refl=100.0):
    """
    Create a filter which removes undesired gates based on moments.

    Creates a gate filter in which the following gates are excluded:

    * Gates where the reflectivity is outside the interval min_refl, max_refl.
    * Gates where the normalized coherent power is below min_ncp.
    * Gates where the cross correlation ratio is below min_rhi.  Using the
      default parameter this filtering is disabled.
    * Gates where any of the above three fields are masked or contain
      invalid values (NaNs or infs).
    * If any of these three fields do not exist in the radar that fields filter
      criteria is not applied.

    Parameters
    ----------
    radar : Radar
        Radar object from which the gate filter will be built.
    refl_field, ncp_field, rhv_field : str
        Names of the radar fields which contain the reflectivity, normalized
        coherent power (signal quality index) and cross correlation ratio
        (RhoHV) from which the gate filter will be created using the above
        criteria.  A value of None for any of these parameters will use the
        default field name as defined in the Py-ART configuration file.
    min_ncp, min_rhv : float
        Minimum values for the normalized coherence power and cross
        correlation ratio.  Gates in these fields below these limits as well as
        gates which are masked or contain invalid values will be excluded and
        not used in calculation which use the filter.  A value of None will
        disable filtering based upon the given field including removing
        masked or gates with an invalid value.  To disable the thresholding
        but retain the masked and invalid filter set the parameter to a value
        below the lowest value in the field.
    min_refl, max_refl : float
        Minimum and maximum values for the reflectivity.  Gates outside
        of this interval as well as gates which are masked or contain invalid
        values will be excluded and not used in calculation which use this
        filter. A value or None for one of these parameters will disable the
        minimum or maximum filtering but retain the other.  A value of None
        for both of these values will disable all filtering based upon the
        reflectivity including removing masked or gates with an invalid value.
        To disable the interval filtering but retain the masked and invalid
        filter set the parameters to values above and below the lowest and
        greatest values in the reflectivity field.

    Returns
    -------
    gatefilter : :py:class:`GateFilter`
        A gate filter based upon the described criteria.  This can be
        used as a gatefilter parameter to various functions in pyart.correct.

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')
    if rhv_field is None:
        rhv_field = get_field_name('cross_correlation_ratio')

    # filter gates based upon field parameters
    gatefilter = GateFilter(radar)
    if (min_ncp is not None) and (ncp_field in radar.fields):
        gatefilter.exclude_below(ncp_field, min_ncp)
        gatefilter.exclude_masked(ncp_field)
        gatefilter.exclude_invalid(ncp_field)
    if (min_rhv is not None) and (rhv_field in radar.fields):
        gatefilter.exclude_below(rhv_field, min_rhv)
        gatefilter.exclude_masked(rhv_field)
        gatefilter.exclude_invalid(rhv_field)
    if refl_field in radar.fields:
        if min_refl is not None:
            gatefilter.exclude_below(refl_field, min_refl)
            gatefilter.exclude_masked(refl_field)
            gatefilter.exclude_invalid(refl_field)
        if max_refl is not None:
            gatefilter.exclude_above(refl_field, max_refl)
            gatefilter.exclude_masked(refl_field)
            gatefilter.exclude_invalid(refl_field)
    return gatefilter


class GateFilter(object):
    """
    A class for building a boolean arrays for filtering gates based on
    a set of condition typically based on the values in the radar fields.
    These filter can be used in various algorithms and calculations within
    Py-ART.

    See :py:func:`pyart.correct.GateFilter.exclude_below` for method
    parameter details.

    Parameters
    ----------
    radar : Radar
        Radar object from which gate filter will be build.
    exclude_based : bool, optional
        True, the default and suggested method, will begin with all gates
        included and then use the exclude methods to exclude gates based on
        conditions.  False will begin with all gates excluded from which
        a set of gates to include should be set using the include methods.

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
    >>> gatefilter = pyart.correct.GateFilter(radar)
    >>> gatefilter.exclude_below('reflectivity', 10)
    >>> gatefilter.exclude_below('normalized_coherent_power', 0.75)

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

    def copy(self):
        """ Return a copy of the gatefilter. """
        a = GateFilter(self._radar)
        a._gate_excluded = self._gate_excluded.copy()
        return a

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

    def exclude_below(self, field, value, exclude_masked=True, op='or',
                      inclusive=False):
        """
        Exclude gates where a given field is below a given value.

        Parameters
        ----------
        field : str
            Name of field compared against the value.
        value : float
            Gates with a value below this value in the specified field will
            be marked for exclusion in the filter.
        exclude_masked : bool, optional
            True to filter masked values in the specified field if the data is
            a masked array, False to include any masked values.
        op : {'and', 'or', 'new'}
            Operation to perform when merging the existing set of excluded
            gates with the excluded gates from the current operation.
            'and' will perform a logical AND operation, 'or' a logical OR,
            and 'new' will replace the existing excluded gates with the one
            generated here. 'or', the default for exclude methods, is
            typically desired when building up a set of conditions for
            excluding gates where the desired effect is to exclude gates which
            meet any of the conditions. 'and', the default for include
            methods, is typically desired when building up a set of conditions
            where the desired effect is to include gates which meet any of the
            conditions.  Note that the 'and' method MAY results in including
            gates which have peviously been excluded because they were masked
            or invalid.
        inclusive : bool
            Indicates whether the specified value should also be excluded.

        """
        if inclusive:
            marked = self._get_fdata(field) <= value
        else:
            marked = self._get_fdata(field) < value
        return self._merge(marked, op, exclude_masked)

    def exclude_above(self, field, value, exclude_masked=True, op='or',
                      inclusive=False):
        """ Exclude gates where a given field is above a given value. """
        if inclusive:
            marked = self._get_fdata(field) >= value
        else:
            marked = self._get_fdata(field) > value
        return self._merge(marked, op, exclude_masked)

    def exclude_inside(self, field, v1, v2, exclude_masked=True, op='or',
                       inclusive=True):
        """ Exclude gates where a given field is inside a given interval. """
        if v2 < v1:
            (v1, v2) = (v2, v1)
        fdata = self._get_fdata(field)
        if inclusive:
            marked = (fdata >= v1) & (fdata <= v2)
        else:
            marked = (fdata > v1) & (fdata < v2)
        return self._merge(marked, op, exclude_masked)

    def exclude_outside(self, field, v1, v2, exclude_masked=True, op='or',
                        inclusive=False):
        """ Exclude gates where a given field is outside a given interval. """
        if v2 < v1:
            (v1, v2) = (v2, v1)
        fdata = self._get_fdata(field)
        if inclusive:
            marked = (fdata <= v1) | (fdata >= v2)
        else:
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

    def include_below(self, field, value, exclude_masked=True, op='and',
                      inclusive=False):
        """ Include gates where a given field is below a given value. """
        if inclusive:
            marked = self._get_fdata(field) <= value
        else:
            marked = self._get_fdata(field) < value
        self._merge(~marked, op, exclude_masked)

    def include_above(self, field, value, exclude_masked=True, op='and',
                      inclusive=False):
        """ Include gates where a given field is above a given value. """
        if inclusive:
            marked = self._get_fdata(field) >= value
        else:
            marked = self._get_fdata(field) > value
        self._merge(~marked, op, exclude_masked)

    def include_inside(self, field, v1, v2, exclude_masked=True, op='and',
                       inclusive=True):
        """ Include gates where a given field is inside a given interval. """
        if v2 < v1:
            (v1, v2) = (v2, v1)
        fdata = self._get_fdata(field)
        if inclusive:
            marked = (fdata >= v1) & (fdata <= v2)
        else:
            marked = (fdata > v1) & (fdata < v2)
        return self._merge(~marked, op, exclude_masked)

    def include_outside(self, field, v1, v2, exclude_masked=True, op='and',
                        inclusive=False):
        """ Include gates where a given field is outside a given interval. """
        if v2 < v1:
            (v1, v2) = (v2, v1)
        fdata = self._get_fdata(field)
        if inclusive:
            marked = (fdata <= v1) | (fdata >= v2)
        else:
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
