"""
pyart.correct.filters
=====================

Functions for creating gate filters (masks) which can be used it various
corrections routines in Py-ART.

.. autosummary::
    :toctree: generated/

    moment_based_gate_filter
    moment_and_texture_based_gate_filter
    temp_based_gate_filter
    iso0_based_gate_filter

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    GateFilter

"""

from copy import deepcopy

import numpy as np

from ..config import get_field_name, get_metadata
from ..util import texture_along_ray


def moment_based_gate_filter(
        radar, ncp_field=None, rhv_field=None, refl_field=None,
        min_ncp=0.5, min_rhv=None, min_refl=-20., max_refl=100.0):
    """
    Create a filter which removes undesired gates based on moments.

    Creates a gate filter in which the following gates are excluded:

    * Gates where the instrument is transitioning between sweeps.
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
    gatefilter.exclude_transition()
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


def moment_and_texture_based_gate_filter(
        radar, zdr_field=None, rhv_field=None, phi_field=None, refl_field=None,
        textzdr_field=None, textrhv_field=None, textphi_field=None,
        textrefl_field=None, wind_size=7, max_textphi=20., max_textrhv=0.3,
        max_textzdr=2.85, max_textrefl=8., min_rhv=0.6):
    """
    Create a filter which removes undesired gates based on texture of moments.

    Creates a gate filter in which the following gates are excluded:

    * Gates where the instrument is transitioning between sweeps.
    * Gates where RhoHV is below min_rhv
    * Gates where the PhiDP texture is above max_textphi.
    * Gates where the RhoHV texture is above max_textrhv.
    * Gates where the ZDR texture is above max_textzdr
    * Gates where the reflectivity texture is above max_textrefl
    * If any of the thresholds is not set or the field (RhoHV, ZDR, PhiDP,
      reflectivity) do not exist in the radar the filter is not applied.

    Parameters
    ----------
    radar : Radar
        Radar object from which the gate filter will be built.
    zdr_field, rhv_field, phi_field, refl_field : str
        Names of the radar fields which contain the differential reflectivity,
        cross correlation ratio, differential phase and reflectivity from
        which the textures will be computed. A value of None for any of these
        parameters will use the default field name as defined in the Py-ART
        configuration file.
    textzdr_field, textrhv_field, textphi_field, textrefl_field : str
        Names of the radar fields given to the texture of the
        differential reflectivity, texture of the cross correlation ratio,
        texture of differential phase and texture of reflectivity. A value
        of None for any of these parameters will use the default field name
        as defined in the Py-ART configuration file
    wind_size : int
        Size of the moving window used to compute the ray texture.
    max_textphi, max_textrhv, max_textzdr, max_textrefl : float
        Maximum value for the texture of the differential phase, texture of
        RhoHV, texture of Zdr and texture of reflectivity. Gates in these
        fields above these limits as well as gates which are masked or contain
        invalid values will be excluded and not used in calculation which use
        the filter.  A value of None will disable filtering based upon the
        given field including removing masked or gates with an invalid value.
        To disable the thresholding but retain the masked and invalid filter
        set the parameter to a value above the highest value in the field.
    min_rhv : float
        Minimum value for the RhoHV. Gates below this limits as well as gates
        which are masked or contain invalid values will be excluded and not
        used in calculation which use the filter. A value of None will disable
        filtering based upon the given field including removing masked or
        gates with an invalid value. To disable the thresholding but retain
        the masked and invalid filter set the parameter to a value below the
        lowest value in the field.

    Returns
    -------
    gatefilter : :py:class:`GateFilter`
        A gate filter based upon the described criteria.  This can be
        used as a gatefilter parameter to various functions in pyart.correct.

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if zdr_field is None:
        zdr_field = get_field_name('differential_reflectivity')
    if rhv_field is None:
        rhv_field = get_field_name('cross_correlation_ratio')
    if phi_field is None:
        phi_field = get_field_name('uncorrected_differential_phase')
    if textrefl_field is None:
        textrefl_field = get_field_name('reflectivity_texture')
    if textzdr_field is None:
        textzdr_field = get_field_name('differential_reflectivity_texture')
    if textrhv_field is None:
        textrhv_field = get_field_name('cross_correlation_ratio_texture')
    if textphi_field is None:
        textphi_field = get_field_name('differential_phase_texture')

    # make deepcopy of input radar (we do not want to modify the original)
    radar_aux = deepcopy(radar)

    # compute the textures of the moments and add them into radar object
    if (max_textphi is not None) and (phi_field in radar_aux.fields):
        textphi = texture_along_ray(radar_aux, phi_field, wind_size=wind_size)
        tphi = get_metadata(textphi_field)
        tphi['data'] = textphi
        radar_aux.add_field(textphi_field, tphi)

    if (max_textrhv is not None) and (rhv_field in radar_aux.fields):
        textrho = texture_along_ray(radar_aux, rhv_field, wind_size=wind_size)
        trhv = get_metadata(textrhv_field)
        trhv['data'] = textrho
        radar_aux.add_field(textrhv_field, trhv)

    if (max_textzdr is not None) and (zdr_field in radar_aux.fields):
        textzdr = texture_along_ray(radar_aux, zdr_field, wind_size=wind_size)
        tzdr = get_metadata(textzdr_field)
        tzdr['data'] = textzdr
        radar_aux.add_field(textzdr_field, tzdr)

    if (max_textrefl is not None) and (refl_field in radar_aux.fields):
        textrefl = texture_along_ray(
            radar_aux, refl_field, wind_size=wind_size)
        trefl = get_metadata(textrefl_field)
        trefl['data'] = textrefl
        radar_aux.add_field(textrefl_field, trefl)

    # filter gates based upon field parameters
    gatefilter = GateFilter(radar_aux)
    gatefilter.exclude_transition()
    if (min_rhv is not None) and (rhv_field in radar_aux.fields):
        gatefilter.exclude_below(rhv_field, min_rhv)
        gatefilter.exclude_masked(rhv_field)
        gatefilter.exclude_invalid(rhv_field)
    if (max_textphi is not None) and (textphi_field in radar_aux.fields):
        gatefilter.exclude_above(textphi_field, max_textphi)
        gatefilter.exclude_masked(textphi_field)
        gatefilter.exclude_invalid(textphi_field)
    if (max_textrhv is not None) and (textrhv_field in radar_aux.fields):
        gatefilter.exclude_above(textrhv_field, max_textrhv)
        gatefilter.exclude_masked(textrhv_field)
        gatefilter.exclude_invalid(textrhv_field)
    if (max_textzdr is not None) and (textzdr_field in radar_aux.fields):
        gatefilter.exclude_above(textzdr_field, max_textzdr)
        gatefilter.exclude_masked(textzdr_field)
        gatefilter.exclude_invalid(textzdr_field)
    if (max_textrefl is not None) and (textrefl_field in radar_aux.fields):
        gatefilter.exclude_above(textrefl_field, max_textrefl)
        gatefilter.exclude_masked(textrefl_field)
        gatefilter.exclude_invalid(textrefl_field)
    return gatefilter


def temp_based_gate_filter(radar, temp_field=None, min_temp=0.,
                           thickness=400., beamwidth=None):
    """
    Create a filter which removes undesired gates based on temperature. Used
    primarily to filter out the melting layer and gates above it.

    Parameters
    ----------
    radar : Radar
        Radar object from which the gate filter will be built.
    temp_field : str
        Name of the radar field which contains the temperature.
        A value of None for will use the default field name as defined in
        the Py-ART configuration file.
    min_temp : float
        Minimum value for the temperature in degrees. Gates below this limits
        as well as gates which are masked or contain invalid values will be
        excluded and not used in calculation which use the filter. A value of
        None will disable filtering based upon the field including removing
        masked or gates with an invalid value. To disable the thresholding but
        retain the masked and invalid filter set the parameter to a value
        below the lowest value in the field.
    thickness : float
        The estimated thickness of the melting layer in m.
    beamwidth : float
        The radar antenna 3 dB beamwidth [deg].

    Returns
    -------
    gatefilter : :py:class:`GateFilter`
        A gate filter based upon the described criteria. This can be
        used as a gatefilter parameter to various functions in pyart.correct.

    """
    # parse the field parameters
    if temp_field is None:
        temp_field = get_field_name('temperature')

    # make deepcopy of input radar (we do not want to modify the original)
    radar_aux = deepcopy(radar)

    # filter gates based upon field parameters
    gatefilter = GateFilter(radar_aux)

    if (min_temp is not None) and (temp_field in radar_aux.fields):
        gatefilter.exclude_below(temp_field, min_temp)
        gatefilter.exclude_masked(temp_field)
        gatefilter.exclude_invalid(temp_field)

    deltar = radar.range['data'][1]-radar.range['data'][0]
    if beamwidth is not None:
        beam_rad = beamwidth*np.pi/180.
    if thickness is not None:
        temp = radar_aux.fields[temp_field]
        temp['data'] = np.ma.masked_where(
            gatefilter.gate_excluded == 1, temp['data'])
        for ray in range(radar_aux.nrays):
            gate_h_ray = radar_aux.gate_altitude['data'][ray, :]
            # index of first excluded gate
            ind_r = np.where(gatefilter.gate_excluded[ray, :] == 1)[0]
            if ind_r.size > 0:
                # some gates are excluded: find the maximum height
                ind_r = ind_r[0]
                if beamwidth is None:
                    hmax = gate_h_ray[ind_r]-thickness
                else:
                    # consider also the radar volume
                    # maximum altitude at the end of the volume
                    if ind_r < radar_aux.ngates-2:
                        hmax = (
                            (gate_h_ray[ind_r] + gate_h_ray[ind_r+1])/2. -
                            thickness)
                    else:
                        hmax = gate_h_ray[ind_r]-thickness
                    beam_radius = (
                        (radar.range['data'][ind_r]+deltar/2.)*beam_rad/2.)
                    delta_h = (
                        beam_radius
                        * np.cos(radar.elevation['data'][ray]*np.pi/180.))
                    hmax -= delta_h

                ind_hmax = np.where(
                    radar_aux.gate_altitude['data'][ray, :] > hmax)[0]
                if ind_hmax.size > 0:
                    ind_hmax = ind_hmax[0]
                    temp['data'][ray, ind_hmax:] = np.ma.masked
        radar_aux.add_field(temp_field, temp, replace_existing=True)
        gatefilter = GateFilter(radar_aux)
        gatefilter.exclude_masked(temp_field)

    return gatefilter


def iso0_based_gate_filter(radar, iso0_field=None, max_h_iso0=0.,
                           thickness=400., beamwidth=None):
    """
    Create a filter which removes undesired gates based height over the iso0.
    Used primarily to filter out the melting layer and gates above it.

    Parameters
    ----------
    radar : Radar
        Radar object from which the gate filter will be built.
    iso0_field : str
        Name of the radar field which contains the height relative to the
        iso0. A value of None for will use the default field name as defined
        in the Py-ART configuration file.
    max_h_iso0 : float
        Maximum height relative to the iso0 in m. Gates below this limits
        as well as gates which are masked or contain invalid values will be
        excluded and not used in calculation which use the filter. A value of
        None will disable filtering based upon the field including removing
        masked or gates with an invalid value. To disable the thresholding but
        retain the masked and invalid filter set the parameter to a value
        below the lowest value in the field.
    thickness : float
        The estimated thickness of the melting layer in m.
    beamwidth : float
        The radar antenna 3 dB beamwidth [deg].

    Returns
    -------
    gatefilter : :py:class:`GateFilter`
        A gate filter based upon the described criteria. This can be
        used as a gatefilter parameter to various functions in pyart.correct.

    """
    # parse the field parameters
    if iso0_field is None:
        iso0_field = get_field_name('height_over_iso0')

    # make deepcopy of input radar (we do not want to modify the original)
    radar_aux = deepcopy(radar)

    # filter gates based upon field parameters
    gatefilter = GateFilter(radar_aux)

    if (max_h_iso0 is not None) and (iso0_field in radar_aux.fields):
        gatefilter.exclude_above(iso0_field, max_h_iso0)
        gatefilter.exclude_masked(iso0_field)
        gatefilter.exclude_invalid(iso0_field)

    deltar = radar.range['data'][1]-radar.range['data'][0]
    if beamwidth is not None:
        beam_rad = beamwidth*np.pi/180.
    if thickness is not None:
        iso0 = radar_aux.fields[iso0_field]
        iso0['data'] = np.ma.masked_where(
            gatefilter.gate_excluded == 1, iso0['data'])
        for ray in range(radar_aux.nrays):
            gate_h_ray = radar_aux.gate_altitude['data'][ray, :]
            # index of first excluded gate
            ind_r = np.where(gatefilter.gate_excluded[ray, :] == 1)[0]
            if ind_r.size > 0:
                # some gates are excluded: find the maximum height
                ind_r = ind_r[0]
                if beamwidth is None:
                    hmax = gate_h_ray[ind_r]-thickness
                else:
                    # consider also the radar volume
                    # maximum altitude at the end of the volume
                    if ind_r < radar_aux.ngates-2:
                        hmax = (
                            (gate_h_ray[ind_r] + gate_h_ray[ind_r+1])/2. -
                            thickness)
                    else:
                        hmax = gate_h_ray[ind_r]-thickness
                    beam_radius = (
                        (radar.range['data'][ind_r]+deltar/2.)*beam_rad/2.)
                    delta_h = (
                        beam_radius
                        * np.cos(radar.elevation['data'][ray]*np.pi/180.))
                    hmax -= delta_h

                ind_hmax = np.where(
                    radar_aux.gate_altitude['data'][ray, :] > hmax)[0]
                if ind_hmax.size > 0:
                    ind_hmax = ind_hmax[0]
                    iso0['data'][ray, ind_hmax:] = np.ma.masked
        radar_aux.add_field(iso0_field, iso0, replace_existing=True)
        gatefilter = GateFilter(radar_aux)
        gatefilter.exclude_masked(iso0_field)

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

    def exclude_transition(self, trans_value=1, exclude_masked=True, op='or'):
        """
        Exclude all gates in rays marked as in transition between sweeps.

        Exclude all gates in rays marked as "in transition" by the
        antenna_transition attribute of the radar used to construct the filter.
        If no antenna transition information is available no gates are
        excluded.

        Parameters
        ----------
        trans_value : int, optional
            Value used in the antenna transition data to indicate that the
            instrument was between sweeps (in transition) during the collection
            of a specific ray. Typically a value of 1 is used to indicate this
            transition and the default can be used in these cases.
        exclude_masked : bool, optional
            True to filter masked values in antenna_transition if the data is
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
            gates which have previously been excluded because they were masked
            or invalid.

        """
        marked = np.zeros_like(self._gate_excluded)
        if self._radar.antenna_transition is not None:
            transition_data = self._radar.antenna_transition['data']
            in_transition = transition_data == trans_value
            marked[in_transition] = True
        return self._merge(marked, op, exclude_masked)

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
            gates which have previously been excluded because they were masked
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

    def exclude_gates(self, mask, exclude_masked=True, op='or'):
        """
        Exclude gates where a given mask is equal True.

        Parameters
        ----------
        mask : numpy array
            Boolean numpy array with same shape as a field array.
        exclude_masked : bool, optional
            True to filter masked values in the specified mask if it is
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
            gates which have previously been excluded because they were masked
            or invalid.

        """
        fdata = next(iter(self._radar.fields.values()))['data']
        if mask.shape != fdata.shape:
            raise ValueError("mask array must be the same size as a field.")
        marked = np.array(mask, dtype='bool')
        return self._merge(marked, op, exclude_masked)

    ####################
    # include_ methods #
    ####################

    def include_not_transition(
            self, trans_value=0, exclude_masked=True, op='and'):
        """
        Include all gates in rays not marked as in transition between sweeps.

        Include all gates in rays not marked as "in transition" by the
        antenna_transition attribute of the radar used to construct the filter.
        If no antenna transition information is available all gates are
        included.

        Parameters
        ----------
        trans_value : int, optional
            Value used in the antenna transition data to indicate that the
            instrument is not between sweeps (in transition) during the
            collection of a specific ray. Typically a value of 0 is used to
            indicate no transition and the default can be used in these cases.
        exclude_masked : bool, optional
            True to filter masked values in antenna_transition if the data is
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
            conditions.  Note that the 'or' method MAY results in excluding
            gates which have previously been included.

        """
        if self._radar.antenna_transition is None:
            include = np.ones_like(self._gate_excluded)  # include all gates
        else:
            include = np.zeros_like(self._gate_excluded)
            transition_data = self._radar.antenna_transition['data']
            not_in_transition = (transition_data == trans_value)
            include[not_in_transition] = True
        return self._merge(~include, op, exclude_masked)

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

    def include_gates(self, mask, exclude_masked=True, op='and'):
        """
        Include gates where a given mask is equal True.

        Parameters
        ----------
        mask : numpy array
            Boolean numpy array with same shape as a field array.
        exclude_masked : bool, optional
            True to filter masked values in the specified mask if it is
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
            conditions.  Note that the 'or' method MAY results in excluding
            gates which have previously been included.

        """
        fdata = next(iter(self._radar.fields.values()))['data']
        if mask.shape != fdata.shape:
            raise ValueError("Mask array must be the same size as a field.")
        marked = ~np.array(mask, dtype='bool')
        return self._merge(marked, op, exclude_masked)
