"""
pyart.correct._common_dealias
=============================

Routines used by multiple dealiasing functions.

.. autosummary::
    :toctree: generated/

    _parse_fields
    _parse_nyquist_vel
    _parse_gatefilter
    _parse_rays_wrap_around

"""

from ..config import get_field_name
from .filters import moment_based_gate_filter, GateFilter


def _parse_fields(vel_field, corr_vel_field):
    """ Parse and return the radar fields for dealiasing. """
    if vel_field is None:
        vel_field = get_field_name('velocity')
    if corr_vel_field is None:
        corr_vel_field = get_field_name('corrected_velocity')
    return vel_field, corr_vel_field


def _parse_nyquist_vel(nyquist_vel, radar, check_uniform):
    """ Parse the nyquist_vel parameter, extract from the radar if needed. """
    if nyquist_vel is None:
        nyquist_vel = [radar.get_nyquist_vel(i, check_uniform) for i in
                       range(radar.nsweeps)]
    else:   # Nyquist velocity explicitly provided
        try:
            len(nyquist_vel)
        except:   # expand single value.
            nyquist_vel = [nyquist_vel for i in range(radar.nsweeps)]
    return nyquist_vel


def _parse_gatefilter(gatefilter, radar, **kwargs):
    """ Parse the gatefilter, return a valid GateFilter object. """
    # parse the gatefilter parameter
    if gatefilter is None:  # create a moment based filter
        gatefilter = moment_based_gate_filter(radar, **kwargs)
    elif gatefilter is False:
        gatefilter = GateFilter(radar)
    else:
        gatefilter = gatefilter.copy()
    return gatefilter


def _parse_rays_wrap_around(rays_wrap_around, radar):
    """ Parse the rays_wrap_around parameter. """
    if rays_wrap_around is None:
        if radar.scan_type == 'ppi':
            rays_wrap_around = True
        else:
            rays_wrap_around = False
    return rays_wrap_around
