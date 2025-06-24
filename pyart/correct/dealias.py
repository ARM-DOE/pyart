"""
Front end to the University of Washington 4DD code for Doppler dealiasing.

"""


def dealias_fourdd(
    radar,
    last_radar=None,
    sonde_profile=None,
    gatefilter=False,
    filt=1,
    rsl_badval=131072.0,
    keep_original=False,
    set_limits=True,
    vel_field=None,
    corr_vel_field=None,
    last_vel_field=None,
    debug=False,
    max_shear=0.05,
    sign=1,
    **kwargs
):
    """
    **Deprecated in Py-ART 2.0. Please use the region-based dealiaser.

    Dealias Doppler velocities using the 4DD algorithm.

    Dealias the Doppler velocities field using the University of Washington
    4DD algorithm utilizing information from a previous volume scan and/or
    sounding data. Either last_radar or sonde_profile must be provided.
    For best results provide both a previous volume scan and sounding data.
    Radar and last_radar must contain the same number of rays per sweep.

    Additional arguments are passed to
    :py:func:`_fourdd_interface.fourdd_dealias`.
    These can be used to fine tune the behavior of the FourDD algorithm.
    See the documentation of Other Parameters for details.  For the default
    values of these parameters see the documentation of
    :py:func:`_fourdd_interface.fourdd_dealias`.

    Parameters
    ----------
    radar : Radar
        Radar object to use for dealiasing. Must have a Nyquist defined in
        the instrument_parameters attribute and have a
        reflectivity_horizontal and mean_doppler_velocity fields.
    last_radar : Radar, optional
        The previous radar volume, which has been successfully
        dealiased. Using a previous volume as an initial condition can
        greatly improve the dealiasing, and represents the final dimension
        in the 4DD algorithm.
    sonde_profile : HorizontalWindProfile, optional
        Profile of horizontal winds from a sonding used for the initial
        condition of the dealiasing.

    Other Parameters
    ----------------
    gatefilter : GateFilter, optional.
        A GateFilter instance which specifies which gates should be
        ignored when performing velocity dealiasing. A value of None will
        create this filter from the radar moments using any additional
        arguments by passing them to :py:func:`moment_based_gate_filter`. The
        default value assumes all gates are valid.
    filt : int, optional
        Flag controlling Bergen and Albers filter, 1 = yes, 0 = no.
    rsl_badval : float, optional
        Value which represents a bad value in RSL.
    keep_original : bool, optional
        True to keep original doppler velocity values when the dealiasing
        procedure fails, otherwise these gates will be masked. NaN values
        are still masked.
    set_limits : bool, optional
        True to set valid_min and valid_max elements in the returned
        dictionary. False will not set these dictionary elements.
    vel_field : str, optional
        Field in radar to use as the Doppler velocities during dealiasing.
        None will use the default field name from the Py-ART configuration
        file.
    corr_vel_field : str, optional
        Name to use for the dealiased Doppler velocity field metadata. None
        will use the default field name from the Py-ART configuration file.
    last_vel_field : str, optional
        Name to use for the dealiased Doppler velocity field metadata in
        last_radar. None will use the corr_vel_field name.
    maxshear : float, optional
        Maximum vertical shear which will be incorporated into the created
        volume from the sounding data. Parameter not used when no
        sounding data is provided.
    sign : int, optional
        Sign convention which the radial velocities in the volume created
        from the sounding data will will. This should match the convention
        used in the radar data. A value of 1 represents when positive values
        velocities are towards the radar, -1 represents when negative
        velocities are towards the radar.
    compthresh : float, optional
        Fraction of the Nyquist velocity to use as a threshold when performing
        continuity (initial) dealiasing. Velocities differences above this
        threshold will not be marked as gate from which to begin unfolding
        during spatial dealiasing.
    compthresh2 : float, optional
        The same as compthresh but the value used during the second pass of
        dealiasing. This second pass is only performed in both a sounding
        and last volume are provided.
    thresh : float, optional
        Fraction of the Nyquist velocity to use as a threshold when performing
        spatial dealiasing. Horizontally adjacent gates with velocities above
        this threshold will count against assigning the gate in question the
        velocity value being tested.
    ckval : float, optional
        When the absolute value of the velocities are below this value they
        will not be marked as gates from which to begin unfolding during
        spatial dealiasing.
    stdthresh : float, optional
       Fraction of the Nyquist velocity to use as a standard deviation
       threshold in the window dealiasing portion of the algorithm.
    epsilon : float, optional
        Difference used when comparing a value to missing value, changing this
        from the default is not recommended.
    maxcount : int, optional
        Maximum allowed number of fold allowed when unfolding velocities.
    pass2 : int, optional
        Controls weather unfolded gates should be removed (a value of 0)
        or retained for unfolding during the second pass (a value of 1) when
        both a sounding volume and last volume are provided.
    rm : int, optional
        Determines what should be done with gates that are left unfolded
        after the first pass of dealiasing. A value of 1 will remove these
        gates, a value of 0 sets these gates to their initial velocity.  If
        both a sounding volume and last volume are provided this parameter is
        ignored.
    proximity : int, optional
        Number of gates and rays to include of either side of the current gate
        during window dealiasing. This value may be doubled in cases where
        a standard sized window does not capture a sufficient number of
        good valued gates.
    mingood : int, optional
        Number of good valued gates required within the window before the
        current gate will be unfolded.
    ba_mincount : int, optional
        Number of neighbors required during Bergen and Albers filter for
        a given gate to be included, must be between 1 and 8, 5 recommended.
    ba_edgecount : int, optional
        Same as ba_mincount but used at ray edges, must be between 1 and 5,
        3 recommended.
    debug : bool, optional
        Set True to return RSL Volume objects for debugging:
        usuccess, radialVelVolume, lastVelVolume, unfoldedVolume, sondVolume

    Returns
    -------
    vr_corr : dict
        Field dictionary containing dealiased Doppler velocities. Dealiased
        array is stored under the 'data' key.

    Notes
    -----
    Due to limitations in the C code do not call with sounding arrays over
    999 elements long.

    References
    ----------
    C. N. James and R. A Houze Jr, A Real-Time Four-Dimensional Doppler
    Dealising Scheme, Journal of Atmospheric and Oceanic Technology, 2001, 18,
    1674.

    """
    # check that FourDD is available (requires TRMM RSL)

    raise DeprecationWarning(
        "FourDD dealiasing has been deprecated in Py-ART 2.0 please use dealias_region_based."
    )
