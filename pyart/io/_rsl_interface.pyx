"""
pyart.io._rsl_interface
=======================

Cython wrapper around the NASA TRMM RSL library.

.. autosummary::
    :toctree: generated/

    copy_volume
    create_volume
    _label_volume
    fourdd_dealias

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    RslFile
    _RslVolume
    _RslSweep
    _RslRay


"""

# TODO
# NULL checking, dealloc
# seperate _rsl_interface.pxd into rsh_h and fourdd_h
# seperate _rsl_interface.pyx into _rsl_interface and _fourdd_interface
# better FourDD dealias
# create volume with correct f/invf routines

cimport _rsl_interface
import numpy as np
cimport numpy as np
from datetime import datetime


cpdef copy_volume(_RslVolume volume):
    """
    copy_volume(volume)

    Return a copy of a _RslVolume object.

    Parameters
    ----------
    volume : _RslVolume
        _RslVolume object to create a copy of.

    Returns
    -------
    nvolume : _RslVolume
        Copy of volume.

    """
    volume_copy = _rsl_interface.RSL_copy_volume(volume._Volume)
    rslvolume = _RslVolume()
    rslvolume.load(volume_copy)
    return rslvolume


cpdef create_volume(np.ndarray[np.float32_t, ndim=3] arr):
    """
    create_volume(arr)

    Create a _RslVolume object from a 3D float32 array.

    No headers parameters except nsweeps, nrays and nbins are not set in the
    resulting _RslVolume object.

    Parameters
    ----------
    arr : array, 3D, float32
        Three dimensional float32 array.

    Returns
    -------
    volumes : _RslVolume
        _RslVolume containing array data.

    """
    cdef _rsl_interface.Volume * volume
    cdef _rsl_interface.Sweep * sweep
    cdef _rsl_interface.Ray * ray

    nsweeps = arr.shape[0]
    nrays = arr.shape[1]
    nbins = arr.shape[2]

    volume = _rsl_interface.RSL_new_volume(nsweeps)
    volume.h.nsweeps = nsweeps

    for nsweep in range(nsweeps):
        sweep = _rsl_interface.RSL_new_sweep(nrays)
        sweep.h.nrays = nrays
        volume.sweep[nsweep] = sweep

        for nray in range(nrays):
            ray = _rsl_interface.RSL_new_ray(nbins)
            sweep.ray[nray] = ray
            ray.h.nbins = nbins
            ray.h.f = _rsl_interface.VR_F   # XXX this should be selectable
            ray.h.invf = _rsl_interface.VR_INVF

            for nbin in range(nbins):
                ray.range[nbin] = ray.h.invf(arr[nsweep, nray, nbin])

    rslvolume = _RslVolume()
    rslvolume.load(volume)
    return rslvolume

cpdef _label_volume(_RslVolume volume, radar):
    """
    _label_volume(volume, radar)

    Add labels for dealiasing to a _RslVolume object from a radar object.

    This function does not set all parameter in the _RslVolume suitable for
    writing out the volume, rather it set those parameters which must be set
    prior to using :py:func:`fourdd_dealias`.

    Parameters
    ----------
    volume : _RslVolume
        Volume object to which parameters will be set as needed prior to
        dealiasing.  Object is manipulated in-place.
    radar : Radar
        Radar object from which parameters are taken.

    """

    cdef _rsl_interface.Sweep * sweep
    cdef _rsl_interface.Ray * ray

    vol = volume._Volume
    nsweeps = vol.h.nsweeps
    nrays = vol.sweep[0].h.nrays
    nbins = vol.sweep[0].ray[0].h.nbins

    gate_size = int(radar.range['meters_between_gates'])
    range_bin1 = int(radar.range['meters_to_center_of_first_gate'])
    if 'shape' in dir(radar.location['altitude']):
        alt = radar.location['altitude']['data'][0]
    else:
        alt = radar.location['altitude']['data']

    nyq_vels = radar.inst_params['nyquist_velocity']['data']
    azimuths = radar.azimuth['data']
    elevs = radar.elevation['data']

    # label the volume
    vol.h.nsweeps = nsweeps

    for nsweep in range(nsweeps):
        sweep = vol.sweep[nsweep]
        for nray in range(nrays):
            ray = sweep.ray[nray]
            ray_index = nsweep * nrays + nray
            ray.h.azimuth = azimuths[ray_index]
            ray.h.elev = elevs[ray_index]
            ray.h.nyq_vel = nyq_vels[ray_index]

            ray.h.range_bin1 = range_bin1
            ray.h.gate_size = gate_size
            ray.h.alt = alt
    return


cpdef fourdd_dealias(_RslVolume DZvolume, _RslVolume radialVelVolume,
                     np.ndarray[np.float32_t, ndim=1] hc,
                     np.ndarray[np.float32_t, ndim=1] sc,
                     np.ndarray[np.float32_t, ndim=1] dc,
                     vad_time, prep, filt):
    """
    fourdd_dealias(DZvolume, radialVelVolume, hc, sc, dc, vad_time, prep,
                   filt)

    Dealias using the FourDD algorithm.

    Parameters
    ----------
    DZvolume : _RslVolume
        Reflectivity to use when thresholding is selected.
    radialVelVolume : _RslVolume
        Radial velocities which will be dealiased.
    hc : array
        Sounding heights in meters.  Must be a contiguous one-dimensional
        float32 array.
    sc : array
        Sounding wind speed in m/s.  Must be a contiguous one-dimensional
        float32 array.
    dc : array
        Sounding wind direction in degrees.  Must be a contiguous
        one-dimensional float32 array.
    vad_time : int
        Time of sounding in YYDDDHHMM format.
    prep : int
        Flag controlling thresholding of DZvolume, 1 = yes, 0 = no.
    filt : int
        Flag controlling Bergen and Albers filter, 1 = yes, 0 = no.

    Returns
    -------
    usuccess : int
        Flag indicating if the unfolding was successful, 1 = yes, 0 = no.
    data : array
        Array of unfolded velocities.

    References
    ----------
    C. N. James and R. A Houze Jr, A Real-Time Four-Dimensional Doppler
    Dealising Scheme, Journal of Atmospheric and Oceanic Technology, 2001, 18,
    1674.

    """
    # TODO version which does not need DZvolume (prep is 0)
    # TODO version which uses
    # See FourDD.c for addition details

    cdef _RslVolume unfoldedVolume = copy_volume(radialVelVolume)
    cdef _RslVolume sondVolume = copy_volume(radialVelVolume)

    cdef float MISSINGVEL = 131072.0
    cdef unsigned short success = 0
    cdef unsigned short usuccess = 0

    # May not always be needed...
    _rsl_interface.firstGuessNoRead(
        sondVolume._Volume, MISSINGVEL, <float *> hc.data, <float *> sc.data,
        <float *> dc.data, <int> len(hc), vad_time, &success)

    if success != 1:
        raise ValueError

    # dealias
    if prep == 1:
        _rsl_interface.prepVolume(DZvolume._Volume, unfoldedVolume._Volume,
                                  MISSINGVEL)
    _rsl_interface.unfoldVolume(
        unfoldedVolume._Volume, sondVolume._Volume, NULL,
        MISSINGVEL, filt, &usuccess)
    data = unfoldedVolume.get_data()
    return usuccess, data


cdef class _RslRay:
    """
    A object for accessing RSL Ray data and header information

    This class should not be initalized from within Python.  _RslRay object are
    returned from the :py:func:`_RslSweep.get_ray` method.

    Attributes
    ----------
    month : int
        Date for this ray, month (1-12).
    day : int
        Date for this ray, day (1-31).
    year : int
        Date for this ray, year (eg. 1993).
    hour : int
        Time for this ray, hour (0-23).
    minute : int
        Time for this ray, minute (0-59).
    sec : float
        Time for this ray, second + fractor of second.
    unam_rng : float
        Unambiguous range in km.
    azimuth : float
        Mean azimuth angle in degrees for the ray, must be positive.
        0 for north, 90 for east, 270 for west.
    ray_num : int
        Ray number within a scan.
    elev : float
        Elevation angle in degrees.
    elev_num : int
        Elevation number within the volume scan.
    range_bin1 : int
        Range to first gate in meters.
    gate_size : int
        Gate size in meters.
    vel_res : float
        Doppler velocity resolution.
    sweep_rate : float
        Sweep rate, full sweeps / minute.
    prf : int
        Pulse repetition frequency in Hz.
    prf2 : int
        Second pulse repition frequenct for dual PRF data.
    azim_rate : float
        Sweep rate in degrees / second.
    fix_angle : float
        Elevation angle for the sweep in degrees.
    pitch : float
        Pitch angle.
    roll : float
        Roll angle.
    heading : float
        Heading.
    pitch_rate : float
        Pitch rate in angle / sec.
    roll_rate : float
        Roll rate in angle / sec.
    heading_rate : float
        Heading rate in angle / sec.
    lat : float
        Latitude in degrees.
    lon : float
        Longitude in degrees.
    alt : int
        Altitude in meters.
    rvs : float
        Radial velocity correction in meters / second.
    vel_east : float
        Platform velocity to the east in meters / second.  Negative values for
        velocity to the west.
    vel_north : float
        Platform velocity to the north in meters / second.  Negative values for
        velocity to the south.
    vel_up : float
        Platform velocity upward in meters / second.  Negative values for
        velocity downward.
    pulse_count : int
        Pulses used in a single dwell time.
    pulse_width : float
        Pulse width in microseconds.
    beam_width : float
        Beamwidth in degrees.
    frequency : float
        Carrier frequency in GHz.
    wavelength : float
        Wavelength in meters.
    nyq_vel : float
        Nyquist velocity in meters / second.
    nbins : int
        Number of array elements in ray data.

    """

    cdef _rsl_interface.Ray * _Ray

    cdef load(self, _rsl_interface.Ray * Ray):
        """ Load the _RslRay object, must be called after creation. """
        self._Ray = Ray

    def get_datetime(self):
        """
        get_datetime()

        Return a datetime describing the date and time of the ray.
        """
        s = self
        return datetime(s.year, s.month, s.day, s.hour, s.minute, int(s.sec))

    def get_data(self):
        """
        get_data()

        Return the one-dimensional data contained in the ray.
        """
        cdef _rsl_interface.Range raw
        cdef np.ndarray[np.float32_t, ndim = 1] data

        shape = (self._Ray.h.nbins)
        data = np.zeros(shape, dtype='float32') + 1.31072000e+05
        for nbin in range(self._Ray.h.nbins):
            raw = self._Ray.range[nbin]
            data[nbin] = self._Ray.h.f(raw)
        return data

    # header properties mapped to class attributes.
    property month:
        def __get__(self):
            return self._Ray.h.month

        def __set__(self, int month):
            self._Ray.h.month = month

    property day:
        def __get__(self):
            return self._Ray.h.day

        def __set__(self, int day):
            self._Ray.h.day = day

    property year:
        def __get__(self):
            return self._Ray.h.year

        def __set__(self, int year):
            self._Ray.h.year = year

    property hour:
        def __get__(self):
            return self._Ray.h.hour

        def __set__(self, int hour):
            self._Ray.h.hour = hour

    property minute:
        def __get__(self):
            return self._Ray.h.minute

        def __set__(self, int minute):
            self._Ray.h.minute = minute

    property sec:
        def __get__(self):
            return self._Ray.h.sec

        def __set__(self, float sec):
            self._Ray.h.sec = sec

    property unam_rng:
        def __get__(self):
            return self._Ray.h.unam_rng

        def __set__(self, float unam_rng):
            self._Ray.h.unam_rng = unam_rng

    property azimuth:
        def __get__(self):
            return self._Ray.h.azimuth

        def __set__(self, float azimuth):
            self._Ray.h.azimuth = azimuth

    property ray_num:
        def __get__(self):
            return self._Ray.h.ray_num

        def __set__(self, int ray_num):
            self._Ray.h.ray_num = ray_num

    property elev:
        def __get__(self):
            return self._Ray.h.elev

        def __set__(self, float elev):
            self._Ray.h.elev = elev

    property elev_num:
        def __get__(self):
            return self._Ray.h.elev_num

        def __set__(self, int elev_num):
            self._Ray.h.elev_num = elev_num

    property range_bin1:
        def __get__(self):
            return self._Ray.h.range_bin1

        def __set__(self, int range_bin1):
            self._Ray.h.range_bin1 = range_bin1

    property gate_size:
        def __get__(self):
            return self._Ray.h.gate_size

        def __set__(self, int gate_size):
            self._Ray.h.gate_size = gate_size

    property vel_res:
        def __get__(self):
            return self._Ray.h.vel_res

        def __set__(self, float vel_res):
            self._Ray.h.vel_res = vel_res

    property sweep_rate:
        def __get__(self):
            return self._Ray.h.sweep_rate

        def __set__(self, float sweep_rate):
            self._Ray.h.sweep_rate = sweep_rate

    property prf:
        def __get__(self):
            return self._Ray.h.prf

        def __set__(self, int prf):
            self._Ray.h.prf = prf

    property prf2:
        def __get__(self):
            return self._Ray.h.prf2

        def __set__(self, int prf2):
            self._Ray.h.prf2 = prf2

    property azim_rate:
        def __get__(self):
            return self._Ray.h.azim_rate

        def __set__(self, float azim_rate):
            self._Ray.h.azim_rate = azim_rate

    property fix_angle:
        def __get__(self):
            return self._Ray.h.fix_angle

        def __set__(self, float fix_angle):
            self._Ray.h.fix_angle = fix_angle

    property pitch:
        def __get__(self):
            return self._Ray.h.pitch

        def __set__(self, float pitch):
            self._Ray.h.pitch = pitch

    property roll:
        def __get__(self):
            return self._Ray.h.roll

        def __set__(self, float roll):
            self._Ray.h.roll = roll

    property heading:
        def __get__(self):
            return self._Ray.h.heading

        def __set__(self, float heading):
            self._Ray.h.heading = heading

    property pitch_rate:
        def __get__(self):
            return self._Ray.h.pitch_rate

        def __set__(self, float pitch_rate):
            self._Ray.h.pitch_rate = pitch_rate

    property roll_rate:
        def __get__(self):
            return self._Ray.h.roll_rate

        def __set__(self, float roll_rate):
            self._Ray.h.roll_rate = roll_rate

    property heading_rate:
        def __get__(self):
            return self._Ray.h.heading_rate

        def __set__(self, float heading_rate):
            self._Ray.h.heading_rate = heading_rate

    property lat:
        def __get__(self):
            return self._Ray.h.lat

        def __set__(self, float lat):
            self._Ray.h.lat = lat

    property lon:
        def __get__(self):
            return self._Ray.h.lon

        def __set__(self, float lon):
            self._Ray.h.lon = lon

    property alt:
        def __get__(self):
            return self._Ray.h.alt

        def __set__(self, int alt):
            self._Ray.h.alt = alt

    property rvc:
        def __get__(self):
            return self._Ray.h.rvc

        def __set__(self, float rvc):
            self._Ray.h.rvc = rvc

    property vel_east:
        def __get__(self):
            return self._Ray.h.vel_east

        def __set__(self, float vel_east):
            self._Ray.h.vel_east = vel_east

    property vel_north:
        def __get__(self):
            return self._Ray.h.vel_north

        def __set__(self, float vel_north):
            self._Ray.h.vel_north = vel_north

    property vel_up:
        def __get__(self):
            return self._Ray.h.vel_up

        def __set__(self, float vel_up):
            self._Ray.h.vel_up = vel_up

    property pulse_count:
        def __get__(self):
            return self._Ray.h.pulse_count

        def __set__(self, int pulse_count):
            self._Ray.h.pulse_count = pulse_count

    property pulse_width:
        def __get__(self):
            return self._Ray.h.pulse_width

        def __set__(self, float pulse_width):
            self._Ray.h.pulse_width = pulse_width

    property beam_width:
        def __get__(self):
            return self._Ray.h.beam_width

        def __set__(self, float beam_width):
            self._Ray.h.beam_width = beam_width

    property frequency:
        def __get__(self):
            return self._Ray.h.frequency

        def __set__(self, float frequency):
            self._Ray.h.frequency = frequency

    property wavelength:
        def __get__(self):
            return self._Ray.h.wavelength

        def __set__(self, float wavelength):
            self._Ray.h.wavelength = wavelength

    property nyq_vel:
        def __get__(self):
            return self._Ray.h.nyq_vel

        def __set__(self, float nyq_vel):
            self._Ray.h.nyq_vel = nyq_vel

    property nbins:
        def __get__(self):
            return self._Ray.h.nbins

        def __set__(self, int nbins):
            self._Ray.h.nbins = nbins


cdef class _RslSweep:
    """
    A object for accessing RSL Sweep data and header information.

    This class should not be initalized from within Python.  _RslSweep objects
    are returned from the :py:func:`_RslVolume.get_sweep` method.

    Attributes
    ----------
    sweep_num : int
        Interger sweep number.
    elev : float
        Mean elevation angle for thr sweep.  -999.0 for RHI sweeps.
    azimuth : float
        Azumuth for the sweep.  -999.0 for PPI scans.
    beam_width : float
        Beam width in degrees.  Can also be found in _RslRay objects.
    vert_half_bw : float
        Vertical beam width divided by 2.
    horz_half_bw : float
        Horizontal beam width divided by 2.
    nrays : int
        Number of rays in the sweep.

    """

    cdef _rsl_interface.Sweep * _Sweep

    cdef load(self, _rsl_interface.Sweep * Sweep):
        """ Load the _RslSweep object, must be called after creation. """
        self._Sweep = Sweep

    def get_ray(self, int ray_number):
        """
        get_ray(ray_number)

        Return a _RslRay for a given ray.

        Parameters
        ----------
        ray_number : int
            Ray number to retrieve

        Returns
        -------
        ray : _RslRay
            _RslRay object containing the requested ray.

        """
        rslray = _RslRay()
        rslray.load(self._Sweep.ray[ray_number])
        return rslray

    def get_data(self):
        """
        get_data()

        Return the two-dimensional data contained in the sweep.
        """
        cdef _rsl_interface.Range raw
        cdef np.ndarray[np.float32_t, ndim = 2] data

        sweep = self._Sweep
        nrays = sweep.h.nrays
        nbins = sweep.ray[0].h.nbins

        shape = (nrays, nbins)
        data = np.zeros(shape, dtype='float32') + 1.31072000e+05
        for nray in range(nrays):
            ray = sweep.ray[nray]
            nbins = ray.h.nbins
            for nbin in range(nbins):
                raw = ray.range[nbin]
                data[nray, nbin] = ray.h.f(raw)
        return data

    # header properties mapped to class attributes.
    property sweep_num:
        def __get__(self):
            return self._Sweep.h.sweep_num

        def __set__(self, int sweep_num):
            self._Sweep.h.sweep_num = sweep_num

    property elev:
        def __get__(self):
            return self._Sweep.h.elev

        def __set__(self, float elev):
            self._Sweep.h.elev = elev

    property azimuth:
        def __get__(self):
            return self._Sweep.h.azimuth

        def __set__(self, float azimuth):
            self._Sweep.h.azimuth = azimuth

    property beam_width:
        def __get__(self):
            return self._Sweep.h.beam_width

        def __set__(self, float beam_width):
            self._Sweep.h.beam_width = beam_width

    property vert_half_bw:
        def __get__(self):
            return self._Sweep.h.vert_half_bw

        def __set__(self, float vert_half_bw):
            self._Sweep.h.vert_half_bw = vert_half_bw

    property horz_half_bw:
        def __get__(self):
            return self._Sweep.h.horz_half_bw

        def __set__(self, float horz_half_bw):
            self._Sweep.h.horz_half_bw = horz_half_bw

    property nrays:
        def __get__(self):
            return self._Sweep.h.nrays

        def __set__(self, int nrays):
            self._Sweep.h.nrays = nrays


cdef class _RslVolume:
    """
    A object for accessing RSL Volume data and header information.

    This class should not be initalized from within Python.  _RslVolume
    objects are returned from the :py:func:`RslFile.get_volume` and other
    functions/methods.

    Attributes
    ----------
    nsweeps : int
        Sweep number.
    calibr_const : float
        Calibration constant.

    """

    cdef _rsl_interface.Volume * _Volume

    cdef load(self, _rsl_interface.Volume * Volume):
        """ Load the _RslVolume object, must be called after creation. """
        self._Volume = Volume

    def get_sweep(self, int sweep_number):
        """
        get_sweep(sweep_numer)

        Return a _RslSweep for a given sweep number.

        Parameters
        ----------
        sweep_number : int
            Sweep number to retrieve

        Returns
        -------
        sweep : _RslSweep
            _RslSweep object containing the requested sweep.

        """
        rslsweep = _RslSweep()
        rslsweep.load(self._Volume.sweep[sweep_number])
        return rslsweep

    def get_nray_list(self):
        """
        get_nray_list()

        Return a list of the number of rays for each sweep.
        """
        return [self._Volume.sweep[i].h.nrays for i in range(self.nsweeps)]

    def get_azimuth_and_elev_array(self):
        """
        get_azimuth_and_elev_array()

        Return azimuth and elevation array for each sweep and ray.
        """
        cdef int nrays = self._Volume.sweep[0].h.nrays
        cdef _rsl_interface.Sweep * sweep
        cdef _rsl_interface.Ray * ray
        azimuth = np.empty([self.nsweeps, nrays], dtype='float32')
        elev = np.empty([self.nsweeps, nrays], dtype='float32')

        for i in range(self.nsweeps):
            sweep = self._Volume.sweep[i]
            for j in range(nrays):
                ray = sweep.ray[j]
                azimuth[i, j] = ray.h.azimuth
                elev[i, j] = ray.h.elev
        return azimuth, elev

    def get_sweep_azimuths(self):
        """
        get_sweep_azimuths()

        Return azimuth array for each sweep.
        """
        azimuth = np.empty((self.nsweeps), dtype='float32')
        for i in range(self.nsweeps):
            azimuth[i] = self._Volume.sweep[i].h.azimuth
        return azimuth

    def get_sweep_elevs(self):
        """
        get_sweep_elevs()

        Return elevation array for each sweep.
        """
        elev = np.empty((self.nsweeps), dtype='float32')
        for i in range(self.nsweeps):
            elev[i] = self._Volume.sweep[i].h.elev
        return elev

    def get_instr_params(self):
        """
        get_instr_params()

        Return instrumental parameter for the volume.

        Returns
        -------
        pm_data : array, (nsweeps)
            Array of prt modes.
        nv_data : array, (nsweeps, nrays)
            Array of nyquist velocities.
        pr_data : array, (nsweeps, nrays)
            Array of pulse repetition frequency in Hz.
        ur_data : array, (nsweeps, nrays)
            Array of unambiguous ranges, in km.

        """
        cdef int nrays = self._Volume.sweep[0].h.nrays
        nyq_vel = self._Volume.sweep[0].ray[0].h.nyq_vel
        cdef _rsl_interface.Sweep * sweep
        cdef _rsl_interface.Ray * ray

        valid_nyq_vel = abs(nyq_vel) > 0.1
        pm_data = np.empty(self.nsweeps, dtype='|S24')
        nv_data = np.empty((self.nsweeps, nrays), dtype='float32')
        pr_data = np.empty((self.nsweeps, nrays), dtype='float32')
        ur_data = np.empty((self.nsweeps, nrays), dtype='float32')

        for i in range(self.nsweeps):
            sweep = self._Volume.sweep[i]
            for j in range(nrays):
                ray = sweep.ray[j]
                if j == 0:
                    pm_data[i] = self._prtmode(ray.h)

                if valid_nyq_vel:
                    nv_data[i, j] = ray.h.nyq_vel
                else:
                    nv_data[i, j] = ray.h.wavelength * ray.h.prf / 4.0

                pr_data[i, j] = 1. / ray.h.prf
                ur_data[i, j] = ray.h.unam_rng * 1000.0

        return pm_data, nv_data, pr_data, ur_data

    def get_data(self):
        """
        get_data()

        Return the three-dimensional data contained in the volume.
        """
        cdef _rsl_interface.Range raw
        cdef np.ndarray[np.float32_t, ndim = 3] data

        vol = self._Volume

        nsweeps = vol.h.nsweeps
        nrays = vol.sweep[0].h.nrays
        nbins = vol.sweep[0].ray[0].h.nbins

        shape = (nsweeps, nrays, nbins)
        data = np.zeros(shape, dtype='float32') + 1.31072000e+05
        for nsweep in range(nsweeps):
            sweep = vol.sweep[nsweep]
            for nray in range(nrays):
                ray = sweep.ray[nray]
                nbins = ray.h.nbins
                for nbin in range(nbins):
                    raw = ray.range[nbin]
                    data[nsweep, nray, nbin] = ray.h.f(raw)
        return data

    cdef _prtmode(self, _rsl_interface.Ray_header h):
        """ Return the prt mode of a given Ray header. """
        # TODO need to add additional logic here
        if h.prf2 != h.prf:
            mode = 'dual                    '
        else:
            mode = 'fixed                   '
        return mode

    # header properties mapped to class attributes.
    property nsweeps:
        def __get__(self):
            return self._Volume.h.nsweeps

        def __set__(self, int nsweeps):
            self._Volume.h.nsweeps = nsweeps

    property calibr_const:
        def __get__(self):
            return self._Volume.h.calibr_const

        def __set__(self, float calibr_const):
            self._Volume.h.calibr_const = calibr_const


cdef class RslFile:
    """
    RslFile(filename)

    A object for accessing Radar data and parameter using the RSL library.

    Parameters
    ----------
    filename : str
        Radar file to read.

    Attributes
    ----------
    month : int
        Date, month (1-12).
    day : int
        Date, day (1-31).
    year : int
        Date, year (eg. 1993).
    hour : int
        Time, hour (0-23).
    minute : int
        Time, minute (0-59).
    sec : float
        Time, second + fractions of second.
    nvolumes : int
        Number of volume slots in the file.
    number : int
        Arbitrary number for this radar site.
    latd, latm, lats : int
        Latitude degrees, minutes and seconds for the site.
    lond, lonm, lons : int
        Longitude degrees, minutes and seconds for the site.
    height : int
        Height of site in meters above sea level.
    spulse : int
        Length of short pulse in ns.
    lpulse : int
        Length of long pulse in ns.
    scan_mode : int
        Scan mode, 0 for PPI, 1 for RHI.
    vcp : int
        Volume coverage pattern, WSR-88D only.

    """
    cdef _rsl_interface.Radar * _Radar
    cdef _rsl_interface.Volume * _Volume
    cdef _rsl_interface.Sweep * _Sweep
    cdef _rsl_interface.Ray * _Ray
    cdef readonly int _first_volume_idx

    def __cinit__(self, filename):
        """ Initalize the _RslFile object. """
        self._Radar = _rsl_interface.RSL_anyformat_to_radar(filename)

    def get_volume(self, int volume_number):
        """
        get_volume(volume_number)

        Return a _RslVolume for a given volume number.

        Parameters
        ----------
        volume_number : int
            Volume number to retrieve

        Returns
        -------
        volume : _RslVolume
            _RslVolume object containing requested volume.

        """
        rslvolume = _RslVolume()
        rslvolume.load(self._Radar.v[volume_number])
        return rslvolume

    def available_moments(self):
        """
        available_moments()

        Return a list of available volume moments.
        """
        av = []
        for i in range(self._Radar.h.nvolumes):
            if self._Radar.v[i] is not NULL:
                av.append(i)
        return av

    def get_radar_header(self):
        """
        get_radar_headers()

        Return a dictionary of radar header parameters.
        """
        return self._Radar.h

    def get_volume_array(self, int volume_num):
        """
        get_volume_array(volume_number)

        Return the three-dimensional data contained in a given volume.

        Parameters
        ----------
        volume_number : int

        Returns
        -------
        volume : array (nsweep, nrays, nbins), float32
            Array containing  data for the given volume.

        """
        cdef _rsl_interface.Range raw
        cdef np.ndarray[np.float32_t, ndim = 3] data

        vol = self._Radar.v[volume_num]

        nsweeps = vol.h.nsweeps
        nrays = vol.sweep[0].h.nrays
        nbins = vol.sweep[0].ray[0].h.nbins

        shape = (nsweeps, nrays, nbins)
        data = np.zeros(shape, dtype='float32') + 1.31072000e+05
        for nsweep in range(nsweeps):
            sweep = vol.sweep[nsweep]
            for nray in range(nrays):
                ray = sweep.ray[nray]
                nbins = ray.h.nbins
                for nbin in range(nbins):
                    raw = ray.range[nbin]
                    data[nsweep, nray, nbin] = ray.h.f(raw)
        return data

    # header properties mapped to class attributes.
    property month:
        def __get__(self):
            return self._Radar.h.month

        def __set__(self, int month):
            self._Radar.h.month = month

    property day:
        def __get__(self):
            return self._Radar.h.day

        def __set__(self, int day):
            self._Radar.h.day = day

    property year:
        def __get__(self):
            return self._Radar.h.year

        def __set__(self, int year):
            self._Radar.h.year = year

    property hour:
        def __get__(self):
            return self._Radar.h.hour

        def __set__(self, int hour):
            self._Radar.h.hour = hour

    property minute:
        def __get__(self):
            return self._Radar.h.minute

        def __set__(self, int minute):
            self._Radar.h.minute = minute

    property sec:
        def __get__(self):
            return self._Radar.h.sec

        def __set__(self, float sec):
            self._Radar.h.sec = sec

    property nvolumes:
        def __get__(self):
            return self._Radar.h.nvolumes

        def __set__(self, int nvolumes):
            self._Radar.h.nvolumes = nvolumes

    property number:
        def __get__(self):
            return self._Radar.h.number

        def __set__(self, int number):
            self._Radar.h.number = number

    property latd:
        def __get__(self):
            return self._Radar.h.latd

        def __set__(self, int latd):
            self._Radar.h.latd = latd

    property latm:
        def __get__(self):
            return self._Radar.h.latm

        def __set__(self, int latm):
            self._Radar.h.latm = latm

    property lats:
        def __get__(self):
            return self._Radar.h.lats

        def __set__(self, int lats):
            self._Radar.h.lats = lats

    property lond:
        def __get__(self):
            return self._Radar.h.lond

        def __set__(self, int lond):
            self._Radar.h.lond = lond

    property lonm:
        def __get__(self):
            return self._Radar.h.lonm

        def __set__(self, int lonm):
            self._Radar.h.lonm = lonm

    property lons:
        def __get__(self):
            return self._Radar.h.lons

        def __set__(self, int lons):
            self._Radar.h.lons = lons

    property height:
        def __get__(self):
            return self._Radar.h.height

        def __set__(self, int height):
            self._Radar.h.height = height

    property spulse:
        def __get__(self):
            return self._Radar.h.spulse

        def __set__(self, int spulse):
            self._Radar.h.spulse = spulse

    property lpulse:
        def __get__(self):
            return self._Radar.h.lpulse

        def __set__(self, int lpulse):
            self._Radar.h.lpulse = lpulse

    property scan_mode:
        def __get__(self):
            return self._Radar.h.scan_mode

        def __set__(self, int scan_mode):
            self._Radar.h.scan_mode = scan_mode

    property vcp:
        def __get__(self):
            return self._Radar.h.vcp

        def __set__(self, int vcp):
            self._Radar.h.vcp = vcp
