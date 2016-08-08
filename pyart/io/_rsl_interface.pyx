"""
pyart.io._rsl_interface
=======================

Cython wrapper around the NASA TRMM RSL library.

.. autosummary::
    :toctree: generated/

    copy_volume
    create_volume
    _label_volume

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    RslFile
    _RslVolume
    _RslSweep
    _RslRay


"""

cimport _rsl_h
# the next line is required so that RSL_F_LIST and RSL_INVF_LIST can be
# properly wrapped as Cython does not export the typedef from _rsl_h
ctypedef unsigned short Range
import numpy as np
cimport numpy as np
from datetime import datetime, timedelta


_RSL_VERSION_STR = _rsl_h._RSL_VERSION_STR


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
    volume_copy = _rsl_h.RSL_copy_volume(volume._Volume)
    rslvolume = _RslVolume()
    rslvolume.load(volume_copy)
    rslvolume._dealloc = 1
    return rslvolume


cpdef create_volume(
        np.ndarray[np.float32_t, ndim=2] arr,
        np.ndarray[np.int32_t, ndim=1] rays_per_sweep,
        int vol_num=1):
    """
    create_volume(arr, rays_per_sweep, vol_num=1)

    Create a _RslVolume object from a 2D float32 array.

    No headers parameters except nsweeps, nrays and nbins are not set in the
    resulting _RslVolume object.

    Parameters
    ----------
    arr : array, 2D, float32
        Two dimensional float32 array.
    rays_per_sweep: array, 1D, int32
        Array listing number of rays in each sweep.
    vol_num : int
        Volume number used to set f and invf in the header.  The default is
        for velocity fields.  Useful values are 0 for reflectivity and
        1 for velocity.

    Returns
    -------
    volumes : _RslVolume
        _RslVolume containing array data.

    """
    # these variables can be moved to the module level if used elsewhere
    cdef (float (*)(_rsl_h.Range)) * RSL_F_LIST = [
        _rsl_h.DZ_F, _rsl_h.VR_F, _rsl_h.SW_F, _rsl_h.CZ_F, _rsl_h.ZT_F,
        _rsl_h.DR_F, _rsl_h.LR_F, _rsl_h.ZD_F, _rsl_h.DM_F, _rsl_h.RH_F,
        _rsl_h.PH_F, _rsl_h.XZ_F, _rsl_h.CD_F, _rsl_h.MZ_F, _rsl_h.MD_F,
        _rsl_h.ZE_F, _rsl_h.VE_F, _rsl_h.KD_F, _rsl_h.TI_F, _rsl_h.DX_F,
        _rsl_h.CH_F, _rsl_h.AH_F, _rsl_h.CV_F, _rsl_h.AV_F, _rsl_h.SQ_F,
        _rsl_h.VS_F, _rsl_h.VL_F, _rsl_h.VG_F, _rsl_h.VT_F, _rsl_h.NP_F,
        _rsl_h.HC_F, _rsl_h.VC_F, _rsl_h.VR_F, _rsl_h.SW_F, _rsl_h.VR_F,
        _rsl_h.SW_F, _rsl_h.DZ_F, _rsl_h.CZ_F, _rsl_h.PH_F, _rsl_h.SD_F,
        _rsl_h.DZ_F, _rsl_h.DZ_F]

    cdef (_rsl_h.Range (*)(float)) * RSL_INVF_LIST = [
        _rsl_h.DZ_INVF, _rsl_h.VR_INVF, _rsl_h.SW_INVF, _rsl_h.CZ_INVF,
        _rsl_h.ZT_INVF, _rsl_h.DR_INVF, _rsl_h.LR_INVF, _rsl_h.ZD_INVF,
        _rsl_h.DM_INVF, _rsl_h.RH_INVF, _rsl_h.PH_INVF, _rsl_h.XZ_INVF,
        _rsl_h.CD_INVF, _rsl_h.MZ_INVF, _rsl_h.MD_INVF, _rsl_h.ZE_INVF,
        _rsl_h.VE_INVF, _rsl_h.KD_INVF, _rsl_h.TI_INVF, _rsl_h.DX_INVF,
        _rsl_h.CH_INVF, _rsl_h.AH_INVF, _rsl_h.CV_INVF, _rsl_h.AV_INVF,
        _rsl_h.SQ_INVF, _rsl_h.VS_INVF, _rsl_h.VL_INVF, _rsl_h.VG_INVF,
        _rsl_h.VT_INVF, _rsl_h.NP_INVF, _rsl_h.HC_INVF, _rsl_h.VC_INVF,
        _rsl_h.VR_INVF, _rsl_h.SW_INVF, _rsl_h.VR_INVF, _rsl_h.SW_INVF,
        _rsl_h.DZ_INVF, _rsl_h.CZ_INVF, _rsl_h.PH_INVF, _rsl_h.SD_INVF,
        _rsl_h.DZ_INVF, _rsl_h.DZ_INVF]

    cdef _rsl_h.Volume * volume
    cdef _rsl_h.Sweep * sweep
    cdef _rsl_h.Ray * ray
    cdef int ray_i

    nsweeps = rays_per_sweep.shape[0]
    nbins = arr.shape[1]

    volume = _rsl_h.RSL_new_volume(nsweeps)
    volume.h.nsweeps = nsweeps

    ray_index = 0
    for nsweep in range(nsweeps):
        nrays = rays_per_sweep[nsweep]
        sweep = _rsl_h.RSL_new_sweep(nrays)
        volume.sweep[nsweep] = sweep
        sweep.h.nrays = nrays
        for nray in range(nrays):
            ray = _rsl_h.RSL_new_ray(nbins)
            sweep.ray[nray] = ray
            ray.h.nbins = nbins
            ray.h.f = RSL_F_LIST[vol_num]
            ray.h.invf = RSL_INVF_LIST[vol_num]
            for nbin in range(nbins):
                ray.range[nbin] = ray.h.invf(arr[ray_index, nbin])
            ray_index += 1

    rslvolume = _RslVolume()
    rslvolume.load(volume)
    rslvolume._dealloc = 1
    return rslvolume


cpdef _label_volume(_RslVolume volume, radar):
    """
    _label_volume(volume, radar)

    Add labels for dealiasing to a _RslVolume object from a radar object.

    This function does not set all parameter in the _RslVolume suitable for
    writing out the volume, rather it set those parameters which must be set
    prior to using :py:func:`pyart.correct._fourdd_interface.fourdd_dealias`.

    Parameters
    ----------
    volume : _RslVolume
        Volume object to which parameters will be set as needed prior to
        dealiasing.  Object is manipulated in-place.
    radar : Radar
        Radar object from which parameters are taken.

    """

    cdef _rsl_h.Sweep * sweep
    cdef _rsl_h.Ray * ray
    cdef int ray_index

    vol = volume._Volume
    nsweeps = vol.h.nsweeps
    nbins = vol.sweep[0].ray[0].h.nbins

    gate_size = int(radar.range['meters_between_gates'])
    range_bin1 = int(radar.range['meters_to_center_of_first_gate'])
    if 'shape' in dir(radar.altitude['data']):
        if radar.altitude['data'].shape == ():
            alt = float(radar.altitude['data'])
        else:
            alt = radar.altitude['data'][0]
    else:
        alt = radar.altitude['data']

    nyq_vels = radar.instrument_parameters['nyquist_velocity']['data']
    azimuths = radar.azimuth['data']
    elevs = radar.elevation['data']

    # label the volume
    ray_index = 0
    vol.h.nsweeps = nsweeps
    for nsweep in range(nsweeps):
        sweep = vol.sweep[nsweep]
        nrays = sweep.h.nrays
        for nray in range(nrays):
            ray = sweep.ray[nray]
            ray.h.azimuth = azimuths[ray_index]
            ray.h.elev = elevs[ray_index]
            ray.h.nyq_vel = nyq_vels[ray_index]
            ray.h.range_bin1 = range_bin1
            ray.h.gate_size = gate_size
            ray.h.alt = alt
            ray_index += 1
    return


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

    cdef _rsl_h.Ray * _Ray
    cdef int _dealloc

    def __dealloc__(self):
        if self._dealloc == 1:
            _rsl_h.RSL_free_ray(self._Ray)

    cdef load(self, _rsl_h.Ray * Ray):
        """ Load the _RslRay object, must be called after creation. """
        if Ray is NULL:
            raise ValueError('cannot load _RslRay with NULL')
        self._Ray = Ray
        self._dealloc = 0

    def get_datetime(self):
        """
        get_datetime()

        Return a datetime describing the date and time of the ray.
        """
        s = self
        full_seconds, fractional_seconds = divmod(s.sec, 1)
        microseconds = int(fractional_seconds * 1e6)
        # Some UF writers incorrectly specify midnight as 24:00:00 rather
        # than 00:00:00.  Handle this case explicitly
        if s.hour == 24:
            s.hour = 23
            s.minute = 0
            full_seconds = 0
            dt = datetime(s.year, s.month, s.day, s.hour, s.minute,
                          full_seconds, microseconds)
            return dt + timedelta(seconds=1)
        return datetime(s.year, s.month, s.day, s.hour, s.minute,
                        int(full_seconds), microseconds)

    def get_data(self):
        """
        get_data()

        Return the one-dimensional data contained in the ray.
        """
        cdef _rsl_h.Range raw
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

    cdef _rsl_h.Sweep * _Sweep
    cdef int _dealloc

    def __dealloc__(self):
        if self._dealloc == 1:
            _rsl_h.RSL_free_sweep(self._Sweep)

    cdef load(self, _rsl_h.Sweep * Sweep):
        """ Load the _RslSweep object, must be called after creation. """
        if Sweep is NULL:
            raise ValueError("cannot load _RslSweep with NULL")
        self._Sweep = Sweep
        self._dealloc = 0

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
        if ray_number < 0 or ray_number >= self._Sweep.h.nrays:
            raise ValueError('invalid ray_number')
        rslray = _RslRay()
        rslray.load(self._Sweep.ray[ray_number])
        return rslray

    def get_data(self):
        """
        get_data()

        Return the two-dimensional data contained in the sweep.

        If a given ray has few bins than the first ray, the missing bins
        will be filled with 131072.0
        """
        cdef _rsl_h.Range raw
        cdef _rsl_h.Ray * ray
        cdef np.ndarray[np.float32_t, ndim = 2] data

        sweep = self._Sweep
        nrays = sweep.h.nrays
        nbins = sweep.ray[0].h.nbins

        shape = (nrays, nbins)
        data = np.zeros(shape, dtype='float32') + 1.31072000e+05
        for nray in range(nrays):
            ray = sweep.ray[nray]
            assert ray is not NULL
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

    def __dealloc__(self):
        if self._dealloc == 1:
            _rsl_h.RSL_free_volume(self._Volume)

    cdef load(self, _rsl_h.Volume * Volume):
        """ Load the _RslVolume object, must be called after creation. """
        if Volume is NULL:
            raise ValueError('cannot load _RslVolume with NULL')
        self._Volume = Volume
        self._dealloc = 0

    def total_rays(self):
        """
        total_rays()

        Return the total number of rays present in all sweeps of the volume.
        """
        return np.sum(self.get_nray_array())

    def get_nray_array(self):
        """
        get_nray_array()

        Return an array of the number of rays for each sweep.
        """
        cdef _rsl_h.Sweep * sweep
        nrays = np.empty((self.nsweeps), dtype='int32')
        for i in range(self.nsweeps):
            sweep = self._Volume.sweep[i]
            assert sweep is not NULL
            nrays[i] = sweep.h.nrays
        return nrays

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
        if sweep_number < 0 or sweep_number >= self._Volume.h.nsweeps:
            raise ValueError('invalid sweep_number')

        rslsweep = _RslSweep()
        rslsweep.load(self._Volume.sweep[sweep_number])
        return rslsweep

    def get_azimuth_and_elev_array(self):
        """
        get_azimuth_and_elev_array()

        Return azimuth and elevation array for each sweep and ray.
        """
        cdef int nrays = self._Volume.sweep[0].h.nrays
        cdef int ray_count
        cdef _rsl_h.Sweep * sweep
        cdef _rsl_h.Ray * ray

        # create empty azimuth and elev output arrays
        total_rays = self.total_rays()
        azimuth = np.empty([total_rays], dtype='float32')
        elev = np.empty([total_rays], dtype='float32')

        # loop over the sweeps and rays storing azimuth and elev
        ray_count = 0
        for i in range(self.nsweeps):
            sweep = self._Volume.sweep[i]
            assert sweep is not NULL
            nrays = sweep.h.nrays
            for j in range(nrays):
                ray = sweep.ray[j]
                assert ray is not NULL
                azimuth[ray_count + j] = ray.h.azimuth
                elev[ray_count + j] = ray.h.elev
            ray_count += nrays
        return azimuth, elev

    def get_sweep_fix_angles(self):
        """
        get_sweep_fix_angles()

        Return array of fix angle for each sweep.

        Angles determined from the first ray in each sweep.
        """
        cdef _rsl_h.Sweep * sweep
        cdef _rsl_h.Ray * ray
        fix_angles = np.empty((self.nsweeps), dtype='float32')
        for i in range(self.nsweeps):
            sweep = self._Volume.sweep[i]
            assert sweep is not NULL
            ray = sweep.ray[0]
            assert ray is not NULL
            fix_angles[i] = ray.h.fix_angle
        return fix_angles

    def get_sweep_azimuths(self):
        """
        get_sweep_azimuths()

        Return azimuth array for each sweep.
        """
        cdef _rsl_h.Sweep * sweep
        azimuth = np.empty((self.nsweeps), dtype='float32')
        for i in range(self.nsweeps):
            sweep = self._Volume.sweep[i]
            assert sweep is not NULL
            azimuth[i] = sweep.h.azimuth
        return azimuth

    def get_sweep_elevs(self):
        """
        get_sweep_elevs()

        Return elevation array for each sweep.
        """
        cdef _rsl_h.Sweep * sweep
        elev = np.empty((self.nsweeps), dtype='float32')
        for i in range(self.nsweeps):
            sweep = self._Volume.sweep[i]
            assert sweep is not NULL
            elev[i] = sweep.h.elev
        return elev

    def get_instr_params(self):
        """
        get_instr_params()

        Return instrumental parameter for the volume.

        Returns
        -------
        pm_data : array, (nsweeps)
            Array of prt modes.
        nv_data : array, (total_rays)
            Array of nyquist velocities.
        pr_data : array, (total_rays)
            Array of pulse repetition frequency in Hz.
        ur_data : array, (total_rays)
            Array of unambiguous ranges, in km.

        """
        cdef int nrays = self._Volume.sweep[0].h.nrays
        cdef int ray_count
        cdef _rsl_h.Sweep * sweep
        cdef _rsl_h.Ray * ray

        # calculate the total number of rays in the volume

        # initalize empty instrument parameter arrays
        total_rays = self.total_rays()
        nyq_vel = self._Volume.sweep[0].ray[0].h.nyq_vel
        valid_nyq_vel = abs(nyq_vel) > 0.1
        pm_data = np.empty(self.nsweeps, dtype='|S24')
        nv_data = np.empty((total_rays), dtype='float32')
        pr_data = np.empty((total_rays), dtype='float32')
        ur_data = np.empty((total_rays), dtype='float32')

        # loop over sweeps and rays storing instrument parameters
        ray_count = 0
        for i in range(self.nsweeps):
            sweep = self._Volume.sweep[i]
            assert sweep is not NULL
            nrays = sweep.h.nrays
            for j in range(nrays):
                ray = sweep.ray[j]
                assert ray is not NULL
                if j == 0:
                    pm_data[i] = self._prtmode(ray.h)

                if valid_nyq_vel:
                    nv_data[ray_count + j] = ray.h.nyq_vel
                else:
                    nv_data[ray_count + j] = (ray.h.wavelength *
                                              ray.h.prf / 4.0)

                if ray.h.prf == 0:
                    pr_data[ray_count + j] = -999.
                else:
                    pr_data[ray_count + j] = 1. / ray.h.prf
                ur_data[ray_count + j] = ray.h.unam_rng * 1000.0
            ray_count += nrays
        return pm_data, nv_data, pr_data, ur_data

    def get_data(self):
        """
        get_data()

        Return the two-dimensional data contained in the volume.

        If a given ray has few bins than the first ray, the missing bins
        will be filled with 131072.0
        """
        cdef _rsl_h.Range raw
        cdef _rsl_h.Sweep * sweep
        cdef _rsl_h.Ray * ray
        cdef int ray_count, nsweeps, nrays, nbins, nray, nsweep, nbin
        cdef np.ndarray[np.float32_t, ndim = 2] data

        vol = self._Volume
        nbins = vol.sweep[0].ray[0].h.nbins
        total_rays = self.total_rays()
        shape = (total_rays, nbins)
        data = np.zeros(shape, dtype='float32') + 1.31072000e+05

        ray_count = 0
        nsweeps = vol.h.nsweeps
        for nsweep in range(nsweeps):
            sweep = vol.sweep[nsweep]
            assert sweep is not NULL
            nrays = sweep.h.nrays
            for nray in range(nrays):
                ray = sweep.ray[nray]
                assert ray is not NULL
                nbins = ray.h.nbins
                for nbin in range(nbins):
                    raw = ray.range[nbin]
                    data[ray_count + nray, nbin] = ray.h.f(raw)
            ray_count += nrays
        return data

    def is_range_bins_uniform(self):
        """
        is_range_bins_uniform()

        Return True is the locations of the range bin are identical for all
        rays, False if locations change in one or more rays.
        """
        cdef int nrays = self._Volume.sweep[0].h.nrays
        cdef _rsl_h.Sweep * sweep
        cdef _rsl_h.Ray * ray

        # loop over the sweeps and rays checking that the gate_size and
        # range_bin1 are the same as the that in the first ray
        sweep = self._Volume.sweep[0]
        assert sweep is not NULL
        ray = sweep.ray[0]
        assert ray is not NULL
        ref_gate_size = ray.h.gate_size
        ref_range_bin1 = ray.h.range_bin1

        for i in range(self.nsweeps):
            sweep = self._Volume.sweep[i]
            assert sweep is not NULL
            nrays = sweep.h.nrays
            for j in range(nrays):
                ray = sweep.ray[j]
                assert ray is not NULL
                if ray.h.gate_size != ref_gate_size:
                    return False
                if ray.h.range_bin1 != ref_range_bin1:
                    return False
        return True

    cdef _prtmode(self, _rsl_h.Ray_header h):
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
    cdef _rsl_h.Radar * _Radar
    cdef _rsl_h.Volume * _Volume
    cdef _rsl_h.Sweep * _Sweep
    cdef _rsl_h.Ray * _Ray

    def __cinit__(self, filename, radar_format=None, callid=None):
        """ Initalize the _RslFile object. """
        if radar_format == 'wsr88d':
            if callid is None:
                raise ValueError('callid must be provided.')
            self._Radar = _rsl_h.RSL_wsr88d_to_radar(filename, callid)
        elif radar_format is None:
            self._Radar = _rsl_h.RSL_anyformat_to_radar(filename)
        else:
            raise ValueError('invalid radar_format:', radar_format)
        if self._Radar is NULL:
            raise IOError('file cannot be read.')

    def __dealloc__(self):
        """ Free memory used by object. """
        _rsl_h.RSL_free_radar(self._Radar)

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
        if volume_number < 0 or volume_number >= self._Radar.h.nvolumes:
            raise ValueError('invalid volume_number')

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
        return self.get_volume(volume_num).get_data()

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
