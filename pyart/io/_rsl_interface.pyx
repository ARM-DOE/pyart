""" Cython wrapper around NASA TRMM RSL Library """

cimport _rsl_interface
import numpy as np
cimport numpy as np
from datetime import datetime

from pyart.io.common import get_metadata


cdef class _RslRay:
    cdef _rsl_interface.Ray * _Ray

    cdef load(self, _rsl_interface.Ray * Ray):
        self._Ray = Ray

    def get_datetime(self):
        s = self
        return datetime(s.year, s.month, s.day, s.hour, s.minute, int(s.sec))

    # header properties
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

    cdef _rsl_interface.Sweep * _Sweep

    cdef load(self, _rsl_interface.Sweep * Sweep):
        self._Sweep = Sweep

    def get_ray(self, int ray_number):
        rslray = _RslRay()
        rslray.load(self._Sweep.ray[ray_number])
        return rslray

    # header properties
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

    cdef _rsl_interface.Volume * _Volume

    cdef load(self, _rsl_interface.Volume * Volume):
        self._Volume = Volume

    def get_sweep(self, int sweep_number):
        rslsweep = _RslSweep()
        rslsweep.load(self._Volume.sweep[sweep_number])
        return rslsweep

    def get_nray_list(self):
        return [self._Volume.sweep[i].h.nrays for i in range(self.nsweeps)]

    def get_azimuth_and_elev_array(self):

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
        azimuth = np.empty((self.nsweeps), dtype='float32')
        for i in range(self.nsweeps):
            azimuth[i] = self._Volume.sweep[i].h.azimuth
        return azimuth

    def get_sweep_elevs(self):
        elev = np.empty((self.nsweeps), dtype='float32')
        for i in range(self.nsweeps):
            elev[i] = self._Volume.sweep[i].h.elev
        return elev

    def get_instr_params(self):
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

    cdef _prtmode(self, _rsl_interface.Ray_header h):
        # TODO need to add additional logic here
        if h.prf2 != h.prf:
            mode = 'dual                    '
        else:
            mode = 'fixed                   '
        return mode

    # header attributes

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

    cdef _rsl_interface.Radar * _Radar
    cdef _rsl_interface.Volume * _Volume
    cdef _rsl_interface.Sweep * _Sweep
    cdef _rsl_interface.Ray * _Ray
    cdef readonly int _first_volume_idx

    def __cinit__(self, s):
        self._Radar = _rsl_interface.RSL_anyformat_to_radar(s)

    def get_volume(self, int volume_number):
        rslvolume = _RslVolume()
        rslvolume.load(self._Radar.v[volume_number])
        return rslvolume

    def available_moments(self):
        av = []
        for i in range(self._Radar.h.nvolumes):
            if self._Radar.v[i] is not NULL:
                av.append(i)
        return av

    def get_radar_header(self):
        return self._Radar.h

    def get_volume_array(self, int volume_num):

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
