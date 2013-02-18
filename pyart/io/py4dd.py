""" Python Wrapper for the University of Washington 4D dealias code

TODO
----
split into fourdd.py, _rslwrap.py, and move others to rsl_utils.py

"""

from ctypes import *
import numpy as N
import os.path
import datetime
import sys


def find_librsl_path():
    """ Return the path of the RSL dynamic library. """
    # assumed to be ./lib/librsl.{so, dylib, DLL}
    ext = {'darwin': 'dylib', 'linux2': 'so', 'win32': 'DLL'}[sys.platform]
    lib_dir = os.path.join(os.path.dirname(__file__), 'lib')
    return os.path.join(lib_dir, 'librsl.' + ext)

librsl = CDLL(find_librsl_path())

STRING = c_char_p
Range = c_ushort


class fieldTypes:

    def __init__(self):
        self.DZ = 0
        self.VR = 1
        self.SW = 2
        self.CZ = 3
        self.ZT = 4
        self.DR = 5
        self.LR = 6
        self.ZD = 7
        self.DM = 8
        self.RH = 9
        self.PH = 10
        self.XZ = 11
        self.CD = 12
        self.MZ = 13
        self.MD = 14
        self.ZE = 15
        self.VE = 16
        self.KD = 17
        self.TI = 18
        self.DX = 19
        self.CH = 20
        self.AH = 21
        self.CV = 22
        self.AV = 23
        self.SQ = 24

        self.list = ['DZ', 'VR', 'SW', 'CZ', 'ZT', 'DR', 'LR', 'ZD', 'DM',
                     'RH', 'PH', 'XZ', 'CD', 'MZ', 'MD', 'ZE', 'VE', 'KD',
                     'TI', 'DX', 'CH', 'AH', 'CV', 'AV', 'SQ', ]


# librsl return and arguments types
DZ_F = librsl.DZ_F
DZ_F.restype = c_float
DZ_F.argtypes = [Range]

VR_F = librsl.VR_F
VR_F.restype = c_float
VR_F.argtypes = [Range]

SW_F = librsl.SW_F
SW_F.restype = c_float
SW_F.argtypes = [Range]

CZ_F = librsl.CZ_F
CZ_F.restype = c_float
CZ_F.argtypes = [Range]

ZT_F = librsl.ZT_F
ZT_F.restype = c_float
ZT_F.argtypes = [Range]

DR_F = librsl.DR_F
DR_F.restype = c_float
DR_F.argtypes = [Range]

LR_F = librsl.LR_F
LR_F.restype = c_float
LR_F.argtypes = [Range]

ZD_F = librsl.ZD_F
ZD_F.restype = c_float
ZD_F.argtypes = [Range]

DM_F = librsl.DM_F
DM_F.restype = c_float
DM_F.argtypes = [Range]

RH_F = librsl.RH_F
RH_F.restype = c_float
RH_F.argtypes = [Range]

PH_F = librsl.PH_F
PH_F.restype = c_float
PH_F.argtypes = [Range]

XZ_F = librsl.XZ_F
XZ_F.restype = c_float
XZ_F.argtypes = [Range]

CD_F = librsl.CD_F
CD_F.restype = c_float
CD_F.argtypes = [Range]

MZ_F = librsl.MZ_F
MZ_F.restype = c_float
MZ_F.argtypes = [Range]

MD_F = librsl.MD_F
MD_F.restype = c_float
MD_F.argtypes = [Range]

ZE_F = librsl.ZE_F
ZE_F.restype = c_float
ZE_F.argtypes = [Range]

VE_F = librsl.VE_F
VE_F.restype = c_float
VE_F.argtypes = [Range]

KD_F = librsl.KD_F
KD_F.restype = c_float
KD_F.argtypes = [Range]

TI_F = librsl.TI_F
TI_F.restype = c_float
TI_F.argtypes = [Range]

DX_F = librsl.DX_F
DX_F.restype = c_float
DX_F.argtypes = [Range]

CH_F = librsl.CH_F
CH_F.restype = c_float
CH_F.argtypes = [Range]

AH_F = librsl.AH_F
AH_F.restype = c_float
AH_F.argtypes = [Range]

CV_F = librsl.CV_F
CV_F.restype = c_float
CV_F.argtypes = [Range]

AV_F = librsl.AV_F
AV_F.restype = c_float
AV_F.argtypes = [Range]

SQ_F = librsl.SQ_F
SQ_F.restype = c_float
SQ_F.argtypes = [Range]

DZ_INVF = librsl.DZ_INVF
DZ_INVF.restype = Range
DZ_INVF.argtypes = [c_float]

VR_INVF = librsl.VR_INVF
VR_INVF.restype = Range
VR_INVF.argtypes = [c_float]

SW_INVF = librsl.SW_INVF
SW_INVF.restype = Range
SW_INVF.argtypes = [c_float]

CZ_INVF = librsl.CZ_INVF
CZ_INVF.restype = Range
CZ_INVF.argtypes = [c_float]

ZT_INVF = librsl.ZT_INVF
ZT_INVF.restype = Range
ZT_INVF.argtypes = [c_float]

DR_INVF = librsl.DR_INVF
DR_INVF.restype = Range
DR_INVF.argtypes = [c_float]

LR_INVF = librsl.LR_INVF
LR_INVF.restype = Range
LR_INVF.argtypes = [c_float]

ZD_INVF = librsl.ZD_INVF
ZD_INVF.restype = Range
ZD_INVF.argtypes = [c_float]

DM_INVF = librsl.DM_INVF
DM_INVF.restype = Range
DM_INVF.argtypes = [c_float]

RH_INVF = librsl.RH_INVF
RH_INVF.restype = Range
RH_INVF.argtypes = [c_float]

PH_INVF = librsl.PH_INVF
PH_INVF.restype = Range
PH_INVF.argtypes = [c_float]

XZ_INVF = librsl.XZ_INVF
XZ_INVF.restype = Range
XZ_INVF.argtypes = [c_float]

CD_INVF = librsl.CD_INVF
CD_INVF.restype = Range
CD_INVF.argtypes = [c_float]

MZ_INVF = librsl.MZ_INVF
MZ_INVF.restype = Range
MZ_INVF.argtypes = [c_float]

MD_INVF = librsl.MD_INVF
MD_INVF.restype = Range
MD_INVF.argtypes = [c_float]

ZE_INVF = librsl.ZE_INVF
ZE_INVF.restype = Range
ZE_INVF.argtypes = [c_float]

VE_INVF = librsl.VE_INVF
VE_INVF.restype = Range
VE_INVF.argtypes = [c_float]

KD_INVF = librsl.KD_INVF
KD_INVF.restype = Range
KD_INVF.argtypes = [c_float]

TI_INVF = librsl.TI_INVF
TI_INVF.restype = Range
TI_INVF.argtypes = [c_float]

DX_INVF = librsl.DX_INVF
DX_INVF.restype = Range
DX_INVF.argtypes = [c_float]

CH_INVF = librsl.CH_INVF
CH_INVF.restype = Range
CH_INVF.argtypes = [c_float]

AH_INVF = librsl.AH_INVF
AH_INVF.restype = Range
AH_INVF.argtypes = [c_float]

CV_INVF = librsl.CV_INVF
CV_INVF.restype = Range
CV_INVF.argtypes = [c_float]

AV_INVF = librsl.AV_INVF
AV_INVF.restype = Range
AV_INVF.argtypes = [c_float]

SQ_INVF = librsl.SQ_INVF
SQ_INVF.restype = Range
SQ_INVF.argtypes = [c_float]


class Ray_header(Structure):
    _fields_ = [('month', c_int),
                ('day', c_int),
                ('year', c_int),
                ('hour', c_int),
                ('minute', c_int),
                ('sec', c_float),
                ('unam_rng', c_float),
                ('azimuth', c_float),
                ('ray_num', c_int),
                ('elev', c_float),
                ('elev_num', c_int),
                ('range_bin1', c_int),
                ('gate_size', c_int),
                ('vel_res', c_float),
                ('sweep_rate', c_float),
                ('prf', c_int),
                ('prf2', c_int),
                ('azim_rate', c_float),
                ('fix_angle', c_float),
                ('pitch', c_float),
                ('roll', c_float),
                ('heading', c_float),
                ('pitch_rate', c_float),
                ('roll_rate', c_float),
                ('heading_rate', c_float),
                ('lat', c_float),
                ('lon', c_float),
                ('alt', c_int),
                ('rvc', c_float),
                ('vel_east', c_float),
                ('vel_north', c_float),
                ('vel_up', c_float),
                ('pulse_count', c_int),
                ('pulse_width', c_float),
                ('beam_width', c_float),
                ('frequency', c_float),
                ('wavelength', c_float),
                ('nyq_vel', c_float),
                ('f', CFUNCTYPE(c_float, c_ushort)),
                ('invf', CFUNCTYPE(Range, c_float)),
                ('nbins', c_int), ]


class Ray(Structure):
    _fields_ = [('h', Ray_header),
                ('range', POINTER(Range)), ]

    def __getattr__(self, attr):
        if attr == 'data':
            nbins = self.h.nbins
            converter = self.h.f
            data = N.array([converter(self.range[i]) for i in range(nbins)])
            return data
        elif attr == 'dists':
            #return range gates in km
            dists = self.h.range_bin1 + self.h.gate_size * N.arange(
                self.h.nbins, dtype=float)
            return dists
        else:
            return self.h.__getattribute__(attr)


class Sweep_header(Structure):
    _fields_ = [('sweep_num', c_int),
                ('elev', c_float),
                ('azimuth', c_float),
                ('beam_width', c_float),
                ('vert_half_bw', c_float),
                ('horz_half_bw', c_float),
                ('nrays', c_int),
                ('f', CFUNCTYPE(c_float, c_ushort)),
                ('invf', CFUNCTYPE(Range, c_float)), ]


class Sweep(Structure):
    _fields_ = [('h', Sweep_header),
                ('ray', POINTER(POINTER(Ray))), ]

    def __getattr__(self, attr):
        if attr == 'rays':
            nrays = self.h.nrays
            rays = []
            for i in range(nrays):
                # There are null pointers where rays don't exist, which
                # throws an error on access
                try:
                    rays.append(self.ray[i].contents)
                except:
                    rays.append(None)
            return rays
        else:
            return self.h.__getattribute__(attr)


class Volume_header(Structure):
    _fields_ = [('type_str', STRING),
                ('nsweeps', c_int),
                ('calibr_const', c_float),
                ('f', CFUNCTYPE(c_float, c_ushort)),
                ('invf', CFUNCTYPE(Range, c_float)), ]


class Volume(Structure):
    _fields_ = [('h', Volume_header),
                ('sweep', POINTER(POINTER(Sweep))), ]

    def __getattr__(self, attr):
        if attr == 'sweeps':
            nsweeps = self.h.nsweeps
            swps = []
            for i in range(nsweeps):
                # There are null pointers where sweeps don't exist,
                # which throws an error on access
                try:
                    swps.append(self.sweep[i].contents)
                except:
                    swps.append(None)
            return swps
        else:
            return self.h.__getattribute__(attr)


class Radar_header(Structure):
    _fields_ = [('month', c_int),
                ('day', c_int),
                ('year', c_int),
                ('hour', c_int),
                ('minute', c_int),
                ('sec', c_float),
                ('radar_type', c_char * 50),
                ('nvolumes', c_int),
                ('number', c_int),
                ('name', c_char * 8),
                ('radar_name', c_char * 8),
                ('project', c_char * 24),
                ('city', c_char * 15),
                ('state', c_char * 3),
                ('country', c_char * 15),
                ('latd', c_int),
                ('latm', c_int),
                ('lats', c_int),
                ('lond', c_int),
                ('lonm', c_int),
                ('lons', c_int),
                ('height', c_int),
                ('spulse', c_int),
                ('lpulse', c_int),
                ('scanmode', c_int),
                ('vcp', c_int), ]


class Radar(Structure):
    _fields_ = [('h', Radar_header),
                ('v', POINTER(POINTER(Volume)))]

    def __getattr__(self, attr):
        if attr == 'volumes':
            nvolumes = self.h.nvolumes
            vols = []
            for i in range(nvolumes):
                # There are null pointers where volumes don't exist,
                # which throws an error on access
                try:
                    vols.append(self.v[i].contents)
                except:
                    vols.append(None)
                #vols = [self.v[i].contents for i in range(nvolumes)]
            return vols
        elif attr == 'datetime':
            datetimeobj = datetime.datetime(
                self.h.year, self.h.month, self.h.day, self.h.hour,
                self.h.minute, int(self.h.sec))
            return datetimeobj
        else:
            return self.h.__getattribute__(attr)


RSL_anyformat_to_radar = librsl.RSL_anyformat_to_radar
RSL_anyformat_to_radar.restype = POINTER(Radar)
RSL_anyformat_to_radar.argtypes = [STRING]

RSL_wsr88d_to_radar = librsl.RSL_wsr88d_to_radar
RSL_wsr88d_to_radar.restype = POINTER(Radar)
RSL_wsr88d_to_radar.argtypes = [STRING, STRING]


RSL_radar_to_uf = librsl.RSL_radar_to_uf
RSL_radar_to_uf.argtypes = [POINTER(Radar), STRING]

RSL_radar_verbose_on = librsl.RSL_radar_verbose_on
RSL_radar_verbose_on.restype = None

RSL_free_radar = librsl.RSL_free_radar
RSL_free_radar.restype = None
RSL_free_radar.argtypes = [POINTER(Radar)]

RSL_new_radar = librsl.RSL_new_radar
RSL_new_radar.restype = POINTER(Radar)
RSL_new_radar.argtypes = [c_int]

RSL_new_volume = librsl.RSL_new_volume
RSL_new_volume.restype = POINTER(Volume)
RSL_new_volume.argtypes = [c_int]

RSL_new_sweep = librsl.RSL_new_sweep
RSL_new_sweep.restype = POINTER(Sweep)
RSL_new_sweep.argtypes = [c_int]

RSL_new_ray = librsl.RSL_new_ray
RSL_new_ray.restype = POINTER(Ray)
RSL_new_ray.argtypes = [c_int]

RSL_copy_volume = librsl.RSL_copy_volume
RSL_copy_volume.restype = POINTER(Volume)
RSL_copy_volume.argtypes = [POINTER(Volume)]

firstGuess = librsl.firstGuess
firstGuess.argtypes = [POINTER(Volume), c_float, STRING,
                       c_int, POINTER(c_ushort)]
firstGuess.restype = None

prepVolume = librsl.prepVolume
prepVolume.argtypes = [POINTER(Volume), POINTER(Volume), c_float]
prepVolume.restype = None

unfoldVolume = librsl.unfoldVolume
unfoldVolume.argtypes = [POINTER(Volume), POINTER(Volume), POINTER(Volume),
                         c_float, c_ushort, POINTER(c_ushort)]
unfoldVolume.restype = None

c_float_p = POINTER(c_float)
firstGuessNoRead = librsl.firstGuessNoRead
firstGuessNoRead.argtypes = [POINTER(Volume), c_float, c_float_p, c_float_p,
                             c_float_p, c_int, c_int, POINTER(c_ushort)]
firstGuessNoRead.restype = None


# functions


# TODO move these to new module as they are not io?

def dealias_radar(radar, vad_time, last_radar=None, sondfile=None,
                  prep=1, DZindex=fieldTypes().DZ, filt=1):
    """
    Dealias a RSL Radar object.

    Parameters
    ----------
    radar : pointer to a Radar object (LP_Radar)
        RSL Radar object to dealias.
    vad_time : int
        Sounding time, format must be YYDDDHHMM.
    last_radar : pointer to Radar object (LP_Radar) or None
        RSL Radar object containg last measurements. Either last_radar or
        sondfile must be supplied.
    sondfile : str or None.
        Filename of sounding file.  Either sondfile or last_radar must be
        supplied.
    prep : int, optional
        Thresholding flag, 1 = yes, 0 = no.
    DZindex : int, optional
        Index of DZ measurements in RSL Radar object.
    filt : int. optional
        Flag controlling Bergen and Albers filter, 1 = yes, 0 = no.

    Returns
    -------
    unfoldedVolume : pointer to Radar object (LP_Radar)
        Dealised RSL Radar object.

    """
    # fixed parameters
    MISSINGVEL = 131072.0
    MAX_RADAR_VOLUMES = 41

    # create a new radar object to hold the dealiased results
    dealiased = RSL_new_radar(MAX_RADAR_VOLUMES)
    dealiased.contents.h = radar.contents.h
    dealiased.contents.h.nvolumes = MAX_RADAR_VOLUMES

    if radar.contents.volumes[DZindex] is not None:
        DZvolume = radar.contents.v[DZindex]
    else:
        print "no DZ field available"
        outDZ = 0
        prep = 0

    if radar.contents.volumes[fieldTypes().VR] is None:
        raise ValueError('No VR field avaialable to unfold')
    radialVelVolume = radar.contents.v[fieldTypes().VR]

    # store last Velocity Volumes from last_radar if present
    if last_radar is not None:
        lastVelVolume = last_radar.contents.volumes[fieldTypes().VE]
        numSweepsLast = lastVelVolume.h.nsweeps
        numSweepsCurrent = radialVelVolume.h.nsweeps
    else:
        lastVelVolume = None
        nolast = 1

    # read in sounding if present
    success = c_ushort(0)
    if os.path.isfile(sondfile):
        sondVolume = RSL_copy_volume(radialVelVolume)
        firstGuess(sondVolume, MISSINGVEL, sondfile, vad_time, success)
    else:
        sondVolume = None
        print "file does not exist"

    if success.value == 1 or lastVelVolume is not None:
        print "dealiasing"
        unfoldedVolume = RSL_copy_volume(radialVelVolume)
        if prep == 1:
            prepVolume(DZvolume, unfoldedVolume, MISSINGVEL)
        usuccess = c_ushort(0)
        unfoldVolume(unfoldedVolume, sondVolume, lastVelVolume, MISSINGVEL,
                     filt, usuccess)
    return unfoldedVolume


def dealias_radar_array(radar, vad_time, height, speed, direction, last_radar,
                        prep=1, DZindex=fieldTypes().DZ, filt=1):
    """
    Dealias a RSL Radar object using sond arrays.

    Parameters
    ----------
    radar : pointer to a Radar object (LP_Radar)
        RSL Radar object to dealias.
    vad_time : int
        Sounding time, format must be YYDDDHHMM.
    height : array_like
        Height of wind profile.
    speed : array_like
        Speed of wind profile.
    direction : array_like
        Direction of wind profile.
    last_radar : pointer to Radar object (LP_Radar) or None, optional.
        RSL Radar object containg last measurements.
    prep : int, optional
        Flag controlling thresholding, 1 = yes, 0 = no.
    DZindex : int, optional
        Index of DZ measurements in RSL Radar object.
    filt : int. optional
        Flag controlling Bergen and Albers filter, 1 = yes, 0 = no.


    Returns
    -------
    unfoldedVolume : pointer to Radar object (LP_Radar)
        Dealised RSL Radar object.
    sondVolume : XXX
        XXX

    """

    # prepare arrays
    hc = N.ascontiguousarray(height, dtype='float32')
    sc = N.ascontiguousarray(speed, dtype='float32')
    dc = N.ascontiguousarray(direction, dtype='float32')

    # fixed parameters
    MISSINGVEL = 131072.0
    MAX_RADAR_VOLUMES = 41

    # create a new radar object to hold dealised results
    dealiased = RSL_new_radar(MAX_RADAR_VOLUMES)
    dealiased.contents.h = radar.contents.h
    dealiased.contents.h.nvolumes = MAX_RADAR_VOLUMES

    if radar.contents.volumes[DZindex] is not None:
        DZvolume = radar.contents.v[DZindex]
    else:
        print "no DZ field available"
        outDZ = 0
        prep = 0

    if radar.contents.volumes[fieldTypes().VR] is None:
        raise ValueError('No VR field avaialable to unfold')
    radialVelVolume = radar.contents.v[fieldTypes().VR]

    # store last velocity volumes from last_data if present
    if last_radar is not None:
        lastVelVolume = last_radar.contents.volumes[fieldTypes().VE]
        numSweepsLast = lastVelVolume.h.nsweeps
        numSweepsCurrent = radialVelVolume.h.nsweeps
    else:
        lastVelVolume = None
        nolast = 1

    # read in sounding if present
    success = c_ushort(0)
    sondVolume = RSL_copy_volume(radialVelVolume)
    firstGuessNoRead(sondVolume, MISSINGVEL, hc.ctypes.data_as(c_float_p),
                     sc.ctypes.data_as(c_float_p),
                     dc.ctypes.data_as(c_float_p), c_int(len(height)),
                     vad_time, success)

    if success.value == 1 or lastVelVolume is not None:
        print "dealiasing"
        unfoldedVolume = RSL_copy_volume(radialVelVolume)
        if prep == 1:
            prepVolume(DZvolume, unfoldedVolume, MISSINGVEL)
        usuccess = c_ushort(0)
        unfoldVolume(unfoldedVolume, sondVolume, lastVelVolume, MISSINGVEL,
                     filt, usuccess)
    return unfoldedVolume, sondVolume


def getAllRays(radar, fieldType=None):
    """"
    Return a list of all rays from a single field in the radar structure.
    Defaults to the reflectivity field
    """

    if fieldType is None:
        fieldType = fieldTypes().DZ

    allrays = []

    swps = radar.contents.volumes[fieldType].sweeps
    for swp in swps:
        rays = swp.rays
        for ray in rays:
            allrays.append(ray)

    return allrays

vol_header_fields = dict(Volume_header._fields_)
f_prototype = vol_header_fields['f']
invf_prototype = vol_header_fields['invf']


# When creating a new radar structure, these need to be set for
# {vol,sweep,ray}.h.{f,invf}
conversion_functions = {
    'DZ': (f_prototype(('DZ_F', librsl)),
           invf_prototype(('DZ_INVF', librsl))),
    'VR': (f_prototype(('VR_F', librsl)),
           invf_prototype(('VR_INVF', librsl))),
    'SW': (SW_F, SW_INVF),
    'CZ': (CZ_F, CZ_INVF),
    'ZT': (ZT_F, ZT_INVF),
    'DR': (DR_F, DR_INVF),
    'LR': (LR_F, LR_INVF),
    'ZD': (ZD_F, ZD_INVF),
    'DM': (DM_F, DM_INVF),
    'RH': (RH_F, RH_INVF),
    'PH': (PH_F, PH_INVF),
    'XZ': (XZ_F, XZ_INVF),
    'CD': (CD_F, CD_INVF),
    'MZ': (MZ_F, MZ_INVF),
    'MD': (MD_F, MD_INVF),
    'ZE': (ZE_F, ZE_INVF),
    'VE': (VE_F, VE_INVF),
    'KD': (KD_F, KD_INVF),
    'TI': (TI_F, TI_INVF),
    'DX': (DX_F, DX_INVF),
    'CH': (CH_F, CH_INVF),
    'AH': (AH_F, AH_INVF),
    'CV': (CV_F, CV_INVF),
    'AV': (AV_F, AV_INVF),
    'SQ': (SQ_F, SQ_INVF),
}


#Some of Scott's helper routines
def does_volume_exist(radar, vol_string, fieldTypes):
    try:
        locn = N.where(N.array(fieldTypes().list) == vol_string)[0][0]
        vol_exists = (radar.contents.volumes[locn] is not None)
    except IndexError:
        print vol_string, " is not a RSL recognized field"
        vol_exists = False
    return vol_exists


def create_PPI_array(sweep):
    ppi = N.zeros([sweep.h.nrays, sweep.rays[0].h.nbins], dtype=float)
    for raynum in range(sweep.h.nrays):
        ppi[raynum, :] = sweep.rays[raynum].data
    return ppi


def return_header_metadata(header):
    retdict = {}
    for key in dir(header):
        if type(getattr(header, key)) in [float, int, str]:
            #print key, " ", getattr(header,key)
            retdict.update({key: getattr(header, key)})
    return retdict


ray_info = {
    'month': 'Time for this ray; month (1-12).',
    'day': 'Time for this ray; day (1-31).',
    'year': 'Time for this ray; year (eg. 1993).',
    'hour': 'Date for this ray; hour (0-23).',
    'minute': 'Date for this ray; minute (0-59).',
    'sec': 'Date for this ray; second + fraction of second.',
    'unam_rng': 'Unambiguous range. (KM).',
    'azimuth': """Azimuth angle. (degrees). Must be positive
    0=North, 90=east, -90/270=west.
    This angle is the mean azimuth for the whole ray.
    Eg. for NSIG the beginning and end azimuths are
    averaged.""",
    'ray_num': 'Ray no. within elevation scan',
    'elev': 'Elevation angle. (degrees).',
    'elev_num': 'Elevation no. within volume scan.',
    'range_bin1': 'Range to first gate.(meters).',
    'gate_size': 'Data gate size (meters).',
    'vel_res': 'Doppler velocity resolution.',
    'sweep_rate': 'Sweep rate. Full sweeps/min.',
    'prf': 'Pulse repetition frequency, in Hz.',
    'prf2': 'Second PRF, for Sigmet dual PRF',
    'azim_rate': 'Sweep rate in degrees/sec.',
    'fix_angle': 'Elevation angle for the sweep. (degrees).',
    'pitch': 'Pitch angle.',
    'roll': 'Roll  angle.',
    'heading': 'Heading.',
    'pitch_rate': '(angle/sec)',
    'roll_rate': '(angle/sec)',
    'heading_rate': '(angle/sec)',
    'lat': 'Latitude (degrees)',
    'lon': 'Longitude (degrees)',
    'alt': 'Altitude (m)',
    'rvc': 'Radial velocity correction (m/sec)',
    'vel_east': 'Platform velocity to the east (negative for west) (m/sec)',
    'vel_north': 'Platform velocity to the north (negative south) (m/sec)',
    'vel_up': 'Platform velocity toward up (negative down) (m/sec)',
    'pulse_count': 'Pulses used in a single dwell time.',
    'pulse_width': 'Pulse width (micro-sec).',
    'beam_width': 'Beamwidth in degrees.',
    'frequency': 'Carrier freq. GHz.',
    'wavelength': 'Wavelength. Meters.',
    'nyq_vel': 'Nyquist velocity. m/s.',
    'f': 'Data conversion function. f(x).',
    'invf': 'Data conversion function. invf(x).',
    'nbins': 'Number of array elements for \'Range\''
}


sweep_info = {
    'sweep_num': """Integer sweep number.  This may be redundant, since
    this will be the same as the Volume.sweep array index.""",
    'elev': 'Elevation angle (mean) for the sweep. Value is -999 for RHI.',
    'azimuth': 'Azimuth for the sweep (RHI). Value is -999 for PPI.',
    'beam_width': 'This is in the ray header too.',
    'vert_half_bw': 'Vertical beam width divided by 2',
    'horz_half_bw': 'Horizontal beam width divided by 2',
    'nrays': 'Number of rays in the sweep',
    'f': 'Data conversion function. f(x)',
    'invf': 'Data conversion function. invf(x).'
}


volume_info = {
    'type_str': "One of:'Reflectivity', 'Velocity' or 'Spectrum width'",
    'nsweeps': 'Number of sweeps in volume',
    'calibr_const': 'Calibration constant.  HDF specific.',
    'f': 'Data conversion function. f(x).',
    'invf': 'Data conversion function. invf(x).'
}


radar_info = {
    'month': 'Month of year', 'day': 'day of month',
    'year': 'year', 'hour': 'hour of day',
    'minute': 'Minute of hour',
    'sec': 'Second plus fractional part.',
    'radar_type': """Type of radar.  Use for QC-ing the data
    Supported types are:
    "wsr88d", "lassen", "uf",
    "nsig", "mcgill",
    "kwajalein", "rsl", "toga",
    "rapic", (rapic is Berrimah Austrailia)
    "radtec", (SPANDAR radar at Wallops Is, VA)
    "EDGE","dorade","south_africa". Set by appropriate ingest routine.""",
    'nvolumes': 'Number of volumes in the file',
    'number': 'arbitrary number of this radar site',
    'name': 'Nexrad site name',
    'radar_name': 'Radar name.',
    'project': 'Project identifier.',
    'city': 'nearest city to  radar site',
    'state': 'state of radar site',
    'country': 'Country of the radar site',
    'latd': 'degrees of latitude of site',
    'latm;': 'minutes of latitude of site',
    'lats': 'seconds of latitude of site',
    'lond': 'degrees of longitude of site',
    'lonm': 'minutes of longitude of site',
    'lons': 'seconds of longitude of site',
    'height': 'height of site in meters above sea level',
    'spulse': 'length of short pulse (ns)',
    'lpulse': 'length of long pulse (ns)',
    'scan_mode': ' 0 = PPI, 1 = RHI,',
    'vcp;': 'Volume Coverage Pattern (for WSR-88D only)'
}


def print_metadata(ray_header, info):
    meta = return_header_metadata(ray_header)
    for key in set(info.keys()) & set(meta.keys()):
        print key, ': ', meta[key], ' ', info[key]
