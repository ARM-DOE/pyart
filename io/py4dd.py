""" Python Wrapper for the University of Washington 4D dealias code

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-------------------------------------------------------------------------------


USE
---
The path to the RSL/4dd lib is hard coded
REQUIREMENTS
------------
RSL v.1.40 with 4dd code wrapped in 

HISTORY
-------
4/28/2010: Updated to work only with RSL v.1.40. 
That version added the vcp field to the Volume_header struct.

11/5/2010:Scott Collis Added extra headers and a few extra links to rsl functions
scollis.acrf@gmail.com

22/10/2011: Scott Collis added autosearch for the RSL library to make it platform independant 

"""

from ctypes import *
import numpy as N
import os.path
#from pylab 
import datetime
import sys

pyart_dir=os.environ.get('PYART_DIR',os.environ['HOME']+'/python')
lib_path=pyart_dir+'/pyart/io/lib/'
ext={'darwin':'dylib', 'linux2':'so', 'win32':'DLL'}[sys.platform]#note win32 not supported at this time...
files_in_dir=os.listdir(lib_path)
if 'librsl.'+ext in files_in_dir:
	print "rsl library found ", lib_path+ 'librsl.'+ext
	librslPath=lib_path+ 'librsl.'+ext



#librslPath = '/home/sc8/python/pyart/io/lib/librsl.so.1.0.41'
#librslPath = '/home/sc8/python/pyart/io/lib/librsl.so.1.0.43'
#librslPath = '/home/titan5/src/temp/lib/librsl.so.1.0.41'
#librslPath = '/bm/gdata/scollis/osra/lib/librsl.so'
# REQUIRES version 1.40 to match the updated volume header (added vcp)
# REQUIRES version 1.43 to fix sigmet issue

_libraries = {}
_libraries[librslPath] = CDLL(librslPath)

STRING = c_char_p
Range = c_ushort


class fieldTypes:
    def __init__(self):
        self.DZ =  0
        self.VR =  1
        self.SW =  2
        self.CZ =  3
        self.ZT =  4
        self.DR =  5
        self.LR =  6
        self.ZD =  7
        self.DM =  8
        self.RH =  9
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
        
        self.list = ['DZ',
                     'VR',
                     'SW',
                     'CZ',
                     'ZT',
                     'DR',
                     'LR',
                     'ZD',
                     'DM',
                     'RH',
                     'PH',
                     'XZ',
                     'CD',
                     'MZ',
                     'MD',
                     'ZE',
                     'VE',
                     'KD',
                     'TI',
                     'DX',
                     'CH',
                     'AH',
                     'CV',
                     'AV',
                     'SQ',]
                     
                     
DZ_F = _libraries[librslPath].DZ_F
DZ_F.restype = c_float
DZ_F.argtypes = [Range]
VR_F = _libraries[librslPath].VR_F
VR_F.restype = c_float
VR_F.argtypes = [Range]
SW_F = _libraries[librslPath].SW_F
SW_F.restype = c_float
SW_F.argtypes = [Range]
CZ_F = _libraries[librslPath].CZ_F
CZ_F.restype = c_float
CZ_F.argtypes = [Range]
ZT_F = _libraries[librslPath].ZT_F
ZT_F.restype = c_float
ZT_F.argtypes = [Range]
DR_F = _libraries[librslPath].DR_F
DR_F.restype = c_float
DR_F.argtypes = [Range]
LR_F = _libraries[librslPath].LR_F
LR_F.restype = c_float
LR_F.argtypes = [Range]
ZD_F = _libraries[librslPath].ZD_F
ZD_F.restype = c_float
ZD_F.argtypes = [Range]
DM_F = _libraries[librslPath].DM_F
DM_F.restype = c_float
DM_F.argtypes = [Range]
RH_F = _libraries[librslPath].RH_F
RH_F.restype = c_float
RH_F.argtypes = [Range]
PH_F = _libraries[librslPath].PH_F
PH_F.restype = c_float
PH_F.argtypes = [Range]
XZ_F = _libraries[librslPath].XZ_F
XZ_F.restype = c_float
XZ_F.argtypes = [Range]
CD_F = _libraries[librslPath].CD_F
CD_F.restype = c_float
CD_F.argtypes = [Range]
MZ_F = _libraries[librslPath].MZ_F
MZ_F.restype = c_float
MZ_F.argtypes = [Range]
MD_F = _libraries[librslPath].MD_F
MD_F.restype = c_float
MD_F.argtypes = [Range]
ZE_F = _libraries[librslPath].ZE_F
ZE_F.restype = c_float
ZE_F.argtypes = [Range]
VE_F = _libraries[librslPath].VE_F
VE_F.restype = c_float
VE_F.argtypes = [Range]
KD_F = _libraries[librslPath].KD_F
KD_F.restype = c_float
KD_F.argtypes = [Range]
TI_F = _libraries[librslPath].TI_F
TI_F.restype = c_float
TI_F.argtypes = [Range]
DX_F = _libraries[librslPath].DX_F
DX_F.restype = c_float
DX_F.argtypes = [Range]
CH_F = _libraries[librslPath].CH_F
CH_F.restype = c_float
CH_F.argtypes = [Range]
AH_F = _libraries[librslPath].AH_F
AH_F.restype = c_float
AH_F.argtypes = [Range]
CV_F = _libraries[librslPath].CV_F
CV_F.restype = c_float
CV_F.argtypes = [Range]
AV_F = _libraries[librslPath].AV_F
AV_F.restype = c_float
AV_F.argtypes = [Range]
SQ_F = _libraries[librslPath].SQ_F
SQ_F.restype = c_float
SQ_F.argtypes = [Range]
DZ_INVF = _libraries[librslPath].DZ_INVF
DZ_INVF.restype = Range
DZ_INVF.argtypes = [c_float]
VR_INVF = _libraries[librslPath].VR_INVF
VR_INVF.restype = Range
VR_INVF.argtypes = [c_float]
SW_INVF = _libraries[librslPath].SW_INVF
SW_INVF.restype = Range
SW_INVF.argtypes = [c_float]
CZ_INVF = _libraries[librslPath].CZ_INVF
CZ_INVF.restype = Range
CZ_INVF.argtypes = [c_float]
ZT_INVF = _libraries[librslPath].ZT_INVF
ZT_INVF.restype = Range
ZT_INVF.argtypes = [c_float]
DR_INVF = _libraries[librslPath].DR_INVF
DR_INVF.restype = Range
DR_INVF.argtypes = [c_float]
LR_INVF = _libraries[librslPath].LR_INVF
LR_INVF.restype = Range
LR_INVF.argtypes = [c_float]
ZD_INVF = _libraries[librslPath].ZD_INVF
ZD_INVF.restype = Range
ZD_INVF.argtypes = [c_float]
DM_INVF = _libraries[librslPath].DM_INVF
DM_INVF.restype = Range
DM_INVF.argtypes = [c_float]
RH_INVF = _libraries[librslPath].RH_INVF
RH_INVF.restype = Range
RH_INVF.argtypes = [c_float]
PH_INVF = _libraries[librslPath].PH_INVF
PH_INVF.restype = Range
PH_INVF.argtypes = [c_float]
XZ_INVF = _libraries[librslPath].XZ_INVF
XZ_INVF.restype = Range
XZ_INVF.argtypes = [c_float]
CD_INVF = _libraries[librslPath].CD_INVF
CD_INVF.restype = Range
CD_INVF.argtypes = [c_float]
MZ_INVF = _libraries[librslPath].MZ_INVF
MZ_INVF.restype = Range
MZ_INVF.argtypes = [c_float]
MD_INVF = _libraries[librslPath].MD_INVF
MD_INVF.restype = Range
MD_INVF.argtypes = [c_float]
ZE_INVF = _libraries[librslPath].ZE_INVF
ZE_INVF.restype = Range
ZE_INVF.argtypes = [c_float]
VE_INVF = _libraries[librslPath].VE_INVF
VE_INVF.restype = Range
VE_INVF.argtypes = [c_float]
KD_INVF = _libraries[librslPath].KD_INVF
KD_INVF.restype = Range
KD_INVF.argtypes = [c_float]
TI_INVF = _libraries[librslPath].TI_INVF
TI_INVF.restype = Range
TI_INVF.argtypes = [c_float]
DX_INVF = _libraries[librslPath].DX_INVF
DX_INVF.restype = Range
DX_INVF.argtypes = [c_float]
CH_INVF = _libraries[librslPath].CH_INVF
CH_INVF.restype = Range
CH_INVF.argtypes = [c_float]
AH_INVF = _libraries[librslPath].AH_INVF
AH_INVF.restype = Range
AH_INVF.argtypes = [c_float]
CV_INVF = _libraries[librslPath].CV_INVF
CV_INVF.restype = Range
CV_INVF.argtypes = [c_float]
AV_INVF = _libraries[librslPath].AV_INVF
AV_INVF.restype = Range
AV_INVF.argtypes = [c_float]
SQ_INVF = _libraries[librslPath].SQ_INVF
SQ_INVF.restype = Range
SQ_INVF.argtypes = [c_float]

class Ray_header(Structure):
    pass
Ray_header._fields_ = [
    ('month', c_int),
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
    ('nbins', c_int),
]

class Ray(Structure):
    #pass
    def __getattr__(self,attr):
        if attr=='data':
            nbins = self.h.nbins
            converter = self.h.f
            data = N.array([converter(self.range[i]) for i in range(nbins)])
            return data
        elif attr == 'dists':
            #return range gates in km
            dists = self.h.range_bin1 + self.h.gate_size * N.arange(self.h.nbins, dtype=float) 
            #dists = [RSL_get_range_of_range_index(ray, i) for i in range(ray.nbins)] --- in km
            return dists
        else:
            return self.h.__getattribute__(attr)
Ray._fields_ = [
    ('h', Ray_header),
    ('range', POINTER(Range)),
]

                     
class Sweep_header(Structure):
    pass
Sweep_header._fields_ = [
    ('sweep_num', c_int),
    ('elev', c_float),
    ('azimuth', c_float),
    ('beam_width', c_float),
    ('vert_half_bw', c_float),
    ('horz_half_bw', c_float),
    ('nrays', c_int),
    ('f', CFUNCTYPE(c_float, c_ushort)),
    ('invf', CFUNCTYPE(Range, c_float)),
]
class Sweep(Structure):
    #pass
    def __getattr__(self,attr):
        if attr=='rays':
            nrays = self.h.nrays
            rays = []
            for i in range(nrays):
                #There are null pointers where rays don't exist, which throws an error on access
                try: rays.append(self.ray[i].contents)
                except: rays.append(None)
            return rays
        else:
            return self.h.__getattribute__(attr)
Sweep._fields_ = [
    ('h', Sweep_header),
    ('ray', POINTER(POINTER(Ray))),
]

class Volume_header(Structure):
    pass
Volume_header._fields_ = [
    ('type_str', STRING),
    ('nsweeps', c_int),
    ('calibr_const', c_float),
    ('f', CFUNCTYPE(c_float, c_ushort)),
    ('invf', CFUNCTYPE(Range, c_float)),
]
class Volume(Structure):
    #pass
    def __getattr__(self,attr):
        if attr=='sweeps':
            nsweeps = self.h.nsweeps
            swps = []
            for i in range(nsweeps):
                #There are null pointers where sweeps don't exist, which throws an error on access
                try: swps.append(self.sweep[i].contents)
                except: swps.append(None)
            return swps
        else:
            return self.h.__getattribute__(attr)
Volume._fields_ = [
    ('h', Volume_header),
    ('sweep', POINTER(POINTER(Sweep))),
]

                    
class Radar_header(Structure):
    pass
Radar_header._fields_ = [
    ('month', c_int),
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
    ('scanmode',c_int),
    ('vcp', c_int),
]
class Radar(Structure):
    #pass
    def __getattr__(self,attr):
        if attr=='volumes':
            nvolumes = self.h.nvolumes
            vols = []
            for i in range(nvolumes):
                #There are null pointers where volumes don't exist, which throws an error on access
                try: vols.append(self.v[i].contents)
                except: vols.append(None)
                #vols = [self.v[i].contents for i in range(nvolumes)]
            return vols
        elif attr=='datetime':
            datetimeobj=datetime.datetime(self.h.year, self.h.month, self.h.day, self.h.hour, self.h.minute, int(self.h.sec))
            return datetimeobj
        else:
            return self.h.__getattribute__(attr)
Radar._fields_ = [
    ('h', Radar_header), ('v', POINTER(POINTER(Volume)))
]


RSL_anyformat_to_radar = _libraries[librslPath].RSL_anyformat_to_radar
RSL_anyformat_to_radar.restype = POINTER(Radar)
RSL_anyformat_to_radar.argtypes = [STRING]
RSL_radar_to_uf = _libraries[librslPath].RSL_radar_to_uf
RSL_radar_to_uf.argtypes = [POINTER(Radar),STRING]
RSL_radar_verbose_on = _libraries[librslPath].RSL_radar_verbose_on
RSL_radar_verbose_on.restype=None
RSL_free_radar = _libraries[librslPath].RSL_free_radar
RSL_free_radar.restype = None
RSL_free_radar.argtypes = [POINTER(Radar)]

RSL_new_radar = _libraries[librslPath].RSL_new_radar
RSL_new_radar.restype=POINTER(Radar)
RSL_new_radar.argtypes=[c_int]

RSL_new_volume = _libraries[librslPath].RSL_new_volume
RSL_new_volume.restype=POINTER(Volume)
RSL_new_volume.argtypes=[c_int]

RSL_new_sweep = _libraries[librslPath].RSL_new_sweep
RSL_new_sweep.restype=POINTER(Sweep)
RSL_new_sweep.argtypes=[c_int]

RSL_new_ray = _libraries[librslPath].RSL_new_ray
RSL_new_ray.restype=POINTER(Ray)
RSL_new_ray.argtypes=[c_int]

RSL_copy_volume = _libraries[librslPath].RSL_copy_volume
RSL_copy_volume.restype=POINTER(Volume)
RSL_copy_volume.argtypes=[POINTER(Volume)]

firstGuess=_libraries[librslPath].firstGuess
firstGuess.argtypes=[POINTER(Volume), c_float, STRING, c_int, POINTER(c_ushort)]
firstGuess.restype=None

prepVolume=_libraries[librslPath].prepVolume
prepVolume.argtypes=[POINTER(Volume), POINTER(Volume), c_float]
prepVolume.restype=None

unfoldVolume=_libraries[librslPath].unfoldVolume
unfoldVolume.argtypes=[POINTER(Volume), POINTER(Volume),POINTER(Volume), c_float, c_ushort, POINTER(c_ushort)]
unfoldVolume.restype=None

c_float_p = POINTER(c_float)
firstGuessNoRead=_libraries[librslPath].firstGuessNoRead
firstGuessNoRead.argtypes=[POINTER(Volume), c_float, c_float_p, c_float_p, c_float_p, c_int, c_int, POINTER(c_ushort)]
firstGuessNoRead.restype=None

#void firstGuessNoRead(Volume* soundVolume, float missingVal,
#          float *height_array, float *speed_array, float *direction_array, int nlevels, int VAD_time, unsigned short* sounding) {

#void unfoldVolume(Volume* rvVolume, Volume* soundVolume, Volume* lastVolume,
#     float missingVal, unsigned short rm, unsigned short* success);

#
#void firstGuess(Volume* soundVolume, float missingVal,
#          char* sounding_name, int VAD_time, unsigned short* sounding)
# __all__ = ['RSL_anyformat_to_radar', Radar, Volume, Sweep, Range,]

def dealias_radar(radar,last_radar, sondfile, vad_time, **kwargs):
    prep=kwargs.get('prep',1)
    LOWDBZ=kwargs.get('LOWDBZ',-1.0)
    HIGHDBZ=kwargs.get('HIGHDZ', 80.0)
    MISSINGVEL=131072.0
    MAX_RADAR_VOLUMES=41
    dealiased=RSL_new_radar(MAX_RADAR_VOLUMES)
    dealiased.contents.h=radar.contents.h
    dealiased.contents.h.nvolumes=MAX_RADAR_VOLUMES
    if radar.contents.volumes[kwargs.get('DZindex', fieldTypes().DZ)] != None:
        DZvolume= radar.contents.v[kwargs.get('DZindex', fieldTypes().DZ)]
    else:
        print "no DZ field available"
        outDZ=0
        prep=0
    if radar.contents.volumes[fieldTypes().VR] != None:
        radialVelVolume = radar.contents.v[fieldTypes().VR]
    else:
        print "No VR field available to unfold; aborting"
        return None
    if last_radar != None:
        lastVelVolume=last_radar.contents.volumes[fieldTypes().VE]
        numSweepsLast=lastVelVolume.h.nsweeps
        numSweepsCurrent = radialVelVolume.h.nsweeps
    else:
        lastVelVolume=None
        nolast=1
    print "Trying to read sounding"
    success=c_ushort(0)
    if os.path.isfile(sondfile):
        sondVolume=RSL_copy_volume(radialVelVolume)
        firstGuess(sondVolume,MISSINGVEL, sondfile, vad_time, success)
    else:
        sondVolume = None
        print "file does not exist"
    if success.value==1 or lastVelVolume!=None:
        print "dealiasing"
        print "High DBZ: ", HIGHDBZ, "Low DBZ:", LOWDBZ
        unfoldedVolume = RSL_copy_volume(radialVelVolume)
        if prep==1: prepVolume(DZvolume, unfoldedVolume, MISSINGVEL)
        usuccess=c_ushort(0)
        unfoldVolume(unfoldedVolume, sondVolume,lastVelVolume,MISSINGVEL, kwargs.get('filt',1), usuccess );
    return unfoldedVolume

def dealias_radar_array(radar,last_radar, height, speed, direction, vad_time, **kwargs):
    hc=N.ascontiguousarray(height, dtype='float32')
    sc=N.ascontiguousarray(speed, dtype='float32')
    dc=N.ascontiguousarray(direction, dtype='float32')
    prep=kwargs.get('prep',1)
    LOWDBZ=kwargs.get('LOWDBZ',-1.0)
    HIGHDBZ=kwargs.get('HIGHDZ', 80.0)
    MISSINGVEL=131072.0
    MAX_RADAR_VOLUMES=41
    dealiased=RSL_new_radar(MAX_RADAR_VOLUMES)
    dealiased.contents.h=radar.contents.h
    dealiased.contents.h.nvolumes=MAX_RADAR_VOLUMES
    if radar.contents.volumes[kwargs.get('DZindex', fieldTypes().DZ)] != None:
        DZvolume= radar.contents.v[kwargs.get('DZindex', fieldTypes().DZ)]
    else:
        print "no DZ field available"
        outDZ=0
        prep=0
    if radar.contents.volumes[fieldTypes().VR] != None:
        radialVelVolume = radar.contents.v[fieldTypes().VR]
    else:
        print "No VR field available to unfold; aborting"
        return None
    if last_radar != None:
        lastVelVolume=last_radar.contents.volumes[fieldTypes().VE]
        numSweepsLast=lastVelVolume.h.nsweeps
        numSweepsCurrent = radialVelVolume.h.nsweeps
    else:
        lastVelVolume=None
        nolast=1
    print "Trying to read sounding"
    success=c_ushort(0)
    #if os.path.isfile(sondfile):
    sondVolume=RSL_copy_volume(radialVelVolume)
    firstGuessNoRead(sondVolume,MISSINGVEL,  hc.ctypes.data_as(c_float_p), sc.ctypes.data_as(c_float_p),dc.ctypes.data_as(c_float_p), c_int(len(height)), vad_time, success)
    #void firstGuessNoRead(Volume* soundVolume, float missingVal,
    #          float *height_array, float *speed_array, float *direction_array, int nlevels, int VAD_time, unsigned short* sounding) {
    #else:
    #    sondVolume = None
    #    print "file does not exist"
    if success.value==1 or lastVelVolume!=None:
        print "dealiasing"
        print "High DBZ: ", HIGHDBZ, "Low DBZ:", LOWDBZ
        unfoldedVolume = RSL_copy_volume(radialVelVolume)
        if prep==1: prepVolume(DZvolume, unfoldedVolume, MISSINGVEL)
        usuccess=c_ushort(0)
        unfoldVolume(unfoldedVolume, sondVolume,lastVelVolume,MISSINGVEL, kwargs.get('filt',1), usuccess );
    return unfoldedVolume, sondVolume



#      if (success==0) {
#      soundVolume=NULL;
#    }
#
#    if (success==1 || lastVelVolume!=NULL) { /* Proceed with unfolding */
#
#      printf("Dealiasing: \n");
#      /* Copy radial velocity data to unfoldedVolume and remove any
#      **   bins where reflectivity is missing (if NODBZRMRV==1) or
#      **   ouside the interval LOWDBZ<= dz <=HIGHDBZ. */      
#     printf("LOWDBZ= %f\n", LOWDBZ);
#     printf("HIGHDBZ= %f\n", HIGHDBZ);
#     unfoldedVolume = RSL_copy_volume(radialVelVolume);
#      if (prep) prepVolume(DZvolume, unfoldedVolume, MISSINGVEL);
#      
#      /* Finally, we call unfoldVolume, which unfolds the velocity
#      ** field in the array unfoldedVolume. The subroutine returns a value
#      ** of 1 if unfolding is performed. */
#      unfoldVolume(unfoldedVolume,
#           soundVolume,lastVelVolume,MISSINGVEL, filt, &usuccess);
#    }
#
#    if (unfoldedVolume==NULL || usuccess==0) {
#      printf("Velocity Volume was NOT unfolded!\n");
#      unfoldedVolume=RSL_copy_volume(radialVelVolume);
#    }
#    else {
#      printf("Velocity Volume was unfolded.\n");
#    }      
#   
#
#    
    
    


def getAllRays(radar, fieldType=None):
    "Return a list of all rays from a single field in the radar structure. Defaults to the reflectivity field"

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
# DZ_F_CFUNC = f_prototype(('DZ_F', _libraries[librslPath]))



# When creating a new radar structure, these need to be set for {vol,sweep,ray}.h.{f,invf}
conversion_functions = {
    'DZ' : (f_prototype(('DZ_F', _libraries[librslPath])), invf_prototype(('DZ_INVF', _libraries[librslPath]))),
    'VR' : (f_prototype(('VR_F', _libraries[librslPath])), invf_prototype(('VR_INVF', _libraries[librslPath]))),
    'SW' : (SW_F, SW_INVF),
    'CZ' : (CZ_F, CZ_INVF),
    'ZT' : (ZT_F, ZT_INVF),
    'DR' : (DR_F, DR_INVF),
    'LR' : (LR_F, LR_INVF),
    'ZD' : (ZD_F, ZD_INVF),
    'DM' : (DM_F, DM_INVF),
    'RH' : (RH_F, RH_INVF),
    'PH' : (PH_F, PH_INVF),
    'XZ' : (XZ_F, XZ_INVF),
    'CD' : (CD_F, CD_INVF),
    'MZ' : (MZ_F, MZ_INVF),
    'MD' : (MD_F, MD_INVF),
    'ZE' : (ZE_F, ZE_INVF),
    'VE' : (VE_F, VE_INVF),
    'KD' : (KD_F, KD_INVF),
    'TI' : (TI_F, TI_INVF),
    'DX' : (DX_F, DX_INVF),
    'CH' : (CH_F, CH_INVF),
    'AH' : (AH_F, AH_INVF),
    'CV' : (CV_F, CV_INVF),
    'AV' : (AV_F, AV_INVF),
    'SQ' : (SQ_F, SQ_INVF),
    }

#Some of Scott's helper routines
def does_volume_exist(radar, vol_string, fieldTypes):
    try:
        locn=N.where(N.array(fieldTypes().list) == vol_string)[0][0]
        vol_exists=radar.contents.volumes[locn]!=None
    except IndexError:
        print vol_string, " is not a RSL recognized field"
        vol_exists=False
    return vol_exists

def create_PPI_array(sweep):
    ppi=N.zeros([sweep.h.nrays, sweep.rays[0].h.nbins], dtype=float)
    for raynum in range(sweep.h.nrays):
        ppi[raynum, :]=sweep.rays[raynum].data
    return ppi

def return_header_metadata(header):
    retdict={}
    for key in dir(header):
        if type(getattr(header,key)) in [float, int, str]:
            #print key, " ", getattr(header,key)
            retdict.update({key:getattr(header,key)})
    return retdict

def ray_info():
    info={'month'   :'Time for this ray; month (1-12).',
          'day'     :'Time for this ray; day (1-31).',
          'year'    :'Time for this ray; year (eg. 1993).',
          'hour'    :'Date for this ray; hour (0-23).',
          'minute'  :'Date for this ray; minute (0-59).',
          'sec'     :'Date for this ray; second + fraction of second.',
          'unam_rng':'Unambiguous range. (KM).',
          'azimuth' :"""Azimuth angle. (degrees). Must be positive
          0=North, 90=east, -90/270=west.
          This angle is the mean azimuth for the whole ray.
          Eg. for NSIG the beginning and end azimuths are
          averaged.""",
          'ray_num' :'Ray no. within elevation scan',
          'elev':'Elevation angle. (degrees).',
          'elev_num':'Elevation no. within volume scan.',
          'range_bin1':'Range to first gate.(meters).',
          'gate_size':'Data gate size (meters).',
          'vel_res':'Doppler velocity resolution.',
          'sweep_rate':'Sweep rate. Full sweeps/min.',
          'prf':'Pulse repetition frequency, in Hz.',
          'prf2':'Second PRF, for Sigmet dual PRF',
          'azim_rate':'Sweep rate in degrees/sec.',
          'fix_angle':'Elevation angle for the sweep. (degrees).',
          'pitch':'Pitch angle.',
          'roll':'Roll  angle.',
          'heading':'Heading.',
          'pitch_rate':'(angle/sec)',
          'roll_rate':'(angle/sec)',
          'heading_rate':'(angle/sec)',
          'lat':'Latitude (degrees)',
          'lon':'Longitude (degrees)',
          'alt':'Altitude (m)',
          'rvc':'Radial velocity correction (m/sec)',
          'vel_east':'Platform velocity to the east (negative for west) (m/sec)',
          'vel_north':'Platform velocity to the north (negative south) (m/sec)',
          'vel_up':'Platform velocity toward up (negative down) (m/sec)',
          'pulse_count':'Pulses used in a single dwell time.',
          'pulse_width':'Pulse width (micro-sec).',
          'beam_width':'Beamwidth in degrees.',
          'frequency':'Carrier freq. GHz.',
          'wavelength':'Wavelength. Meters.',
          'nyq_vel':'Nyquist velocity. m/s.',
          'f':'Data conversion function. f(x).',
          'invf':'Data conversion function. invf(x).',
          'nbins':'Number of array elements for \'Range\''}
    return info 

def sweep_info():
    retdict={'sweep_num':"""Integer sweep number.  This may be redundant, since
             this will be the same as the Volume.sweep array index.""",
            'elev':'Elevation angle (mean) for the sweep. Value is -999 for RHI.',
            'azimuth':'Azimuth for the sweep (RHI). Value is -999 for PPI.',
            'beam_width':'This is in the ray header too.',
            'vert_half_bw':'Vertical beam width divided by 2',
            'horz_half_bw':'Horizontal beam width divided by 2',
            'nrays':'Number of rays in the sweep',
            'f':'Data conversion function. f(x)',
            'invf':'Data conversion function. invf(x).'}
    return retdict

def volume_info():
    retdict={'type_str': "One of:'Reflectivity', 'Velocity' or 'Spectrum width'",
             'nsweeps':'Number of sweeps in volume',
             'calibr_const':'Calibration constant.  HDF specific.',
             'f':'Data conversion function. f(x).',
             'invf':' Data conversion function. invf(x).'}
    return retdict

def radar_info():
    retdict={'month':'Month of year', 'day':'day of month', 'year':'year','hour':'hour of day', 'minute':'Minute of hour',
             'sec':'Second plus fractional part.',
             'radar_type':"""Type of radar.  Use for QC-ing the data
             Supported types are:
             "wsr88d", "lassen", "uf",
             "nsig", "mcgill",
              "kwajalein", "rsl", "toga",
              "rapic", (rapic is Berrimah Austrailia)
              "radtec", (SPANDAR radar at Wallops Is, VA)
              "EDGE","dorade","south_africa". Set by appropriate ingest routine.""",
              'nvolumes':'Number of volumes in the file',
              'number':'arbitrary number of this radar site',
              'name':'Nexrad site name',
              'radar_name':'Radar name.',
              'project':'Project identifier.',
              'city':'nearest city to  radar site',
              'state':'state of radar site',
              'country':'Country of the radar site',
              'latd':'degrees of latitude of site',
              'latm;':'minutes of latitude of site',
              'lats':'seconds of latitude of site',
              'lond':'degrees of longitude of site',
              'lonm':'minutes of longitude of site',
              'lons':'seconds of longitude of site',
              'height':'height of site in meters above sea level',
              'spulse':'length of short pulse (ns)',
              'lpulse':'length of long pulse (ns)',
              'scan_mode':' 0 = PPI, 1 = RHI,',
              'vcp;':'Volume Coverage Pattern (for WSR-88D only)'}
    return retdict

def print_metadata(ray_header, info):
    meta=return_header_metadata(ray_header)
    for key in set(info().keys())&set(meta.keys()):
        print key, ': ', meta[key], ' ', info()[key]



