"""
Cython wrapper around University of Washington 4DD code.
"""

cimport pyart.io._rsl_h as _rsl_h

cdef extern from "FourDD.h":

    cdef int MAXRAYS
    cdef int MAXBINS
    cdef float HIGHDZBTHRESHOLD
    cdef int VERBOSE
    cdef int PROXIMITY
    cdef float COMPTHRESH
    cdef float COMPTHRESH2
    cdef float THRESH
    cdef int MINGOOD
    cdef float STDTHRESH
    cdef float LOWDBZ
    cdef float HIGHDBZ
    cdef int NODBZRMRV
    cdef int RM
    cdef int PASS2
    cdef int DELNUM
    cdef float CKVAL
    cdef int SIGN
    cdef int MAXCOUNT
    cdef float MAXSHEAR

    void firstGuess(_rsl_h.Volume* soundVolume, float missingVal,
                    char* sounding_name, int VAD_time,
                    unsigned short* sounding)

    void firstGuessNoRead(_rsl_h.Volume* soundVolume, float missingVal,
                          float *height_array, float *speed_array,
                          float *direction_array, int nlevels, int VAD_time,
                          unsigned short* sounding)

    void unfoldVolume(_rsl_h.Volume* rvVolume, _rsl_h.Volume* soundVolume,
                      _rsl_h.Volume* lastVolume, float missingVal,
                      unsigned short rm, unsigned short* success)

    void prepVolume(_rsl_h.Volume* DBZVolume, _rsl_h.Volume* rvVolume,
                    float missingVal)

    float window(_rsl_h.Volume* rvVolume, int sweepIndex, int startray,
                 int endray, int firstbin, int lastbin, float* std,
                 float missingVal, unsigned short* success)

    int findRay(_rsl_h.Volume* rvVolume1, _rsl_h.Volume* rvVolume2,
                int sweepIndex1, int sweepIndex2, int currIndex1,
                float missingVal)

    float previousVal(_rsl_h.Volume* rvVolume, _rsl_h.Volume* lastVolume,
                      int sweepIndex, int currIndex, int rangeIndex,
                      float missingVal)
