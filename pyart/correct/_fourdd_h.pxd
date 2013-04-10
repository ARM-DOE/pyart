""" Cython wrapper around NASA TRMM RSL Library """

# data structures
cimport pyart.io._rsl_h as _rsl_h


cdef extern from "FourDD.h":
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
