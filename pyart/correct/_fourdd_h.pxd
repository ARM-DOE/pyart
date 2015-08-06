"""
Cython wrapper around University of Washington 4DD code.
"""

cimport pyart.io._rsl_h as _rsl_h

cdef extern from "dealias_fourdd.h":

    int dealias_fourdd(_rsl_h.Volume* rvVolume, _rsl_h.Volume* soundVolume,
                       _rsl_h.Volume* lastVolume, float missingVal,
                       float compthresh, float compthresh2, float thresh,
                       float ckval, float stdthresh, float epsilon, 
                       int maxcount, int pass2, int rm, int proximity,
                       int mingood, int filt, int ba_mincount, 
                       int ba_edgecount)


cdef extern from "sounding_to_volume.h":
    
    int sounding_to_volume(_rsl_h.Volume* soundVolume, float missingVal,
                           float *height_array, float *speed_array,
                           float *direction_array, int nlevels, 
                           float maxshear, int sign)
