int dealias_fourdd(
    Volume* rvVolume, Volume* soundVolume, Volume* lastVolume,
    float missingVal, float compthresh, float compthresh2, float thresh,
    float ckval, float stdthresh, float epsilon,
    int maxcount, int pass2, int rm, int proximity, int mingood, int filt,
    int ba_mincount, int ba_edgecount);
