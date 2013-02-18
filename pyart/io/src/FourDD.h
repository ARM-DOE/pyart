/*
**
**  UW Radial Velocity Unfolding Algorithm
**  Four-Dimensional Dealiasing (4DD)
**
**  DESCRIPTION:
**     This algorithm unfolds a volume of single Doppler radial velocity data.
**  The algorithm uses a previously unfolded volume (or VAD if previous volume
**  is unavailable) and a higher elevation sweep to unfold some of the bins
**  in each sweep. Then, it completes the unfolding, assuming spatial
**  continuity around each bin.
**
**  DEVELOPER:
**     Curtis N. James            25 Jan 99
**
*/

#ifndef FDD_H
#define FDD_H
#define MAXRAYS   375      /* added by SRB 980310 */
#define MAXBINS 2048
#define HIGHDBZTHRESHOLD (float) 80.0
#define VERBOSE 0  /* Verbose=1 for detailed printout during execution */
#define PROXIMITY 5 /* For unfolding using windowing.*/
#define COMPTHRESH 0.25 /* The threshold for performing initial dealiasing 
			** using a previously unfolded volume. */
#define COMPTHRESH2 0.49 /* The threshold for performing initial dealiasing 
			** using sounding (or VAD). */
#define THRESH 0.4 /* The unfolding threshold for unfolding using horizontal
		   ** continuity. */
#define MINGOOD 5 /* Number of good values required within unfolding window
		  **  to unfold the current bin. */
#define STDTHRESH 0.8 /* Fraction of the Nyquist velocity to use as a standard
		      **  deviation threshold when windowing. */
#define LOWDBZ 0.0 /* All VR bins with DZ values less than LOWDBZ will be
		     **  deleted. */
#define HIGHDBZ 80.0 /* All bins with DZ values greater than HIGHDBZ will
		     ** be deleted. */
#define NODBZRMRV 0 /* If NODBZRMRV==1, all VR bins with DZ values missing
		    ** will be deleted. */
#define RM 0 /* If soundvolume is not available, remove cells left over after
	     **  first pass. */
#define PASS2 1 /* Flag specifying the use of a second pass using only the
	        **   sounding (or VAD).*/
#define DELNUM 0 /* The first DELNUM velocity bins will be deleted along each
		 **  ray (should be between 0 and 5). */
#define CKVAL 1.0 /* If absolute value of the radial velocity gate is less 
		  ** than this value, it will not be used as a PRELIM gate. */
#define SIGN 1 /* Sign convention: if SIGN=-1, negative radial velocity is
	       ** towards the radar, if SIGN=1 positive value towards radar. */
#define MAXCOUNT 10 /* Maximum number of folds. */
#define A 6372.0 /* Radius of earth in kilometers. */
#define PI 3.1415927 
#define MAXSHEAR 0.05 /* Maximum vertical shear allowed in input sounding */
#include <rsl.h> /* Sweep */ 

void firstGuess(Volume* soundVolume, float missingVal, char* sounding_name,
     int VAD_time, unsigned short* sounding);

void unfoldVolume(Volume* rvVolume, Volume* soundVolume, Volume* lastVolume,
     float missingVal, unsigned short rm, unsigned short* success);

float window(Volume* rvVolume, int sweepIndex, int startray, int endray,
     int firstbin, int lastbin, float* std, float missingVal, unsigned 
     short* success);

void prepVolume (Volume* DBZVolume, Volume* rvVolume, float missingVal);

int findRay (Volume* rvVolume1, Volume* rvVolume2, int sweepIndex1, int 
	sweepIndex2, int currIndex1, float missingVal);

float previousVal (Volume* rvVolume, Volume* lastVolume, int sweepIndex, int
	currIndex, int rangeIndex, float missingVal);

#endif /* DEALIAS_H */










