/*
 * Unfold the doppler velocity data for a single radar volume using 
 * sounding data and/or a priviously unfolded volume.
 *  
 * Adapted from routines f* 
 * UW Radial Velocity Dealiasing Algorithm
 * Four-Dimensional Dealiasing (4DD)
 * 
 * Developer by Curtis N. James     Jan-Feb 1999
 *
 * See UWASHINGTON_4DD_README file for license
 *
 * Adapted for use in Py-ART by Jonathan J. Helmus May 2014
 */


/*
**
**
**  DESCRIPTION:
**     This algorithm unfolds a volume of single Doppler radial velocity data.
**  The algorithm uses a previously unfolded volume (or VAD if previous volume
**  is unavailable) and a higher elevation sweep to unfold some of the bins
**  in each sweep. Then, it completes the unfolding, assuming spatial
**  continuity around each bin.
**
**
*/

#include "helpers.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <rsl.h> /* Sweep */ 

/******************************
 * Private data and functions *
 ******************************/

typedef struct DealiasParams {
    float missingVal;   /* value used to mark missing values in sweeps */ 
    float epsilon;      /* minimum distance to use when checking for missingVal */
    float compthresh;   /* The threshold for performing initial dealiasing 
                         * using a previously unfolded volume. */ 
    float ckval;        /* If absolute value of the radial velocity gate is less 
                         * than this value, it will not be used as a PRELIM gate. */
    float thresh;       /* The unfolding threshold for unfolding using horizontal
                         * continuity. between 0.0 and 1.0*/
    int maxcount;       /* Maximum number of folds. */
    int pass2;          /* Flag specifying the use of a second pass using only the
                         *  sounding (or VAD).*/
    int rm;             /* If soundvolume is not available, remove cells left
                         * over after first pass. */
    int proximity;      /* For unfolding using windowing.*/
    float compthresh2;  /* The threshold for performing initial dealiasing 
                         * using sounding (or VAD). */
    int mingood;        /* Number of good values required within unfolding window
                         *  to unfold the current bin. */
    float stdthresh;    /* Fraction of the Nyquist velocity to use as a standard
                         *  deviation threshold when windowing. */
    int ba_mincount;    /* Minimum number of non-missing neighbors to require
                         * during Berger and Albers 3x3 filter, a value of
                         * 5 will include points where a majority of neighbors
                         * are non-missing (MUST be between 1 and 8). */
    int ba_edgecount;   /* Same as ba_mincount except used at ray edges
                         * which have only 5 neightbors, a value of 3 will 
                         * keep those gates where a majority of neightbors
                         * are non-missing (MUST be between 1 and 5). */
} DealiasParams;

/* Structure for storing data on sweeps being examined for unfolding */
typedef struct SweepCollection {
    Sweep *vals;
    Sweep *rv;
    Sweep *last;
    Sweep *sound;
    Sweep *above;
    float NyqVelocity;
    float NyqInterval;
    int nrays;
    int nbins;
} SweepCollection;

/*
 * Marks structure and getter and setter.
 * All interactions with marks should be done through these functions.
 * 
 */

typedef struct Marks {
    short *mark_data;       /* Array holding gate marks */
    int maxrays;            /* Maximum number of rays in any sweep */
    int maxbins;            /* Maximum number of bins in any sweep */
} Marks;

/* set mark */
static void set_mark(Marks *marks, int r_idx, int b_idx, short val)
{
    marks->mark_data[b_idx * marks->maxrays + r_idx] = val;
    return;
}

/* get mark */
static short get_mark(Marks *marks, int r_idx, int b_idx)
{
    return marks->mark_data[b_idx * marks->maxrays + r_idx];
}

/* Set a mark and the corresponding ray value in a sweep */
static void set_mark_ray(Marks *marks, Sweep *sweep, 
    int ray_idx, int bin_idx, short mark_val, float ray_val)
{
    ray_set(sweep->ray[ray_idx], bin_idx, ray_val);
    set_mark(marks, ray_idx, bin_idx, mark_val);
    return;
}


/* Check a value against missing value */
static int is_missing(float val, DealiasParams *dp)
{
    //if (val == dp->missingVal)
    if (fabs(val - dp->missingVal) < dp->epsilon)
        return 1;
    else
        return 0;
}

static int not_missing(float val, DealiasParams *dp)
{
    //if (val == dp->missingVal)
    if (fabs(val - dp->missingVal) < dp->epsilon)
        return 0;
    else
        return 1;
}

/* 
 * Functions used by multiple other functions
 */

/* Minize the difference between val and cval by adding or substracting 
 * the Nyquist Interval, update val. */
static float min_diff(float *val, float cval, float NyqVelocity, int maxcount)
{
    int numtimes; 
    numtimes=0;
    while (fabs(cval-*val) > 0.99999 * NyqVelocity && numtimes <= maxcount) {
        numtimes++;
        if (*val > cval)
            *val-=2 * NyqVelocity;
        else
            *val+=2 * NyqVelocity;
    }
    return fabs(cval-*val);
}

/*
 * Berger and Albers filter functions
 */

/* Count the number of neighbors which are not missing value */
static int count_nonmissing_neighbors(
    Sweep *sweep, int ray_index, int i, DealiasParams *dp)
{
    int countindex = 0;
    int left, right, next, prev;

    if (ray_index == 0) 
        left = sweep->h.nrays-1;
    else 
        left = ray_index - 1;
    if (ray_index == sweep->h.nrays-1) 
        right = 0;
    else 
        right = ray_index+1;
    next = i+1;
    prev = i-1;
    /* Look at all bins adjacent to current bin in question: */
    if (i != 0) {
        if (not_missing(ray_val(sweep->ray[left], prev), dp))
            countindex++;
        if (not_missing(ray_val(sweep->ray[ray_index], prev), dp))
            countindex++;
        if (not_missing(ray_val(sweep->ray[right], prev), dp))
            countindex++;
    }
    if (not_missing(ray_val(sweep->ray[left], i), dp))
        countindex++;
    if (not_missing(ray_val(sweep->ray[right], i), dp))
        countindex++;
    if (i < sweep->ray[0]->h.nbins-1) {  
        if (not_missing(ray_val(sweep->ray[left], next), dp))
            countindex++;
        if (not_missing(ray_val(sweep->ray[ray_index], next), dp))
            countindex++;
        if (not_missing(ray_val(sweep->ray[right], next), dp))
            countindex++;
    }
    return countindex;
}

/* Perform a 3x3 missing value Berger and Albers filter */
static void berger_albers_filter(
    Sweep *sweep, DealiasParams *dp, Marks *marks)
{
    int ray_index, i, count;
    float val;
    for (ray_index=0; ray_index < sweep->h.nrays; ray_index++) {
        for (i=0; i < sweep->ray[0]->h.nbins; i++) {
            val = ray_val(sweep->ray[ray_index], i); 
            if (is_missing(val, dp)) {
                set_mark(marks, ray_index, i, -1);
            } else {
                count = count_nonmissing_neighbors(
                    sweep, ray_index, i, dp);
                if (count >= dp->ba_mincount)
                    set_mark(marks, ray_index, i, 0);
                else if (i == (sweep->ray[0]->h.nbins - 1) && 
                         count >= dp->ba_edgecount)
                    set_mark(marks, ray_index, i, 0);
                else
                    set_mark(marks, ray_index, i, -1);
            }
        }
    }
}

/*
 * Zero-ing marks array function
 */

/* Mark with zeros all location where the value in sweep is not missing, mark
 * with zero if missing */
static void init_marks(Sweep *sweep, DealiasParams *dp, Marks *marks)
{
    int ray_index, i;
    float val;
    for (ray_index=0; ray_index < sweep->h.nrays; ray_index++) {
        for (i=0; i < sweep->ray[0]->h.nbins; i++) {
            val = ray_val(sweep->ray[ray_index], i); 
            if (is_missing(val, dp))
                set_mark(marks, ray_index, i, -1);
            else
                set_mark(marks, ray_index, i, 0);
        }
    }
}

/*
 * Sweep initialization
 */

/* Initalize (fill) a sweep with a given value */
static void init_sweep(Sweep *sweep, float val)
{
    int ray_index, i;
    for (ray_index=0; ray_index < sweep->h.nrays; ray_index++) {
        for (i=0; i<sweep->ray[0]->h.nbins; i++) {
            ray_set(sweep->ray[ray_index], i, val);
        }
    }
}

/*
 * Continuity Dealising (initial dealias) functions 
 */

/* Finds the rayindex of the nearest ray in sweep2 to a ray in sweep1. */
static int findRay(Sweep *sweep1, Sweep *sweep2, int ray_index) 
{
    int numRays, rayIndex1;
    int minraynum; 
    int i;
    float az0, az1, diffaz, currmaxdiffaz;
    currmaxdiffaz=99999999;
    
    numRays = sweep2->h.nrays;
    az0 = sweep1->ray[ray_index]->h.azimuth;

    for(i=0; i<numRays; i++){
        az1 = sweep2->ray[i]->h.azimuth;
        diffaz = fabs(az1-az0);
        if(diffaz>180.0){
            diffaz = 360 - diffaz;
        }
        if(diffaz<currmaxdiffaz){
            currmaxdiffaz = diffaz;
            minraynum = i;
        }


    }    
    return minraynum;
}



/* Unfold bins where vertical and temporal continuity
 * produces the same value (i.e., bins where the radial velocity
 * agrees with both the previous sweep and previous volume within
 * a COMPTHRESH of the Nyquist velocity). This produces a number
 * of good data points from which to begin unfolding assuming spatial
 * continuity. */
static void continuity_dealias(
    SweepCollection *sweepc, DealiasParams *dp, Marks *marks)
{
    int ray_index;
    int i, prevIndex = 0, abIndex = 0; 
    float val;
    float prevval, soundval, abval, diff, thresh;
    float spacing;
    short valid_above_ray;
    valid_above_ray = 1;

    thresh = dp->compthresh * sweepc->NyqVelocity;

    for (ray_index=0; ray_index < sweepc->nrays; ray_index++) { 
        
        if (sweepc->last != NULL){
           prevIndex = findRay(sweepc->rv, sweepc->last, ray_index);
        }

        if (sweepc->above != NULL){
            valid_above_ray = 0;
            //Checking to make sure that the above ray is actually *above*.
            //This takes an average of the azimuthal spacing over the first four rays
            spacing = fabs(sweepc->above->ray[0]->h.azimuth - sweepc->above->ray[4]->h.azimuth);
            if (spacing>180.0)
                spacing = 360.0 - spacing;
            spacing = spacing / 4.0;

            abIndex = findRay(sweepc->rv, sweepc->above, ray_index);
            float aznew = sweepc->above->ray[abIndex]->h.azimuth;
            float azold = sweepc->rv->ray[ray_index]->h.azimuth;
            float diffaz_newold = fabs(aznew-azold);
            if(diffaz_newold > 180.0)
                diffaz_newold = 360.0-diffaz_newold;
            if(fabs(aznew-azold)>spacing/2){
                valid_above_ray = 1;
            }
        }
        for (i=0; i < sweepc->nbins; i++) { 


            if (get_mark(marks, ray_index, i) != 0) continue;
            val = ray_val(sweepc->vals->ray[ray_index], i);
            if (is_missing(val, dp)) continue;
            if (fabs(val) <= dp->ckval) continue;
            
            /* find velocities in related sweeps */
            prevval = soundval = abval = dp->missingVal;
            if (sweepc->last!=NULL)
                prevval=ray_val(sweepc->last->ray[prevIndex], i);
            if (sweepc->sound != NULL && sweepc->last == NULL)
                soundval=ray_val(sweepc->sound->ray[ray_index], i);
            if (sweepc->above != NULL && valid_above_ray==0)
                abval=ray_val(sweepc->above->ray[abIndex], i);
            
            if (sweepc->last == NULL) {
                if (is_missing(soundval, dp)) continue;
                if (is_missing(abval, dp)) {
                    diff = min_diff(&val, soundval, sweepc->NyqVelocity, 
                                    dp->maxcount);
                    if (diff < thresh)
                        set_mark_ray(marks, sweepc->rv, ray_index, i, 1, val);
                } else {
                    diff = min_diff(&val, abval, sweepc->NyqVelocity,
                                    dp->maxcount);
                    if (diff < thresh && fabs(soundval-val) < thresh)
                        set_mark_ray(marks, sweepc->rv, ray_index, i, 1, val);
                }
            } else {
                if (not_missing(prevval, dp) && not_missing(abval, dp)) {
                    diff = min_diff(&val, prevval, sweepc->NyqVelocity,
                                    dp->maxcount);
                    if (diff < thresh && fabs(abval-val) < thresh)
                        set_mark_ray(marks, sweepc->rv, ray_index, i, 1, val); 
                }
            }
        }
    }
}


/* Mark gates where all neighbors are marked with a 0 as -1 */
static void mark_neighborless(
    Marks *marks, int ray_index, int i, int nrays, int nbins)
{
    int next, prev, left, right;
    int countbins = 0;
    if (ray_index == 0) left = nrays-1;
    else left = ray_index-1;
    if (ray_index == nrays-1) right = 0;
    else right = ray_index+1;
    next = i+1;
    prev = i-1;
 
    if (i != 0) {
        if (get_mark(marks, left, prev) == 0)
            countbins++;
        if (get_mark(marks, ray_index, prev) == 0)
            countbins++;
        if (get_mark(marks, right, prev) == 0)
            countbins++;
    }
    if (get_mark(marks, left, i) == 0)
        countbins++;
    if (get_mark(marks, right, i) == 0)
        countbins++;
    if (i < nbins-1) {
        if (get_mark(marks, left, next) == 0)
            countbins++;
        if (get_mark(marks, ray_index, next) == 0)
            countbins++;
        if (get_mark(marks, right, next) == 0)
            countbins++;
    }
    if (countbins == 0)
        set_mark(marks, ray_index, i, -1);
    return;
}

/* Record the location of a bin and ray index in the arrays */
static void record_gate(int idx, int b_idx, int r_idx, int b[8], int r[8])
{
    b[idx] = b_idx;
    r[idx] = r_idx;
    return;
}

/* Count the number of neighbors marked with a 1 */
static int count_neighbors(
    SweepCollection *sweepc, Marks *marks,
    int ray_index, int i, int binindex[8], int rayindex[8]
    )
{
    int next, prev, left, right;
    int countindex = 0;
    
    if (ray_index == 0) left = sweepc->nrays - 1;
    else left = ray_index - 1;
    if (ray_index == sweepc->nrays - 1) right=0;
    else right = ray_index + 1;
    next = i + 1;
    prev = i - 1;
 
    /* Look at all 8 bins adjacent to current bin in question: */
    if (i!=0) { 
        if (get_mark(marks, left, prev) == 1)
            record_gate(countindex++, prev, left, binindex, rayindex);
        if (get_mark(marks, ray_index, prev) == 1)
            record_gate(countindex++, prev, ray_index, binindex, rayindex);
        if (get_mark(marks, right, prev) == 1)
            record_gate(countindex++, prev, right, binindex, rayindex);
    }
    if (get_mark(marks, left, i) == 1)
        record_gate(countindex++, i, left, binindex, rayindex);
    if (get_mark(marks, right, i) == 1)
        record_gate(countindex++, i, right, binindex, rayindex);
    if (i < sweepc->nbins-1) {  
        if (get_mark(marks, left, next) == 1)
            record_gate(countindex++, next, left, binindex, rayindex);
        if (get_mark(marks, ray_index, next) == 1)
            record_gate(countindex++, next, ray_index, binindex, rayindex);
        if (get_mark(marks, right, next) == 1)
            record_gate(countindex++, next, right, binindex, rayindex);
    }
    return countindex;
}


/* Unfold against all adjacent values with marks of 1 */
static void unfold_adjacent(SweepCollection *sweepc, DealiasParams *dp,
    Marks *marks, int i, int ray_index, int countindex, int loopcount,
    int binindex[8], int rayindex[8])
{
    int in, out, numneg, numpos, l, numtimes;
    float diff, val;
    val = ray_val(sweepc->vals->ray[ray_index], i); 
    numtimes=0;
    while(not_missing(val, dp) && get_mark(marks, ray_index, i) == 0) {
        numtimes++;
        in = out = numneg = numpos = 0;
        for (l=0; l<countindex; l++) {
            diff = ray_val(sweepc->rv->ray[rayindex[l]], binindex[l]) - val;
            if (fabs(diff) < dp->thresh * sweepc->NyqVelocity)
                in++;
            else
                out++;
            if (diff > sweepc->NyqVelocity)
                numpos++;
            else if (diff < -sweepc->NyqVelocity)
                numneg++; 
        }
        if (in>0 && out==0) {
            set_mark_ray(marks, sweepc->rv, ray_index, i, 1, val);
        } else {
            if (numtimes <= dp->maxcount) {
                if ((numpos+numneg)<(in+out-(numpos+numneg))) {
                    if (loopcount>2) {
                        /* Keep the value after two passes through data. */
                        set_mark_ray(marks, sweepc->rv, ray_index, i, 1, val);
                    }
                } else if (numpos>numneg) {
                    val = val + sweepc->NyqInterval;
                } else if (numneg>numpos) {
                    val = val - sweepc->NyqInterval;
                } else {
                    /* Save bin for windowing if unsuccessful after
                    ** four passes: */
                    if (loopcount>4) 
                        set_mark(marks, ray_index, i, -2);
                }
            } else {
                /* Remove bin: */
                set_mark(marks, ray_index, i, -2);
            }
        }
    }
}

/* similar to unfold_adjacent but slighly different, used 
 * during second pass of dealiasing */
static void unfold_adjacent2(SweepCollection *sweepc, DealiasParams *dp,
    Marks *marks, int i, int ray_index, int countindex, int loopcount,
    int binindex[8], int rayindex[8]
    )
{
    int in, out, numneg, numpos, l, numtimes;
    float diff, val;
    val = ray_val(sweepc->vals->ray[ray_index], i); 
    numtimes=0;
    while(not_missing(val, dp) && get_mark(marks, ray_index, i) == 0) {
        numtimes++;
        in = out = numneg = numpos = 0;
        for (l=0; l<countindex; l++) {
            diff = ray_val(sweepc->rv->ray[rayindex[l]], binindex[l]) - val;
            if (fabs(diff) < dp->thresh * sweepc->NyqVelocity) {
                in++;
            } else {
                out++;
            }
            if (diff > sweepc->NyqVelocity) {
                numpos++;
            } else if (diff< -sweepc->NyqVelocity) {
                numneg++;
            }
        }
        if (in>out) {
            set_mark_ray(marks, sweepc->rv, ray_index, i, 1, val);
        } else {
            if (numtimes <= dp->maxcount) {
                if (numpos+numneg<(in+out-(numpos+numneg))) {
                    if (loopcount<=2) {
                        val = dp->missingVal; /* Try later */
                    } else {
                        /* Keep the value after two passes */
                        set_mark_ray(marks, sweepc->rv, ray_index, i, 1, val); 
                    }
                } else if (numpos>numneg) {
                    val = val + sweepc->NyqInterval;
                } else if (numneg>numpos) {
                    val = val - sweepc->NyqInterval;
                } else {
                    /* Remove bin after four passes through data: */
                    if (loopcount>4) 
                        set_mark(marks, ray_index, i, -1);
                }
            } else {
                /* Remove bin: */
                set_mark(marks, ray_index, i, -1);
            }
        }
    }
}

/* Unfold bins mark with 0 assuming spatial continuity: */
static void spatial_dealias(
    SweepCollection *sweepc, DealiasParams *dp, 
    Marks *marks, int *pstep, int pass)
{
    int loopcount, start, end, i, countindex, ray_index;
    int binindex[8], rayindex[8];
    unsigned short flag;

    loopcount=0;
    flag=1;
    while (flag==1) {
        loopcount++;
        flag=0;
        if (*pstep==1) {
            *pstep=-1;
            start=sweepc->nrays-1;
            end=-1;
        } else {
            *pstep=1;
            start=0;
            end=sweepc->nrays;
        }
        if (pass == 1) {
            for (i=0; i < sweepc->nbins; i++) {
                for (ray_index=start; ray_index != end; ray_index=ray_index+*pstep) {
                    if (get_mark(marks, ray_index, i) != 0)
                        continue;
                
                    countindex = count_neighbors(sweepc, marks, ray_index, i,
                                                 binindex, rayindex); 
                    if (countindex>=1) {
                        flag=1;
                        unfold_adjacent(sweepc, dp, marks,
                            i, ray_index, countindex, loopcount,
                            binindex, rayindex);
                    }
                
                    /* only performed on the first pass of the spatial
                     * dealiasing */
                    if (loopcount==1 && countindex<1) {
                    mark_neighborless(marks, ray_index, i, 
                                      sweepc->nrays, sweepc->nbins);
                    }
                }
            }
        } else {  /* 2nd pass, reverse ray/bin loop order, different unfold 
                     and no marking of neighborless gates */
            for (ray_index=start;ray_index!=end;ray_index=ray_index+*pstep) {
                for (i=0; i < sweepc->nbins; i++) {
                    if (get_mark(marks, ray_index, i) != 0)
                        continue;
                
                    countindex = count_neighbors(sweepc,
                        marks, ray_index, i, binindex, rayindex);
                
                    if (countindex>=1) {
                        flag=1;
                        unfold_adjacent2(sweepc, dp, marks,
                            i, ray_index, countindex, loopcount,
                            binindex, rayindex);
                    }
                }
            }
        }    
    }
}

/* Calculate the mean for gates within a window of a sweep */
static float window(
    DealiasParams *dp, Sweep *rv_sweep, int startray, int endray, 
    int firstbin, int lastbin, unsigned short* success) 
{
    int num, ray_index, rangeIndex, numRays, numBins;
    float val, sum, sumsq, ave, NyqVelocity, std;
     
    *success=0;
    NyqVelocity = rv_sweep->ray[0]->h.nyq_vel;
    numRays = rv_sweep->h.nrays;
    numBins = rv_sweep->ray[0]->h.nbins;

    /* Now, sum the data in the window region between startray, 
     **  endray, firstbin, lastbin. */
    num=0;
    sum=0.0;
    sumsq=0.0;
       
    if (firstbin>=numBins || lastbin>=numBins || firstbin<0 || lastbin<0)
        return dp->missingVal;
	for (rangeIndex=firstbin; rangeIndex<=lastbin; rangeIndex++) {
        if (startray>endray){
            for (ray_index=startray; ray_index<numRays; ray_index++) {
	            val = ray_val(rv_sweep->ray[ray_index], rangeIndex);
	            if (not_missing(val, dp)) {
	                num++;
	                sum+=val;
	                sumsq+=(val*val);
	            }
	        }
            for (ray_index=0; ray_index<=endray; ray_index++) {
	            val = ray_val(rv_sweep->ray[ray_index], rangeIndex);
	            if (not_missing(val, dp)) {
	                num++;
	                sum+=val;
	                sumsq+=(val*val);
	            }
	        }
        } else {
            for (ray_index=startray; ray_index<=endray; ray_index++) { 
	            val = ray_val(rv_sweep->ray[ray_index], rangeIndex);
                if (not_missing(val, dp)) {
	                num++;
	                sum+=val;
	                sumsq+=(val*val);
	            }
	        }
        }
    }
    if (num>=dp->mingood) {
        ave=sum/num;
        std=sqrt(fabs((sumsq-(sum*sum)/num)/(num-1)));
        if (std<= dp->stdthresh*NyqVelocity) *success=1;
    } else {
        ave = dp->missingVal;
        std=0.0;
        *success=1;
    }
    return ave; 
}

/* Detemine the window mean and if windowing is succeeded */
static float find_window(
    DealiasParams *dp, Sweep *rv_sweep, int ray_index, int i, 
    unsigned short *wsuccess)
{
    int startray, endray, firstbin, lastbin;
    float winval;
    int numRays = rv_sweep->h.nrays;
    int numBins = rv_sweep->ray[0]->h.nbins;

    startray = ray_index - dp->proximity;
    endray = ray_index + dp->proximity;
    firstbin = i - dp->proximity;
    lastbin = i + dp->proximity;
    if (startray < 0) startray = numRays + startray;
    if (endray > numRays - 1) endray = endray - numRays;
    if (firstbin < 0) firstbin = 0;
    if (lastbin > numBins-1) lastbin = numBins - 1;
    winval=window(dp, rv_sweep, startray, endray, 
                   firstbin, lastbin, wsuccess);
    if (is_missing(winval, dp) && *wsuccess == 1) {     
        /* Expand the window: */
        startray = ray_index - 2 * dp->proximity;
        endray = ray_index + 2 * dp->proximity;
        firstbin = i - 2 * dp->proximity;
        lastbin = i + 2 * dp->proximity;
        if (startray < 0) startray = numRays + startray;
        if (endray > numRays - 1) endray = endray - numRays;
        if (firstbin < 0) firstbin = 0;
        if (lastbin > numBins - 1) lastbin = numBins - 1;
        winval=window(dp, rv_sweep, startray, endray, 
                       firstbin, lastbin, wsuccess);
    }
    return winval;
}

/* Unfold remote bins or those that were previously unsuccessful
 * using a window with dimensions 2(PROXIMITY)+1 x 2(PROXIMITY)+1:
 * if still no luck delete data (or unfold against VAD if PASS2). */
static void unfold_remote(
    SweepCollection *sweepc, DealiasParams *dp, Marks *marks)
{
    int i, j, ray_index; 
    unsigned short wsuccess;
    float val, diff, winval;
    
    for (i=0; i < sweepc->nbins; i++) {
        for (ray_index=0; ray_index < sweepc->nrays; ray_index++) { 
            
            if (get_mark(marks, ray_index, i) !=0 && get_mark(marks, ray_index, i) != -2)
                continue;
            val = ray_val(sweepc->vals->ray[ray_index], i); 
            winval = find_window(dp, sweepc->rv, ray_index, i, &wsuccess);
            if (not_missing(winval, dp)) {
                
                /* Add or subtract Nyquist intervals until within range */
                for (j=0; j < dp->maxcount; j++) {
                    if (fabs(winval - val) <= 0.99999 * sweepc->NyqVelocity)
                        break;
                    if (val > winval)
                        val -= sweepc->NyqInterval;
                    else
                        val += sweepc->NyqInterval;
                }
                diff = fabs(winval - val);
                
                if (diff < dp->thresh * sweepc->NyqVelocity) {
                    /* Return the value. */
                    set_mark_ray(marks, sweepc->rv, ray_index, i, 1, val);
                } else if (diff<(1.0 - (1.0 - dp->thresh)/2.0)*sweepc->NyqVelocity) {
                    /* If within relaxed threshold, then return value, but
                    **   do not use to dealias other bins. */
                    set_mark_ray(marks, sweepc->rv, ray_index, i, -1, val);
                } else {
                    /* Remove bin */
                    set_mark(marks, ray_index, i, -1);
                }
            } else {
                if (wsuccess==0) {
                    /* Remove bin */
                    set_mark(marks, ray_index, i, -1);
                    /* I don't think this ever happens -jjh */
                } else if (sweepc->sound==NULL || sweepc->last==NULL) {
                    if (get_mark(marks, ray_index, i) == 0 && dp->rm != 1) {
                        /* Leave bin untouched. */
                        /* Don't use to unfold other bins*/
                        val = ray_val(sweepc->vals->ray[ray_index], i);
                        set_mark_ray(marks, sweepc->rv, ray_index, i, -1, val);
                    } else {
                        /* Remove bin */
                        set_mark(marks, ray_index, i, -1);
                    }
                } else if (get_mark(marks, ray_index, i) == 0 && dp->pass2 &&
                           sweepc->sound != NULL && sweepc->last != NULL) {
                    /* Leave bins marked with a 0 for a second pass. */
                } else {
                    /* Remove bin */
                    set_mark(marks, ray_index, i, -1);
                }
            }
        }
    }
}

/* Second pass dealiasing using only sounding data */
/* This is similar to continuity_dealias */
static void second_pass(
    SweepCollection *sweepc, DealiasParams *dp, Marks *marks)
{
    int i, ray_index; 
    float val, diff, soundval;
    
    for (ray_index=0; ray_index < sweepc->nrays; ray_index++) {
        for (i=0; i < sweepc->nbins; i++) {
    
            if (get_mark(marks, ray_index, i) != 0) continue;
            val = ray_val(sweepc->vals->ray[ray_index], i);
            if (is_missing(val, dp)) continue; 
            if (fabs(val) <= dp->ckval) continue;
            
            soundval = ray_val(sweepc->sound->ray[ray_index], i);
            if (is_missing(soundval, dp)) continue; 

            diff = min_diff(&val, soundval, sweepc->NyqVelocity, dp->maxcount); 
            if (diff < dp->compthresh2 * sweepc->NyqVelocity)
                set_mark_ray(marks, sweepc->rv, ray_index, i, 1, val);
        }
    }
}

/********************
 * Public functions *
*********************/

/*
**     This algorithm unfolds a volume of single Doppler radial velocity data.
**  The algorithm uses a previously unfolded volume (or VAD if previous volume
**  is unavailable) and the previous elevation sweep to unfold some of gates 
**  in each sweep. Then, it spreads outward from the 'good' gates, completing
**  the unfolding using gate-to-gate continuity. Gates that still remain
**  unfolded are compared to an areal average of neighboring dealiased gates.
**  Isolated echoes that still remain uncorrected are dealiased against a VAD
**  (as a last resort).
**
** This routine performs a preliminary unfold on the data using the bin in the
** next highest sweep (already unfolded), and the last unfolded volume. If the 
** last unfolded volume is not available, the VAD is used (or sounding if VAD
** not available). If this is successful, the bin is considered GOOD. Every 
** other bin is set to GOOD=0 or GOOD=-1 (if the bin is missing or bad).
**
** Then, the algorithm scans azimuthally and radially. If the majority of the
** GOOD=1 bins (up to 8) adjacent to a particular GOOD=0 bin are within 
** THRESH*Vnyq it is saved and considered GOOD as well. Otherwise, Nyquist 
** intervals are added/subtracted until the condition is met. If unsuccessful, 
** then GOOD is set to -2 for that bin.
**
** When all GOOD=0 bins that are next to GOOD=1 bins have been examined, then
** GOOD=0 and GOOD=-2 bins are unfolded against a window of GOOD=1 values. If
** there are too few GOOD=1 values inside the window, then the data are
** unfolded against the VAD, if available.
** 
** PARAMETERS:
** rvVolume: The radial velocity field to be unfolded.
** soundVolume: A first guess radial velocity field using VAD (or sounding) 
** lastVolume: The last radial velocity field unfolded (time of lastVolume 
**    should be within 10 minutes prior to rvVolume)
** missingVal: The value for missing radial velocity data.
** filt: A flag that specifies whether or not to use Bergen/Albers filter.
** success: flag indicating whether or not unfolding is possible.
*/
int dealias_fourdd(
    Volume* rvVolume, Volume* soundVolume, Volume* lastVolume,
    float missingVal, float compthresh, float compthresh2, float thresh,
    float ckval, float stdthresh, float epsilon,
    int maxcount, int pass2, int rm, int proximity, int mingood, int filt,
    int ba_mincount, int ba_edgecount)
{
    int sweepIndex, numSweeps, maxrays, maxbins;
    int step = -1;
    short *mark_data;
    Volume* VALS;
    SweepCollection sweepc;
    DealiasParams dp; 
    Marks marks;

    /* Either a sounding or last volume must be provided */
    if (soundVolume==NULL && lastVolume==NULL) {
        printf("First guess not available.\n");
        return 0;
    }

    /* Find the maximum number of rays and bins in all sweeps */
    numSweeps = rvVolume->h.nsweeps;
    maxrays = maxbins = 0;
    for (sweepIndex=numSweeps-1; sweepIndex>=0; sweepIndex--) {
        if (rvVolume->sweep[sweepIndex]->h.nrays > maxrays)
            maxrays = rvVolume->sweep[sweepIndex]->h.nrays;
        if (rvVolume->sweep[sweepIndex]->ray[0]->h.nbins > maxbins)
            maxbins = rvVolume->sweep[sweepIndex]->ray[0]->h.nbins;
    }
    
    /* Populate marks */
    marks.maxrays = maxrays;
    marks.maxbins = maxbins;
    if (NULL == (mark_data = malloc(maxbins * maxrays * sizeof(short)))) {
        printf("MemoryError: Error allocating memory for mark_data\n");
        return 0;
    }
    marks.mark_data = mark_data;

    /* Fill in Dealiasing parameters from arguments */
    dp.missingVal = missingVal;
    dp.epsilon = epsilon;
    dp.compthresh = compthresh;
    dp.compthresh2 = compthresh2;
    dp.thresh = thresh;
    dp.ckval = ckval;
    dp.stdthresh = stdthresh;
    dp.maxcount = maxcount;
    dp.pass2 = pass2;
    dp.rm = rm;
    dp.proximity = proximity;
    dp.mingood = mingood;
    dp.ba_mincount = ba_mincount;
    dp.ba_edgecount = ba_edgecount;

    /* Loop over sweeps, unfolding each one */
    VALS=RSL_copy_volume(rvVolume); 
    for (sweepIndex=numSweeps-1; sweepIndex>=0; sweepIndex--) {
        
        /* load sweeps and sweep structure information */
        sweepc.last = sweepc.sound = sweepc.above = NULL;
        if (lastVolume!=NULL) 
            sweepc.last = lastVolume->sweep[sweepIndex];
        if (soundVolume!=NULL) 
            sweepc.sound = soundVolume->sweep[sweepIndex];
        if (sweepIndex < numSweeps-1) 
            sweepc.above = rvVolume->sweep[sweepIndex+1];
        sweepc.rv = rvVolume->sweep[sweepIndex];
        sweepc.vals = VALS->sweep[sweepIndex];
        sweepc.NyqVelocity =rvVolume->sweep[sweepIndex]->ray[0]->h.nyq_vel;
        sweepc.NyqInterval = 2.0 * sweepc.NyqVelocity;
        sweepc.nrays = sweepc.rv->h.nrays;
        sweepc.nbins = sweepc.rv->ray[0]->h.nbins;

        /* Mark gates with -1 missing value or 0 where not missing.
         * If requested the Berger Albers 3x3 filter is applied. */
        if (filt == 1)
            berger_albers_filter(sweepc.vals, &dp, &marks);
        else
            init_marks(sweepc.vals, &dp, &marks);
    
        /* Initialize the output sweep with missingVal */
        init_sweep(sweepc.rv, dp.missingVal);

        /* Find non-missing value gates where the velocity either agrees with
         * or can be unfolded to agree with the sounding, the last volume, 
         * or the sweep directly above the current sweep. Record the unfolded
         * velocity in these cases and mark with a 1. */
        continuity_dealias(&sweepc, &dp, &marks);
        
        spatial_dealias(&sweepc, &dp, &marks, &step, 1);
        
        unfold_remote(&sweepc, &dp, &marks);

        if (sweepc.last != NULL && sweepc.sound != NULL) {
            second_pass(&sweepc, &dp, &marks); 
            spatial_dealias(&sweepc, &dp, &marks, &step, 2);
        }
    } // end of loop over sweeps
    free(mark_data);
    return 1;
}
