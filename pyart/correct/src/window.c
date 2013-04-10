/*
**
**  UW Radial Velocity Dealiasing Algorithm
**  Four-Dimensional Dealiasing (4DD)
**
**  DESCRIPTION:
**      This routine averages the values in a range and azimuth window of a
**      sweep and computes the standard deviation.
**
**  DEVELOPER:
**	Curtis N. James    26 Jan 1999
**
**
**
*/
#include "FourDD.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <rsl.h> /* Sweep */

float window(Volume* rvVolume, int sweepIndex, int startray, 
	int endray, int firstbin, int lastbin, float* std, float
        missingVal, unsigned short* success) {

     int num, currIndex, rangeIndex, numRays, numBins;
     float val, sum, sumsq, ave, NyqVelocity;
     
     *success=0;
     NyqVelocity = rvVolume->sweep[sweepIndex]->ray[0]->
	 h.nyq_vel;
     numRays = rvVolume->sweep[sweepIndex]->h.nrays;
     numBins = rvVolume->sweep[sweepIndex]->ray[0]->h.nbins;

     /* Now, sum the data in the window region between startray, 
     **  endray, firstbin, lastbin. */
     *std=0.0;
     ave=0.0;
     num=0;
     sum=0.0;
     sumsq=0.0;
       
     if (firstbin>=numBins || lastbin>=numBins || firstbin<0 || lastbin<0)
       return missingVal;
     if (startray>endray){
       for (currIndex=startray; currIndex<numRays; currIndex++) {
	 for (rangeIndex=firstbin; rangeIndex<=lastbin; rangeIndex++) {
	   val=(float)rvVolume->sweep[sweepIndex]->ray[currIndex]->
	     h.f(rvVolume->sweep[sweepIndex]->ray[currIndex]->
		 range[rangeIndex]);
	   if (val!=missingVal) {
	     num=num+1;
	     sum=sum+val;
	     sumsq=sumsq+val*val;
	   }
	 }
       }
       for (currIndex=0; currIndex<=endray; currIndex++) {
	 for (rangeIndex=firstbin; rangeIndex<=lastbin; rangeIndex++) {
	   val=(float)rvVolume->sweep[sweepIndex]->ray[currIndex]->
	     h.f(rvVolume->sweep[sweepIndex]->ray[currIndex]->
		 range[rangeIndex]);
	   if (val!=missingVal) {
	     num=num+1;
	     sum=sum+val;
	     sumsq=sumsq+val*val;
	   }
	 }
       }
     } else {
       for (currIndex=startray; currIndex<=endray; currIndex++) { 
	 for (rangeIndex=firstbin; rangeIndex<=lastbin; rangeIndex++) {
	   val=(float)rvVolume->sweep[sweepIndex]->ray[currIndex]->
	     h.f(rvVolume->sweep[sweepIndex]->ray[currIndex]->
		 range[rangeIndex]);
	   if (val!=missingVal) {
	     num=num+1;
	     sum=sum+val;
	     sumsq=sumsq+val*val;
	   }
	 }
       }
     }
     if (num>=MINGOOD) {
       ave=sum/num;
       *std=sqrt(fabs((sumsq-(sum*sum)/num)/(num-1)));
       if (*std<=STDTHRESH*NyqVelocity) *success=1;
       /* printf("ave=%0.2f, std=%0.2f, sum=%0.2f\n", ave, *std, sum); */
     } else {
       ave=missingVal;
       *std=0.0;
       *success=1;
     }
     return ave; 
}

