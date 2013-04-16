/*
**
**  UW Radial Velocity Dealiasing Algorithm
**  Four-Dimensional Dealiasing (4DD)
**
**  DESCRIPTION:
**      This routine finds the rayindex of the nearest ray in sweepIndex2 of 
**  rvVolume2 to sweepIndex1 in rvVolume1.
**
**  DEVELOPER:
**	Curtis N. James    1 Feb 1999
**
**
**
*/
#include "FourDD.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <rsl.h>

int findRay (Volume* rvVolume1, Volume* rvVolume2, int sweepIndex1, int
     sweepIndex2, int currIndex, float missingVal) {

     int numSweeps, numRays, numBins, gatesz0, gatesz1, rayIndex1, i1;
     float prevval, az0, az1, diffaz, range0, range1;
     float startrange0,startrange1,spacing;
     short direction, lastdir;
     
     numRays = rvVolume2->sweep[sweepIndex2]->h.nrays;

     az0 = rvVolume1->sweep[sweepIndex1]->ray[currIndex]->h.azimuth;
     if (currIndex<numRays) rayIndex1=currIndex;
     else rayIndex1=numRays-1;
     az1 = rvVolume2->sweep[sweepIndex2]->ray[rayIndex1]->h.azimuth;
     if (az0==az1) {
       return rayIndex1;
     } else {
       /* Since the beamwidth is not necessarily the spacing between rays: */
       spacing = fabs(rvVolume2->sweep[sweepIndex2]->ray[0]->h.azimuth-
		      rvVolume2->sweep[sweepIndex2]->ray[50]->h.azimuth); 
       if (spacing>180) spacing=360.0-spacing;
       spacing=spacing/50.0;

       /* Compute the difference in azimuth between the two rays: */
       diffaz=az0-az1;
       if (diffaz>=180.0) diffaz=diffaz-360.0;
       else if (diffaz<-180.0) diffaz=diffaz+360.0;
       
       /* Get close to the correct index: */
       rayIndex1=rayIndex1+(int) (diffaz/spacing);
       if (rayIndex1>=numRays) rayIndex1=rayIndex1-numRays;
       if (rayIndex1<0) rayIndex1=numRays+rayIndex1;
       az1=rvVolume2->sweep[sweepIndex2]->ray[rayIndex1]->h.azimuth;
       diffaz=az0-az1;
       if (diffaz>=180.0) diffaz=diffaz-360.0;
       else if (diffaz<-180.0) diffaz=diffaz+360.0;

       /* Now add or subtract indices until the nearest ray is found: */
       if (diffaz>=0) lastdir=1;
       else lastdir=-1;
       while (fabs(diffaz)>spacing/2.0) {
	 if (diffaz>=0) {
	   rayIndex1++;
	   direction=1;
	 } else {
	   rayIndex1--;
	   direction=-1;
	 }
	 if (rayIndex1>=numRays) rayIndex1=rayIndex1-numRays;
	 if (rayIndex1<0) rayIndex1=numRays+rayIndex1;
	 az1=rvVolume2->sweep[sweepIndex2]->ray[rayIndex1]->h.azimuth;
	 diffaz=az0-az1;
	 if (diffaz>=180.0) diffaz=diffaz-360.0;
	 else if (diffaz<-180.0) diffaz=diffaz+360.0;
	 if (direction!=lastdir) break;
	 else lastdir=direction;
       }
       return rayIndex1;
     }
}

