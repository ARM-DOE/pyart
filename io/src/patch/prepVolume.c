/*
**
**  UW Radial Velocity Dealiasing Algorithm
**  Four-Dimensional Dealiasing (4DD)
**
**  DESCRIPTION:
**      This routine thresholds VR bins where reflectivity is below
**      LOWDBZ or missing (if NODBZRMRV==1).
**
**  DEVELOPER:
**	Curtis N. James 20 Mar 98
**      Modified by C. James 25 Jan 99
**
**
*/
#include "FourDD.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <rsl.h> /* Sweep */

void testme(){
	printf("chho hoo");
}

void prepVolume(Volume* DBZVolume, Volume* rvVolume, float missingVal) {

     int currIndex, sweepIndex, i, j, DBZIndex, numRays, numBins, numDBZRays,
       numDBZBins, numSweepsRV, numSweepsDZ;
     float val, dzval, finalval, minpossDBZ = -50.0, startaz;
     int DBZGateSize = 0, DBZRangeBin1 = 0, DBZfactor, step, gap, limit;
     int rvGateSize = 0, rvRangeBin1 = 0, remainder;
     div_t ratioGateSize;
     
     numSweepsRV = rvVolume->h.nsweeps;
     numSweepsDZ = DBZVolume->h.nsweeps;
     if (numSweepsRV!=numSweepsDZ) {
       return;
     } else {
       for (sweepIndex=0;sweepIndex<numSweepsRV;sweepIndex++) {
	 numRays = rvVolume->sweep[sweepIndex]->h.nrays;
	 numBins = rvVolume->sweep[sweepIndex]->ray[0]->h.nbins;
	 if (DBZVolume!=NULL) {
	   numDBZRays = DBZVolume->sweep[sweepIndex]->h.nrays;
	   numDBZBins = DBZVolume->sweep[sweepIndex]->ray[0]->h.nbins;
	 }
     
	 /* Get the bin geometry for both DBZ and rv: */

	 if (VERBOSE) printf(" ...preparing sweep...\n");
	 if (DBZVolume!=NULL) {
	   DBZGateSize=DBZVolume->sweep[sweepIndex]->ray[0]->
	     h.gate_size;
	   DBZRangeBin1=DBZVolume->sweep[sweepIndex]->ray[0]->
	     h.range_bin1;
	 }
	 rvGateSize=rvVolume->sweep[sweepIndex]->ray[0]->
	   h.gate_size;
	 rvRangeBin1=rvVolume->sweep[sweepIndex]->ray[0]->
	   h.range_bin1;

	 /* Terminate routine if the locations of the first range bins
	 ** in DBZ and rv do not coincide. */
       
	 if ((DBZRangeBin1!=rvRangeBin1 || numDBZRays!=numRays) &&
             DBZVolume!=NULL) {
	   return;
	 }
	 /* Compare the bin geometry between DBZ and rv to know which
	 ** DBZ bin(s) corresponds spatially to the rv bin(s)*/
       
	 if (DBZVolume!=NULL) {
	   ratioGateSize=div(rvGateSize,DBZGateSize);
	   remainder=ratioGateSize.rem;
	   if (remainder!=0) {
	     return;
	   } else if (remainder==0 && rvGateSize!=0 && DBZGateSize !=0) {
	     DBZfactor=ratioGateSize.quot;
	   } else {
	     return;
	   }
	 }

	 /* Erase the first DELNUM bins of DBZ and radial velocity, as per
	 **  communication with Peter Hildebrand. */
	 for (currIndex=0;currIndex<numRays;currIndex++) {
	   for (i = 0; i < DELNUM; i++) {
	     rvVolume->sweep[sweepIndex]->ray[currIndex]->range[i]=
	       (unsigned short) (missingVal);
	     if (DBZVolume!=NULL) {
	       for (j=0; j<DBZfactor; j++) {	       
		 DBZIndex=i*DBZfactor+j;
		 DBZVolume->sweep[sweepIndex]->ray[currIndex]->range[DBZIndex]=
		   (unsigned short) (missingVal);
	       }
	     }
	   }
	 }
	 
	 /* Now that we know the bin geometry, we assign missingVal to any
	 ** velocity bins where the reflectivity is missing (if NODBZRMRV==1)
	 ** or outside LOWDBZ and HIGHDBZ. */
	 

	 limit=numBins;
	 for (currIndex=0; currIndex<numRays; currIndex++) {
	   for (i=DELNUM; i<limit; i++) {
 	     if (DBZVolume!=NULL) {
	       val=(float) rvVolume->sweep[sweepIndex]->ray[currIndex]
		 ->h.f(rvVolume->sweep[sweepIndex]->ray[currIndex]->
		      range[i]);
	       for (j=0; j<DBZfactor; j++) {	       
		 DBZIndex=i*DBZfactor+j;
		 if (DBZIndex>=numDBZBins) {
		   printf("Invalid bin geometry!\n");
		   return;
		 }
		 dzval=(float) DBZVolume->sweep[sweepIndex]->ray[currIndex]
		   ->h.f(DBZVolume->sweep[sweepIndex]->ray[currIndex]->
			 range[DBZIndex]);
		 if ((dzval>=minpossDBZ && dzval<LOWDBZ)||(((dzval<
		     minpossDBZ)||(dzval>HIGHDBZ))&&NODBZRMRV==1)) {
		    rvVolume->sweep[sweepIndex]->ray[currIndex]->
		      range[i]=(unsigned short) (missingVal);
		    val=missingVal;
		 } 
	       }
	     }
	   }
	 }
       }
     }
     if (VERBOSE) printf(" ...volume prepped.\n");      
     return;
}
