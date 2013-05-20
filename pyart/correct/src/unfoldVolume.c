/*
**
**  UW Radial Velocity Dealiasing Algorithm
**  Four-Dimensional Dealiasing (4DD)
**
**  DESCRIPTION:
**     This algorithm unfolds a volume of single Doppler radial velocity data.
**  The algorithm uses a previously unfolded volume (or VAD if previous volume
**  is unavailable) and the previous elevation sweep to unfold some of gates 
**  in each sweep. Then, it spreads outward from the 'good' gates, completing
**  the unfolding using gate-to-gate continuity. Gates that still remain
**  unfolded are compared to an areal average of neighboring dealiased gates.
**  Isolated echoes that still remain uncorrected are dealiased against a VAD
**  (as a last resort).
**
**  DEVELOPER:
**	Curtis N. James     25 Jan 99
**
**
*/
#include "FourDD.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <rsl.h> /* Sweep */

/* This routine performs a preliminary unfold on the data using the bin in the
** next highest sweep (already unfolded), and the last unfolded volume. If the ** last unfolded volume is not available, the VAD is used (or sounding if VAD
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


void unfoldVolume(Volume* rvVolume, Volume* soundVolume, Volume* lastVolume,
      float missingVal, unsigned short filt, unsigned short* success) {

     int sweepIndex, currIndex, i, l, m, n, direction, numSweeps, numRays,
       numBins, numdo, left, right, next, prev, rayindex[8], binindex[8],
       countindex, numneg, numpos, in, out, startray, endray, firstbin, 
       lastbin, step = -1, startindex, endindex, prevIndex, abIndex, loopcount,
       countbins;
    
     unsigned short numtimes, dcase, flag=1, wsuccess;
     short GOOD[MAXBINS][MAXRAYS];
     float NyqVelocity, NyqInterval, val, diff, fraction, finalval, initval, 
       valcheck, goodval, winval, vdiff, diffs[8], fraction2;
     float prevval, abval, pfraction, backval, cval, soundval, std;
     Volume* VALS;

     numSweeps = rvVolume->h.nsweeps;
     VALS=RSL_copy_volume(rvVolume);
     if (COMPTHRESH>1.0 || COMPTHRESH<=0.0 ) fraction=0.25;
     else fraction=COMPTHRESH;
     if (COMPTHRESH2>1.0 || COMPTHRESH2<=0.0 ) fraction2=0.25;
     else fraction2=COMPTHRESH2;
     if (THRESH>1.0 || THRESH<=0.0 ) pfraction=0.5;
     else pfraction=THRESH;

     if (soundVolume!=NULL || lastVolume!=NULL) {
       for (sweepIndex=numSweeps-1;sweepIndex>=0;sweepIndex--) {
         NyqVelocity =rvVolume->sweep[sweepIndex]->ray[0]->h.nyq_vel;
	 if (NyqVelocity == 0.0) NyqVelocity=9.8;
         NyqInterval = 2.0 * NyqVelocity;
         printf("NyqVelocity %f \n", NyqVelocity);
	 numRays = rvVolume->sweep[sweepIndex]->h.nrays;
	 numBins = rvVolume->sweep[sweepIndex]->ray[0]->h.nbins;

	 /* First, unfold bins where vertical and temporal continuity
	 ** produces the same value (i.e., bins where the radial velocity
	 ** agrees with both the previous sweep and previous volume within
	 ** a COMPTHRESH of the Nyquist velocity). This produces a number
	 ** of good data points from which to begin unfolding assuming spatial
         ** continuity. */
	 printf("Sweep: %d\n", sweepIndex);
	 if (VERBOSE) printf("NyqVelocity: %f, missingVal: %f\n",
			     NyqVelocity, missingVal);
	 flag=1;

	 for (currIndex=0;currIndex<numRays;currIndex++) {
	   if (lastVolume!=NULL) prevIndex=findRay(rvVolume, lastVolume,
	      sweepIndex, sweepIndex,currIndex, missingVal);
	   if (sweepIndex<numSweeps-1) abIndex=findRay(rvVolume, rvVolume,
	      sweepIndex, sweepIndex+1, currIndex, missingVal);
 	   for (i=DELNUM;i<numBins;i++) {
	     /* Initialize Output Sweep with missing values: */
	     initval=(float)rvVolume->sweep[sweepIndex]->
	       ray[currIndex]->h.invf(missingVal);
	     rvVolume->sweep[sweepIndex]->ray[currIndex]->
	       range[i]=(unsigned short) (initval);

	     /* Assign uncorrect velocity bins to the array VALS: */
	     val=(float) VALS->sweep[sweepIndex]->
                 ray[currIndex]->h.f(VALS->sweep[sweepIndex]->ray
                 [currIndex]->range[i]);
	     valcheck=val;
	     if (val==missingVal) GOOD[i][currIndex]=-1;
	     else {
	       if (filt==1) {
		 /*Perform a 3x3 filter, as proposed by Bergen & Albers 1988*/
		 countindex=0;
		 if (currIndex==0) left=numRays-1;
		 else left=currIndex-1;
		 if (currIndex==numRays-1) right=0;
		 else right=currIndex+1;
		 next=i+1;
		 prev=i-1;
		 /* Look at all bins adjacent to current bin in question: */
		 if (i>DELNUM) {
		   if (((float) VALS->sweep[sweepIndex]->
			ray[left]->h.f(VALS->sweep[sweepIndex]->ray
	       	        [left]->range[prev]))!=missingVal) {
		     countindex=countindex+1;
		   }
		   if (((float) VALS->sweep[sweepIndex]->
			ray[currIndex]->h.f(VALS->sweep[sweepIndex]->ray
		        [currIndex]->range[prev]))!=missingVal) {
		     countindex=countindex+1;
		   }
		   if (((float) VALS->sweep[sweepIndex]->
			ray[right]->h.f(VALS->sweep[sweepIndex]->ray
			[right]->range[prev]))!=missingVal) {
		     countindex=countindex+1;
		   }
		 }
		 if (((float) VALS->sweep[sweepIndex]->
		      ray[left]->h.f(VALS->sweep[sweepIndex]->ray
		      [left]->range[i]))!=missingVal) {
		   countindex=countindex+1;
		 }
		 if (((float) VALS->sweep[sweepIndex]->
		      ray[right]->h.f(VALS->sweep[sweepIndex]->ray
		      [right]->range[i]))!=missingVal) {
		   countindex=countindex+1;
		 }
		 if (i<numBins-1) {  
		   if (((float) VALS->sweep[sweepIndex]->
			ray[left]->h.f(VALS->sweep[sweepIndex]->ray
                        [left]->range[next]))!=missingVal) {
		     countindex=countindex+1;
		   }
		   if (((float) VALS->sweep[sweepIndex]->
			ray[currIndex]->h.f(VALS->sweep[sweepIndex]->ray
			[currIndex]->range[next]))!=missingVal) {
		     countindex=countindex+1;
		   }
		   if (((float) VALS->sweep[sweepIndex]->
			ray[right]->h.f(VALS->sweep[sweepIndex]->ray
                        [right]->range[next]))!=missingVal) {
		     countindex=countindex+1;
		   }
		 }
		 if (((i==numBins-1 || i==DELNUM) && countindex>=3)||
		     (countindex>=5)) {
		   /* Save the bin for dealiasing: */
		   GOOD[i][currIndex]=0;
		 } else {
		   /* Assign missing value to the current bin. */
		   GOOD[i][currIndex]=-1;
		 }   
	       /* End 3 x 3 filter */
	       } else {
		 /* If no filter is being applied save bin for dealiasing: */
		 GOOD[i][currIndex]=0;
	       }
	     }
	     
	     if (GOOD[i][currIndex]==0) {
	       /* Try to dealias the bin using vertical and temporal 
	       **   continuity (this is initial dealiasing). */
	       if (val!=missingVal && lastVolume!=NULL) {
		 prevval=(float) lastVolume->sweep[sweepIndex]->ray[prevIndex]
		   ->h.f(lastVolume->sweep[sweepIndex]->ray[prevIndex]->
		   range[i]);
	       } else {
		 prevval=missingVal;
	       }
	       if (val!=missingVal && soundVolume!=NULL && lastVolume==NULL) {
		 soundval=(float) soundVolume->
		   sweep[sweepIndex]->ray[currIndex]->h.f(soundVolume->
		   sweep[sweepIndex]->ray[currIndex]->range[i]);
	       } else {
		 soundval=missingVal;
	       }
	       if (val!=missingVal && sweepIndex<numSweeps-1) {
		 abval=(float) rvVolume->sweep[sweepIndex+1]->ray[abIndex]->
		   h.f(rvVolume->sweep[sweepIndex+1]->ray[abIndex]->range[i]);
	       } else {
		 abval=missingVal;
	       }

	       if (val!=missingVal&&lastVolume==NULL&&soundval!=
		   missingVal&&abval==missingVal) {
		 cval=soundval;
		 dcase=1;
	       }
	       else if (val!=missingVal&&lastVolume==NULL&&soundval!=
			missingVal&&abval!=missingVal) {
		 cval=abval;
		 dcase=2;
	       }
	       else if (val!=missingVal&&prevval!=missingVal&&abval!=
			missingVal) {
		 cval=prevval;
		 dcase=3;
	       } 
	       else {
		 cval=missingVal;
		 dcase=0;
	       }
	       if (dcase>0) {
		 diff=cval-val;	     
		 if (diff<0.0) {
		   diff=-diff;
		   direction=-1;
		 } else direction=1;
		 numtimes=0;
		 while (diff>0.99999*NyqVelocity && numtimes<=MAXCOUNT) {
		   val=val+NyqInterval*direction;
		   numtimes=numtimes+1;
		   diff=cval-val;
		   if (diff<0.0) {
		     diff=-diff;
		     direction=-1;
		   } else direction=1;
		 }
		 if (dcase==1&&diff<fraction*NyqVelocity&&
		     fabs(valcheck)>CKVAL) {
		   if (VERBOSE) printf("GOOD1: %f\n", val);
		   finalval=(float)rvVolume->sweep[sweepIndex]->
		     ray[currIndex]->h.invf(val);
		   rvVolume->sweep[sweepIndex]->ray[currIndex]->
		     range[i]=(unsigned short) (finalval);
		   GOOD[i][currIndex]=1;
		 }
		 else if (dcase==2&&diff<fraction*NyqVelocity&&fabs(soundval
		   -val)<fraction*NyqVelocity&&fabs(valcheck)>CKVAL) {
		   if (VERBOSE) printf("GOOD2: %f\n", val);
		   finalval=(float)rvVolume->sweep[sweepIndex]->
		     ray[currIndex]->h.invf(val);
		   rvVolume->sweep[sweepIndex]->ray[currIndex]->
		     range[i]=(unsigned short) (finalval);
		   GOOD[i][currIndex]=1;
		 }
		 else if (dcase==3&&diff<fraction*NyqVelocity&&fabs(abval-val)<
			  fraction*NyqVelocity&&fabs(valcheck)>CKVAL) {
		   if (VERBOSE) printf("GOOD3: %f\n",val);
		   finalval=(float)rvVolume->sweep[sweepIndex]->
		     ray[currIndex]->h.invf(val);
		   rvVolume->sweep[sweepIndex]->ray[currIndex]->
		     range[i]=(unsigned short) (finalval); 
		   GOOD[i][currIndex]=1;
		 }
	       }
	     }
	   }
	 }

	 /* Now, unfold GOOD=0 bins assuming spatial continuity: */
	 loopcount=0;
	 while (flag==1) {
	   loopcount=loopcount+1;
	   flag=0;
	   if (step==1) {
	     step=-1;
	     startindex=numRays-1;
	     endindex=-1;
	   } else {
	     step=1;
	     startindex=0;
	     endindex=numRays;
	   }
	   for (i=DELNUM;i<numBins;i++) {
	     for (currIndex=startindex;currIndex!=endindex;currIndex=currIndex
		+step) {
	       if (GOOD[i][currIndex]==0) {
		 countindex=0;
		 countbins=0;
		 val=(float) VALS->sweep[sweepIndex]->ray[currIndex]->
		   h.f(VALS->sweep[sweepIndex]->ray[currIndex]->range[i]);
		 if (VERBOSE) printf("\n");
	         if (VERBOSE) printf("Startval: %f\n", val);
       		 if (currIndex==0) left=numRays-1;
		 else left=currIndex-1;
		 if (currIndex==numRays-1) right=0;
		 else right=currIndex+1;
		 next=i+1;
		 prev=i-1;
		 /* Look at all bins adjacent to current bin in question: */
		 if (val != missingVal) {  // don't look for neighbors of missingVal
         if (i>DELNUM) {
		   if (GOOD[prev][left]==1) {
		     if (flag==0) flag=1;
		     countindex=countindex+1;
		     binindex[countindex-1]=prev;
		     rayindex[countindex-1]=left;
		     if (VERBOSE) printf("pl ");
		   }
		   if (GOOD[prev][left]==0) {
		     countbins=countbins+1;
		   }
		   if (GOOD[prev][currIndex]==1) {
		     if (flag==0) flag=1;
		     countindex=countindex+1;
		     binindex[countindex-1]=prev;
		     rayindex[countindex-1]=currIndex;
		     if (VERBOSE) printf("pc ");
		   }
		   if (GOOD[prev][currIndex]==0) {
		     countbins=countbins+1;
		   }
		   if (GOOD[prev][right]==1) {
		     if (flag==0) flag=1;
		     countindex=countindex+1;
		     binindex[countindex-1]=prev;
		     rayindex[countindex-1]=right;
		     if (VERBOSE) printf("pr ");
		   }
	           if (GOOD[prev][right]==0) {
		     countbins=countbins+1;
		   }
		 }
		 if (GOOD[i][left]==1) {
		   if (flag==0) flag=1;
		   countindex=countindex+1;
		   binindex[countindex-1]=i;
		   rayindex[countindex-1]=left;
		   if (VERBOSE) printf("il ");
		 }
		 if (GOOD[i][left]==0) {
		   countbins=countbins+1;
		 }
		 if (GOOD[i][right]==1) {
		   if (flag==0) flag=1;
		   countindex=countindex+1;
		   binindex[countindex-1]=i;
		   rayindex[countindex-1]=right;
		   if (VERBOSE) printf("ir "); 	   
		 }
		 if (GOOD[i][right]==0) {
		   countbins=countbins+1;
		 }
		 if (i<numBins-1) {  
		   if (GOOD[next][left]==1) {
		     if (flag==0) flag=1;
		     countindex=countindex+1;
		     binindex[countindex-1]=next;
		     rayindex[countindex-1]=left;
		     if (VERBOSE) printf("nl "); 
		   }
		   if (GOOD[next][left]==0) {
		     countbins=countbins+1;
		   }
		   if (GOOD[next][currIndex]==1) {
		     if (flag==0) flag=1;
		     countindex=countindex+1;
		     binindex[countindex-1]=next;
		     rayindex[countindex-1]=currIndex;
		     if (VERBOSE) printf("nc ");
		   }
		   if (GOOD[next][currIndex]==0) {
		     countbins=countbins+1;
		   }
		   if (GOOD[next][right]==1) {
		     if (flag==0) flag=1;
		     countindex=countindex+1;
		     binindex[countindex-1]=next;
		     rayindex[countindex-1]=right;
		     if (VERBOSE) printf("nr ");
		   }
		   if (GOOD[next][right]==0) {
		     countbins=countbins+1;
		   }
		 }
         }
		 /* Perform last step of Bergen and Albers filter: */
		 if (loopcount==1 && countbins+countindex<1) 
		       GOOD[i][currIndex]=-1;  

		 if (countindex>=1) {
		   /* Unfold against all adjacent values where GOOD==1 */
		   numtimes=0;
		   while(val!=missingVal&&GOOD[i][currIndex]==0) {
		     numtimes=numtimes+1;
		     in=0;
		     out=0;
		     numneg=0;
		     numpos=0;
		     if (VERBOSE) printf("countindex: %d\n", countindex); 
		     for (l=0;l<countindex;l++) {
		       if (VERBOSE) printf("%d %d %f\n",binindex[l],
			  rayindex[l],goodval);
		       goodval=(float)rvVolume->sweep[sweepIndex]->ray
			 [rayindex[l]]->h.f(rvVolume->sweep[sweepIndex]->
			 ray[rayindex[l]]->range[binindex[l]]);
	               diffs[l]=goodval-val;
		       if (fabs(diffs[l])<pfraction*NyqVelocity) in=in+1;
		       else {
			 out=out+1;
			 if (diffs[l]>NyqVelocity) {
			   numpos=numpos+1;
			 } else if (diffs[l]<-NyqVelocity){
			   numneg=numneg+1;
			 }
		       }
		     }
		     if (in>0 && out==0) {
		       finalval=(float)rvVolume->sweep[sweepIndex]->
			 ray[currIndex]->h.invf(val);
		       rvVolume->sweep[sweepIndex]->ray[currIndex]->
			 range[i]=(unsigned short) (finalval);
		       GOOD[i][currIndex]=1;
		       if (VERBOSE) printf("Finalval: %4.2f\n", val);
		     } else {
		       if (numtimes<=MAXCOUNT) {
			 if ((numpos+numneg)<(in+out-(numpos+numneg))) {
			   if (loopcount>2)
			   { /* Keep the value after two passes through
			     ** data. */
			     finalval=(float)rvVolume->sweep[sweepIndex]->
			       ray[currIndex]->h.invf(val);
			     rvVolume->sweep[sweepIndex]->ray[currIndex]->
			       range[i]=(unsigned short) (finalval);
			     GOOD[i][currIndex]=1;
			     if (VERBOSE) printf(
				     "Keeping after two passes: %f\n",val);
			   }
			 } else if (numpos>numneg) {
			   val=val+NyqInterval;
			 } else if (numneg>numpos) {
			   val=val-NyqInterval;
			 } else {
			   /* Save bin for windowing if unsuccessful after
			   ** four passes: */
			   if (loopcount>4) GOOD[i][currIndex]=-2;
			   if (VERBOSE) printf("Saving for windowing\n");
			 }
		       } else {
			 /* Remove bin: */
			 GOOD[i][currIndex]=-2;
			 if (VERBOSE) printf("Saving for windowing\n");
		       }
		     }
		   }
		 }
	       }
	     }
	   }
	 }

	 /* Unfold remote bins or those that were previously unsuccessful
	 **   using a window with dimensions 2(PROXIMITY)+1 x 2(PROXIMITY)+1:
	 **   if still no luck delete data (or unfold against VAD if PASS2). */
	 for (i=DELNUM;i<numBins;i++) {
	   for (currIndex=0;currIndex<numRays;currIndex++) { 
	     if (GOOD[i][currIndex]==0 || GOOD[i][currIndex]==-2) {
	       val=(float) VALS->sweep[sweepIndex]->ray[currIndex]->
		   h.f(VALS->sweep[sweepIndex]->ray[currIndex]->range[i]);
	       startray=currIndex-PROXIMITY;
	       endray=currIndex+PROXIMITY;
	       firstbin=i-PROXIMITY;
	       lastbin=i+PROXIMITY;
	       if (startray<0) startray=numRays+startray;
	       if (endray>numRays-1) endray=endray-numRays;
	       if (firstbin<0) firstbin=0;
	       if (lastbin>numBins-1) lastbin=numBins-1;
	       winval=window(rvVolume, sweepIndex, startray, endray, 
			     firstbin, lastbin, &std, missingVal, &wsuccess);
	       if (winval==missingVal && wsuccess==1){ /* Expand the window: */
		 startray=currIndex-2*PROXIMITY;
		 endray=currIndex+2*PROXIMITY;
		 firstbin=i-2*PROXIMITY;
		 lastbin=i+2*PROXIMITY;
		 if (startray<0) startray=numRays+startray;
		 if (endray>numRays-1) endray=endray-numRays;
		 if (firstbin<0) firstbin=0;
		 if (lastbin>numBins-1) lastbin=numBins-1;
		 winval=window(rvVolume, sweepIndex, startray, endray, 
			       firstbin, lastbin, &std, missingVal, &wsuccess);
	       }
	       if (winval!=missingVal) {
		 diff=winval-val;
		 if (diff<0.0) {
		   diff=-diff;
		   direction=-1;
		 } else direction=1;
		 numtimes=0;
		 while (diff>0.99999*NyqVelocity && numtimes<=MAXCOUNT) {
		   val=val+NyqInterval*direction;
		   numtimes=numtimes+1;
		   diff=winval-val;
		   if (diff<0.0) {
		     diff=-diff;
		     direction=-1;
		   } else direction=1;
		 }
		 if (diff<pfraction*NyqVelocity) {
		   /* Return the value. */
		   finalval=(float)rvVolume->sweep[sweepIndex]->
		     ray[currIndex]->h.invf(val);
		   rvVolume->sweep[sweepIndex]->ray[currIndex]->
		     range[i]=(unsigned short) (finalval);
		   GOOD[i][currIndex]=1;
		 } else if (diff<(1.0 - (1.0 - pfraction)/2.0)*NyqVelocity) {
		   /* If within relaxed threshold, then return value, but
		   **   do not use to dealias other bins. */
		   finalval=(float)rvVolume->sweep[sweepIndex]->
		     ray[currIndex]->h.invf(val);
		   rvVolume->sweep[sweepIndex]->ray[currIndex]->
		     range[i]=(unsigned short) (finalval);
		   GOOD[i][currIndex]=-1;		   
		 } else {
		   /* Remove bin */
		   GOOD[i][currIndex]=-1;
		 }
	       } else { 
		 if (wsuccess==0) {
		   /* Remove bin */
		   GOOD[i][currIndex]=-1;
		 } else if (soundVolume==NULL || lastVolume==NULL) {
		   if (GOOD[i][currIndex]==0 && RM!=1) {
		     /* Leave bin untouched. */
		     val=(float) VALS->sweep[sweepIndex]->ray[currIndex]->
		       h.f(VALS->sweep[sweepIndex]->ray[currIndex]->range[i]);
		     finalval=(float)rvVolume->sweep[sweepIndex]->
		       ray[currIndex]->h.invf(val);
		     rvVolume->sweep[sweepIndex]->ray[currIndex]->
		       range[i]=(unsigned short) (finalval);
		     GOOD[i][currIndex]=-1; /* Don't use to unfold other bins*/
		   } else {
		     /* Remove bin */
		     GOOD[i][currIndex]=-1;
		   }
		 } else if (GOOD[i][currIndex]==0 && PASS2 &&
		     soundVolume!=NULL && lastVolume!=NULL) {
		   /* Leave GOOD[i][currIndex]=0 bins for a second pass.
		   ** In the second pass, we repeat unfolding, except this
		   ** time we use soundVolume for comparison instead of
		   ** lastVolume. */
		 } else {
		   /* Remove bin */
		     GOOD[i][currIndex]=-1;
		 }
	       }
	     }
	   }
	 }
	 
	 /* Beginning second pass through the data, this time using 
	 **  soundVolume only: */

	 if (lastVolume!=NULL && soundVolume!=NULL) {
	   flag=1;
	   for (currIndex=0;currIndex<numRays;currIndex++) {
	     if (sweepIndex<numSweeps-1) abIndex=findRay(rvVolume, rvVolume,
	         sweepIndex, sweepIndex+1, currIndex, missingVal);
	     for (i=DELNUM;i<numBins;i++) {
	       if (GOOD[i][currIndex]==0) {
		 val=(float) VALS->sweep[sweepIndex]->ray[currIndex]->
		   h.f(VALS->sweep[sweepIndex]->ray[currIndex]->range
		       [i]);
		 valcheck=val;
		 soundval=(float) soundVolume->
		   sweep[sweepIndex]->ray[currIndex]->h.f(soundVolume->
		    sweep[sweepIndex]->ray[currIndex]->range[i]);
	         
		 if (soundval!=missingVal&&val!=missingVal) {
		   diff=soundval-val;	     
		   if (diff<0.0) {
		     diff=-diff;
		     direction=-1;
		   } else direction=1;
		   numtimes=0;
		   while (diff>0.99999*NyqVelocity&&numtimes<=MAXCOUNT) {
		     val=val+NyqInterval*direction;
		     numtimes=numtimes+1;
		     diff=soundval-val;
		     if (diff<0.0) {
		       diff=-diff;
		       direction=-1;
		     } else direction=1;
		   }
		   if (diff<fraction2*NyqVelocity&&fabs(valcheck)>CKVAL) {
		     /* Return the good value. */
		     finalval=(float)rvVolume->sweep[sweepIndex]->
		       ray[currIndex]->h.invf(val);
		     rvVolume->sweep[sweepIndex]->ray[currIndex]->
		       range[i]=(unsigned short) (finalval);
		     GOOD[i][currIndex]=1;
		   }
		 }
	       }
	     }
	   }
		   
	   /* Now, try to unfold the rest of the GOOD=0 bins, assuming spatial
	      continuity: */
	   loopcount=0;
	   while (flag==1) {
	     loopcount=loopcount+1;
	     flag=0;
	     if (step==1) {
	       step=-1;
	       startindex=numRays-1;
	       endindex=-1;
	     } else {
	       step=1;
	       startindex=0;
	       endindex=numRays;
	     }
	     for (currIndex=startindex;currIndex!=endindex;currIndex=
		      currIndex+step) {
	       for (i=DELNUM;i<numBins;i++) {
		 if (GOOD[i][currIndex]==0) {
		   countindex=0;
		   val=(float) VALS->sweep[sweepIndex]->ray[currIndex]
		     ->h.f(VALS->sweep[sweepIndex]->ray[currIndex]->
			   range[i]);
		   valcheck=val;
		   if (currIndex==0) left=numRays-1;
		   else left=currIndex-1;
		   if (currIndex==numRays-1) right=0;
		   else right=currIndex+1;
		   next=i+1;
		   prev=i-1;
		   /* Look at all bins adjacent to current bin in
		      question: */
		   if (i>DELNUM) {
		     if (GOOD[prev][left]==1) {
		       if (flag==0) flag=1;
		       countindex=countindex+1;
		       binindex[countindex-1]=prev;
		       rayindex[countindex-1]=left;
		       if (VERBOSE) printf("pl ");
		     }
		     if (GOOD[prev][currIndex]==1) {
		       if (flag==0) flag=1;
		       countindex=countindex+1;
		       binindex[countindex-1]=prev;
		       rayindex[countindex-1]=currIndex;
		       if (VERBOSE) printf("pc ");
		     }
		     if (GOOD[prev][right]==1) {
		       if (flag==0) flag=1;
		       countindex=countindex+1;
		       binindex[countindex-1]=prev;
		       rayindex[countindex-1]=right;
		       if (VERBOSE) printf("pr ");
		     }
		   }
		   if (GOOD[i][left]==1) {
		     if (flag==0) flag=1;
		     countindex=countindex+1;
		     binindex[countindex-1]=i;
		     rayindex[countindex-1]=left;
		     if (VERBOSE) printf("il ");
		   }
		   if (GOOD[i][right]==1) {
		     if (flag==0) flag=1;
		     countindex=countindex+1;
		     binindex[countindex-1]=i;
		     rayindex[countindex-1]=right;
		     if (VERBOSE) printf("ir "); 	   
		   }
		   if (i<numBins-1) {  
		     if (GOOD[next][left]==1) {
		       if (flag==0) flag=1;
		       countindex=countindex+1;
		       binindex[countindex-1]=next;
		       rayindex[countindex-1]=left;
		       if (VERBOSE) printf("nl "); 
		     }
		     if (GOOD[next][currIndex]==1) {
		       if (flag==0) flag=1;
		       countindex=countindex+1;
		       binindex[countindex-1]=next;
		       rayindex[countindex-1]=currIndex;
		       if (VERBOSE) printf("nc ");
		     }
		     if (GOOD[next][right]==1) {
		       if (flag==0) flag=1;
		       countindex=countindex+1;
		       binindex[countindex-1]=next;
		       rayindex[countindex-1]=right;
		       if (VERBOSE) printf("nr ");
		     }
		   }
		   /* Unfold against all adjacent values with GOOD==1*/
		   if (countindex>=1) {
		     numtimes=0;
		     while(val!=missingVal&&GOOD[i][currIndex]==0) {
		       numtimes=numtimes+1;
		       in=0;
		       out=0;
		       numneg=0;
		       numpos=0;
		       if (VERBOSE) printf("%d: ", countindex); 
		       for (l=0;l<countindex;l++) {
			 goodval=(float)rvVolume->sweep[sweepIndex]->
			   ray[rayindex[l]]->h.f(rvVolume->sweep
			    [sweepIndex]->ray[rayindex[l]]->range
			    [binindex[l]]);
			 diffs[l]=goodval-val;
			 if (fabs(diffs[l])<pfraction*NyqVelocity) in=in+1;
			 else {
			   out=out+1;
			   if (diffs[l]>NyqVelocity) {
			     numpos=numpos+1;
			   } else if (diffs[l]<-NyqVelocity){
			     numneg=numneg+1;
			   }
			 }
		       }
		       if (in>out) {
			 finalval=(float)rvVolume->sweep[sweepIndex]->
			   ray[currIndex]->h.invf(val);
			 rvVolume->sweep[sweepIndex]->ray[currIndex]->
			   range[i]=(unsigned short) (finalval);
			 GOOD[i][currIndex]=1;
			 if (VERBOSE) printf("Val: %4.2f\n", val);
		       } else {
			 if (numtimes<=MAXCOUNT) {
			   if (numpos+numneg<(in+out-(numpos+numneg))) {
			     if (loopcount<=2) val=missingVal; /* Try later */
			     else {
			       /* Keep the value after two passes through
			       ** data */
			       finalval=(float)rvVolume->sweep[sweepIndex]->
				 ray[currIndex]->h.invf(val);
			       rvVolume->sweep[sweepIndex]->ray[currIndex]->
				 range[i]=(unsigned short) (finalval);
			       GOOD[i][currIndex]=1;
			     }
			   } else if (numpos>numneg) {
			     val=val+NyqInterval;
			   } else if (numneg>numpos) {
			     val=val-NyqInterval;
			   } else {
			     /* Remove bin after four passes through data: */
			     if (loopcount>4) GOOD[i][currIndex]=-1;
			   }
			 } else {
			   /* Remove bin: */
			   GOOD[i][currIndex]=-1;
			 }
		       }
		     }
		   }
		 }
	       }
	     }
	   }
	 }
       }
       *success=1;
     } else {
       printf("First guess not available.\n");
       *success=0;
     }
     return;
}
