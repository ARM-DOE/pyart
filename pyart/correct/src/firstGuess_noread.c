/*
**
**  UW Radial Velocity Dealiasing Algorithm
**  Four-Dimensional Dealiasing (4DD)
**
**  DESCRIPTION:
**      This routine creates a firstGuess radial velocity field given a
**      sounding or VAD. Assumes standard atmosphere refraction (4Rearth/3) 
**      and extrapolates sounding data to all radar bins (assuming the wind
**      is horizontally uniform).
**
**  DEVELOPER:
**	Curtis N. James     08 Dec 98
**
**
**  HISTORY:
**	An elaboration of the NSSL-Eilts algorithm.
**
**
**
*/
#include "FourDD.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <rsl.h> /* Sweep */


void firstGuessNoRead(Volume* soundVolume, float missingVal,
          float *height_array, float *speed_array, float *direction_array, int nlevels, int VAD_time, unsigned short* sounding) {

     int numLevs, alt, i,idx, sweepIndex, currIndex, index, numBins, numRays, 
       numSweeps, time, naz;
     unsigned short sond_type = 0; /* sond_type: 0=sounding and 1=VAD */
     unsigned short flag = 0;
     float ua_data[1000][5],data;
     float ke,dRdz,height,rnge,elev,az,start_range,h_range,gate_size,val
       ,wind, wind_val_rv,dir,offset,ang, U, V, H1, H2, meanShearU=0.0,
       meanShearV=0.0;   
     /*float sumz, sumz2, sumz3, sumz4, N, sumU, sumUz, sumUz2, sumV, sumVz,
     **  sumVz2, DET, DETAU0, DETAU1, DETAU2, DETAV0, DETAV1, DETAV2;
     **     Summations and determinants for least-squares quadratic fit
     ** float AU0, AU1, AU2, AV0, AV1, AV2; Final coefficients of quadratic
     **					    for U and V least-squares fit */
     char line[100];
     FILE *sond_file;

	 /* Sounding */
	 i=1;
	 ua_data[0][0]=0.0;
	 ua_data[0][1]=0.0;
	 ua_data[0][2]=0.0;
	 ua_data[0][3]=0.0;
	 ua_data[0][4]=0.0;

	 for (idx = 1; idx < nlevels; ++idx) {
			ua_data[i][0]=height_array[idx]; /* Height (m) */
			ua_data[i][1]=sin(direction_array[idx]*PI/180.0)*speed_array[idx]; /* U REMOVED - (negative)*/
			ua_data[i][2]=cos(direction_array[idx]*PI/180.0)*speed_array[idx]; /* V REMOVED - (negative)*/
			ua_data[i][3]=(ua_data[i][1]-ua_data[i-1][1])/
			     (ua_data[i][0]-ua_data[i-1][0]); /* Shear U */
			   ua_data[i][4]=(ua_data[i][2]-ua_data[i-1][2])/
			     (ua_data[i][0]-ua_data[i-1][0]); /* Shear V */
	                   if (fabs(ua_data[i][3])<=MAXSHEAR&&fabs(ua_data[i][4])<=
					     MAXSHEAR) {
	                	   printf("height %f u %f v %f dir %f speed %f\n", ua_data[i][0], ua_data[i][1], ua_data[i][2], direction_array[i], speed_array[i]);
			     i++;
			     meanShearU=meanShearU+ua_data[i][3];
			     meanShearV=meanShearV+ua_data[i][4];
	                   	   	   	   }
	 }
	 numLevs=i-1;
	 /* Force wind at ground to be same as that of first level: */
	 ua_data[0][1]=ua_data[1][1];
	 ua_data[0][2]=ua_data[1][2];
	 ua_data[1][3]=0.0;
	 ua_data[1][4]=0.0;
	 /* Calculate mean shear: */
	 if (numLevs>1) {
	   meanShearU=meanShearU/(ua_data[numLevs-1][0]-ua_data[1][0]);
	   meanShearV=meanShearV/(ua_data[numLevs-1][0]-ua_data[1][0]);
	 } else {
	   meanShearU=0.0;
	   meanShearV=0.0;
	 }
	 printf("Number of sounding levels used: %d\n\n", numLevs);

       numSweeps = soundVolume->h.nsweeps;
       dRdz=-39.2464; /* Standard Atmosphere refractivity gradient in km^-1 */
       alt=soundVolume->sweep[0]->ray[0]->h.alt;
/*       fprintf(stderr,"Radar altitude: %d\n",alt);*/

       for(sweepIndex = 0; sweepIndex < numSweeps; sweepIndex++) {
	
	 numRays=soundVolume->sweep[sweepIndex]->h.nrays;
	 numBins=soundVolume->sweep[sweepIndex]->ray[0]->h.nbins;
	 start_range=soundVolume->sweep[sweepIndex]->ray[0]->
	    h.range_bin1;
	 gate_size=soundVolume->sweep[sweepIndex]->ray[0]->
	    h.gate_size;
	 elev=PI*(soundVolume->sweep[sweepIndex]->ray[0]->h.elev)/180.0;
	 index=0;
	 /* Now we create a first-guess velocity field using the sounding
	 ** or VAD profile. (We assume standard atmospheric refraction).*/
	 for(i = 0; i < numBins; i++) {
	   /* To print out a range circle of radial velocity values: */
	   rnge=start_range+i*gate_size+gate_size/2.0;
	   ke=1/(1+A*dRdz*pow(10,-6)); /* Doviak and Zrnic, 1993 */
	   height=sqrt(pow(rnge,2)+pow(ke*A*1000,2)+2*rnge*ke*A*1000
              *sin(elev))-ke*A*1000+alt;
	   h_range=ke*A*1000*asin(rnge*cos(elev)/(ke*A*1000+height-alt));
	   ang=atan(cos(elev)*sin(elev+h_range/ke/A/1000)*pow(cos(elev+      		    h_range/ke/A/1000),-2)); /* atan (dh/ds) */
	 /*fprintf(stderr,"Sweep= %d, Elev= %2f, ", sweepIndex, elev*180/PI);
           fprintf(stderr,"Height= %0.2f, Angle= %0.2f, ", height, ang/PI*
              180.0);
           fprintf(stderr,"ke= %0.2f, h_range= %0.2f, Range= %0.2f\n",
	      ke, h_range, range);*/
	   if (height>=ua_data[0][0] && height<=ua_data[numLevs-1][0]) {
	     while(1) {
 	       if (height>=ua_data[index][0] && height<ua_data[index+1][0]) {
		 U=ua_data[index][1]+ua_data[index+1][3]*(height-ua_data
		    [index][0]);
		 V=ua_data[index][2]+ua_data[index+1][4]*(height-ua_data
		    [index][0]);
		 wind=sqrt(pow(U,2)+pow(V,2));
		 if (SIGN<0) offset=0.0;
		 else offset=PI;
		 if (U>=0) dir=(acos(V/wind)+offset)*180/PI;
		 else dir=(offset-acos(V/wind))*180/PI;
		 break;
	       } 
	       else if (height<ua_data[index][0]) index--;
	       else index++;
	     }
	   } else if (height>ua_data[numLevs-1][0]) {
	     U=ua_data[numLevs-1][1]+meanShearU*(height-ua_data
						      [numLevs-1][0]);
	     V=ua_data[numLevs-1][2]+meanShearV*(height-ua_data
						      [numLevs-1][0]);
	     wind=sqrt(pow(U,2)+pow(V,2));
	     if (SIGN<0) offset=0.0;
	     else offset=PI;
	     if (U>=0) dir=(acos(V/wind)+offset)*180/PI;
	     else dir=(offset-acos(V/wind))*180/PI;
           } 
	   for(currIndex = 0; currIndex < numRays; currIndex++) {
	     val=(float) soundVolume->sweep[sweepIndex]->ray[currIndex]
                ->h.f(soundVolume->sweep[sweepIndex]->ray[currIndex]->
                range[i]);
	     /* if (val!=missingVal) { */
	       if (wind>=0.0 && dir>=0.0) {
		 az=PI*(soundVolume->sweep[sweepIndex]->ray[currIndex]->
		   h.azimuth)/180.0;
		 wind_val_rv=wind*cos(ang);
		 wind_val_rv=wind_val_rv*cos(PI*dir/180.0-az);
                 val=(float)soundVolume->sweep[sweepIndex]->ray[currIndex]
                   ->h.invf(wind_val_rv);
                 soundVolume->sweep[sweepIndex]->ray[currIndex]->range[i]=
                   (unsigned short) (val); 
		 /* wind_val_az=wind*cos(ang);
		 ** wind_val_az=-wind_val_az*sin(PI*dir/180.0-az);
                 ** val=(float)soundVolumeaz->sweep[sweepIndex]->ray[currIndex]
                 **  ->h.invf(wind_val_az);
                 ** soundVolumeaz->sweep[sweepIndex]->ray[currIndex]->range[i]=
                 **  (unsigned short) (val); */
		 if (flag==0) flag=1;
	       } else {
		 wind_val_rv=missingVal;
                 val=(float)soundVolume->sweep[sweepIndex]->ray[currIndex]
                   ->h.invf(wind_val_rv);
                 soundVolume->sweep[sweepIndex]->ray[currIndex]->range[i]=
                   (unsigned short) (val);
	       } 	 	       
	       /* } */ 
	   }
	 }
       }
     if (numLevs==0) {
       flag=0;
     }
     if (flag) *sounding=1;
     return;
}









