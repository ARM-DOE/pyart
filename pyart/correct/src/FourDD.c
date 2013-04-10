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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <rsl.h>
#include "FourDD.h"
#define MISSINGVEL 131072.0

main(int argc, char** argv) 
{
    char lastFN[256];
    char inputFN[256];
    char outputFN[256];
    char sounding_name[256];
    Radar* lastVolume;
    Radar* radarData;
    Radar* radarOut;
    Volume* DZvolume;
    Volume* lastVelVolume;
    Volume* radialVelVolume;
    Volume* unfoldedVolume;
    Volume* soundVolume;
    FILE* sond_file;

    unsigned short prep = 1, outDZ, outVR, filt = 0, nolast=0, success=0,
      usuccess=0;
    int VAD_time;
    int numSweepsLast;
    int numSweepsCurrent;

    if (argc < 10) {
	printf("Not enough arguments\n");
	printf("usage: 4DD lastfile inputfile outputfile VADfile VADtime prep filt outputDZ outputVR\n");
	printf("  lastfile = filename of previously dealiased volume\n");
	printf("  inputfile = filename of volume you want to dealias\n");
	printf("  outputfile = filename of output volume\n");
	printf("  VADfile = name of output file\n");
	printf("  VADtime = time of VAD in (YYDDDHHMM)\n");
	printf("  prep = thresholding flag (1=yes, 0=no)\n");
	printf("  filt = Bergen and Albers filter (1=yes, 0=no)\n");
	printf("  outputDZ = output flag for DZ field (1=yes, 0=no)\n");
	printf("  outputVR = output flag for VR field (1=yes, 0=no)\n");
	exit(-1);
    } 

    strcpy(lastFN,argv[1]);
    printf("Last file = %s \n", lastFN);
    strcpy(inputFN,argv[2]);
    printf("Input file = %s \n", inputFN); 
    strcpy(outputFN,argv[3]);
    printf("output file = %s \n", outputFN); 
    strcpy(sounding_name,argv[4]);
    printf("VAD file = %s", sounding_name);
    VAD_time=atoi(argv[5]);
    prep=atoi(argv[6]);
    printf(", time = %d\n", VAD_time);
    filt=atoi(argv[7]);
    outDZ=atoi(argv[8]);
    outVR=atoi(argv[9]);
 
    lastVolume = RSL_uf_to_radar(lastFN);
    radarData = RSL_uf_to_radar(inputFN);
    /*lastVolume = RSL_lassen_to_radar(lastFN);
    radarData = RSL_lassen_to_radar(inputFN);*/


    if (radarData == NULL) {
	printf("Conversion of uf file %s to radar format failed; aborting\n",
	       inputFN);
	exit(-1);
    }

    /* Create a new radar structure for the output */
    radarOut             = RSL_new_radar(MAX_RADAR_VOLUMES); 
    radarOut->h          = radarData->h;
    radarOut->h.nvolumes = MAX_RADAR_VOLUMES; /* hardwired */

    if (prep || outDZ) DZvolume = radarData->v[CZ_INDEX];

    /* error handling if there is no DZ field */
    if (DZvolume == NULL && outDZ){
	printf("No DZ field available\n");
	outDZ=0;
	prep=0;
    }
   
    if (radarData->v[VR_INDEX] != NULL)  {
        radialVelVolume = radarData->v[VR_INDEX];
    } else {
        printf("No VR field available to unfold; aborting\n");
        exit(-1);
    }
    if (lastVolume !=NULL) {
      lastVelVolume=lastVolume->v[VE_INDEX]; /* Loading edited radial velocity
					     ** field from previous volume. */
      numSweepsLast=lastVelVolume->h.nsweeps;
      numSweepsCurrent=radialVelVolume->h.nsweeps;
    } else {
      lastVelVolume=NULL;
      nolast=1;
    }
    printf("Trying to read sounding");
    /* Check to see if VAD or sounding available. If so, create first guess
    **    field: */
    if ((sond_file = fopen(sounding_name, "r"))==NULL) {
       soundVolume = NULL;
       printf("%s not readable -- unfolding using radial continuity\n"
	  ,sounding_name);
    } else {
      soundVolume = RSL_copy_volume(radialVelVolume);
      /* Here, if lastVelVolume is NULL, we compute a first guess VAD field
      ** in the array soundVolume. */
      firstGuess(soundVolume, MISSINGVEL, sounding_name, VAD_time, &success);
    }
    close(sond_file);
    if (success==0) {
      soundVolume=NULL;
    }

    if (success==1 || lastVelVolume!=NULL) { /* Proceed with unfolding */

      printf("Dealiasing: \n");
      /* Copy radial velocity data to unfoldedVolume and remove any
      **   bins where reflectivity is missing (if NODBZRMRV==1) or
      **   ouside the interval LOWDBZ<= dz <=HIGHDBZ. */      
     printf("LOWDBZ= %f\n", LOWDBZ);
     printf("HIGHDBZ= %f\n", HIGHDBZ);
     unfoldedVolume = RSL_copy_volume(radialVelVolume);
      if (prep) prepVolume(DZvolume, unfoldedVolume, MISSINGVEL);
      
      /* Finally, we call unfoldVolume, which unfolds the velocity
      ** field in the array unfoldedVolume. The subroutine returns a value
      ** of 1 if unfolding is performed. */
      unfoldVolume(unfoldedVolume,
		   soundVolume,lastVelVolume,MISSINGVEL, filt, &usuccess);
    }

    if (unfoldedVolume==NULL || usuccess==0) {
      printf("Velocity Volume was NOT unfolded!\n");
      unfoldedVolume=RSL_copy_volume(radialVelVolume);
    }
    else {
      printf("Velocity Volume was unfolded.\n");
    }      
   
    /* now finish making the new radar field: */
    
    
    if (outDZ) radarOut->v[CZ_INDEX] = DZvolume;
    radarOut->v[VE_INDEX] = unfoldedVolume;
    if (outVR) radarOut->v[VR_INDEX] = radialVelVolume;
    /*
    ** Uncomment the following line if you want to output soundVolume, 
    ** which is an artificial volume created by the VAD or sounding. */
    if (soundVolume!=NULL) radarOut->v[TI_INDEX] = soundVolume; 
    printf("Writing out radar ... \n");
    RSL_radar_to_uf(radarOut,outputFN);
    printf("freeing memory ... \n");
    RSL_free_radar(radarOut);
    return 0;
} /* main */

