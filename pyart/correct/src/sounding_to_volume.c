/*
 * Create a radial velocity RSL Volume from sounding data.
 *  
 * Adapted from routines from:
 * 
 * UW Radial Velocity Dealiasing Algorithm
 * Four-Dimensional Dealiasing (4DD)
 * 
 * Developer by Curtis N. James     08 Dec 98
 *
 * See UWASHINGTON_4DD_README file for license
 *
 * Adapted for use in Py-ART by Jonathan J. Helmus May 2014
 */


#include "helpers.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <rsl.h>

#define A 6372.0 /* Radius of earth in kilometers. */
#define PI 3.1415927 

/******************************
 * Private data and functions *
 ******************************/

typedef struct Sounding {
    float data[1000][5];    /* sounding data; height, U, V, shearU, shearV */
    int numLevs;            /* number of levels recorded in data */
    float offset;           /* offset to use when computing direction */
    float meanShearU;       /* mean U shear */
    float meanShearV;       /* mean V shear */
} Sounding;

/* Load sounding data from arrays into a Sounding structure.
 * Filter out data with shear greater than maxshear. */
static void load_sounding(
    float *height_array, float *speed_array, float *direction_array, 
    int nlevels, float maxshear, Sounding *sound)
{
    int i, idx;
    i=1;
    sound->data[0][0]=0.0;
    sound->data[0][1]=0.0;
    sound->data[0][2]=0.0;
    sound->data[0][3]=0.0;
    sound->data[0][4]=0.0;
    sound->meanShearU = 0;
    sound->meanShearV = 0;

    /* store sounding data which meets shear conditions */
    for (idx = 1; idx < nlevels; ++idx) {
        /* Height (m) */
        sound->data[i][0] = height_array[idx]; 
        /* U REMOVED - (negative)*/
        sound->data[i][1] = sin(direction_array[idx] * PI / 180.0) * 
                        speed_array[idx];
        /* V REMOVED - (negative)*/
        sound->data[i][2] = cos(direction_array[idx] * PI / 180.0) *
                        speed_array[idx];
        /* Shear U */
        sound->data[i][3] = (sound->data[i][1] - sound->data[i-1][1]) /
                        (sound->data[i][0] - sound->data[i-1][0]);
        /* Shear V */
        sound->data[i][4] = (sound->data[i][2] - sound->data[i-1][2]) /
                        (sound->data[i][0] - sound->data[i-1][0]);
        if (fabs(sound->data[i][3]) <= maxshear && 
            fabs(sound->data[i][4]) <= maxshear) {
            i++;
            sound->meanShearU = sound->meanShearU + sound->data[i][3];
            sound->meanShearV = sound->meanShearV + sound->data[i][4];
        }
    } 
    /* Force wind at ground to be same as that of first level: */
    sound->data[0][1] = sound->data[1][1];
    sound->data[0][2] = sound->data[1][2];
    sound->data[1][3] = 0.0;
    sound->data[1][4] = 0.0;
    
    /* Calculate mean shear: */
    if (i>2) {
        sound->meanShearU = sound->meanShearU / 
                            (sound->data[i][0] - sound->data[1][0]);
        sound->meanShearV = sound->meanShearV / 
                            (sound->data[i][0] - sound->data[1][0]);
    } else {
        sound->meanShearU=0.0;
        sound->meanShearV=0.0;
    }
    sound->numLevs = i-1;
    return;
}

/* Calculate direction and speed from wind parameters */
static void calc_wind_dir(
    Sounding *sound, int index, float height, float *wind, float *dir)
{
    float U, V, shearU, shearV;
    if (index >= sound->numLevs-1) {
        shearU = sound->meanShearU;
        shearV = sound->meanShearV;
    } else {
        shearU = sound->data[index+1][3];
        shearV = sound->data[index+1][4];
    }
    U = sound->data[index][1] + shearU * (height - sound->data[index][0]);
    V = sound->data[index][2] + shearV * (height - sound->data[index][0]);
    *wind = sqrt(pow(U,2) + pow(V,2));
    if (U>=0)
        *dir = (acos(V / *wind) + sound->offset) * 180 / PI;
    else 
        *dir = (sound->offset - acos(V / *wind)) * 180 / PI;
    return;
}

/* Calculate the wind speed and direction by interpolating from the 
 * nearest measurements in the sounding data. */
static void interpolate_wind(
    Sounding *sound, float height, float *wind, float *dir)
{
    static int index=0; 
    if (height >= sound->data[0][0] && 
        height <= sound->data[sound->numLevs-1][0]) {
        /* height is between two soundings measurments, 
         * find and interpolate from lower measurement */
        while(1) {
            if (height >= sound->data[index][0] &&
                height < sound->data[index+1][0]) {
                calc_wind_dir(sound, index, height, wind, dir);
                break;
            } else if (height < sound->data[index][0]) {
                index--;
            } else {
                index++;
            }
        }    
    } else if (height > sound->data[sound->numLevs-1][0]) {
        /* Height is above the heightest sounding measurement, extrapolate
         * from the heighest measurments using mean shear. */
        calc_wind_dir(sound, sound->numLevs-1, height, wind, dir);
    }
}

/********************
 * Public functions *
*********************/

/*
 * Creates a radial velocity RSL volume from a sounding. 
 *
 * Assumes standard atmosphere refraction (4 * R_earth / 3) and 
 * interpolates sounding data to all radar bins (assuming the wind
 * is horizontally uniform). An elaboration of the NSSL-Eilts algorithm.
 *
 * maxshear : Maximum vertical shear allowed in input sounding 
 * sign : Sign convention, if sign=-1, negative radial velocity is
 *        towards the radar, if sign=1 positive value towards radar.
 */
int sounding_to_volume(
    Volume* soundVolume, float missingVal,
    float *height_array, float *speed_array, float *direction_array, 
    int nlevels, float maxshear, int sign) 
{
    int alt, i, sweepIndex, currIndex, flag = 0;
    float ke, dRdz, height, rnge, elev, az, start_range, h_range, gate_size,
          wind, wind_val_rv, dir, ang;
    Sweep *sweep;
    Sounding sound;

    load_sounding(height_array, speed_array, direction_array, nlevels, 
                  maxshear, &sound);
    if (sign < 0) 
        sound.offset=0.0;
    else 
        sound.offset=PI;
 
    dRdz=-39.2464; /* Standard Atmosphere refractivity gradient in km^-1 */
    ke = 1/(1 + A * dRdz * pow(10,-6)); /* Doviak and Zrnic, 1993 */
    alt = soundVolume->sweep[0]->ray[0]->h.alt;
    wind = dir = 0.0;
    
    /* Create a first-guess velocity field using the sounding
     * data , assuming standard atmospheric refraction. */
    for(sweepIndex = 0; sweepIndex < soundVolume->h.nsweeps; sweepIndex++) {
        
        sweep = soundVolume->sweep[sweepIndex]; 
        start_range = sweep->ray[0]->h.range_bin1;
        gate_size = sweep->ray[0]->h.gate_size;
        elev = PI * (sweep->ray[0]->h.elev) / 180.0;
        
        for(i = 0; i < (sweep->ray[0]->h.nbins); i++) {
            
            /* To print out a range circle of radial velocity values: */
            rnge = start_range + i * gate_size + gate_size / 2.0;
            height = sqrt(pow(rnge,2) + pow(ke*A*1000, 2) + 
                          2*rnge*ke*A*1000*sin(elev)) - ke*A*1000 + alt;
            h_range = ke*A*1000*asin(rnge*cos(elev) / 
                      (ke*A*1000 + height - alt));
            ang = atan(cos(elev) * sin(elev+h_range/ke/A/1000) * 
                pow(cos(elev + h_range/ke/A/1000), -2) ); /* atan (dh/ds) */
            
            for(currIndex = 0; currIndex < sweep->h.nrays; currIndex++) {
                
                interpolate_wind(&sound, height, &wind, &dir);
                if (wind>=0.0 && dir>=0.0) {
                    az = PI * (sweep->ray[currIndex]->h.azimuth) / 180.0;
                    wind_val_rv = wind * cos(ang) * cos(PI*dir/180.0 - az); 
                    ray_set(sweep->ray[currIndex], i, wind_val_rv);
                    flag=1;
                } else {
                    ray_set(sweep->ray[currIndex], i, missingVal);
                }               
            }
        }
    }
    if (sound.numLevs==0) {
        flag=0;
    }
    return flag;
}





