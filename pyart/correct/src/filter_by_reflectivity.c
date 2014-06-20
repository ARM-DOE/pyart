/*
 * Filter a velocity volume based on corresponding reflectivities.
 *  
 * Adapted from routines from:
 * 
 * UW Radial Velocity Dealiasing Algorithm
 * Four-Dimensional Dealiasing (4DD)
 * 
 * Developer by Curtis N. James     20 Mar 98, 25 Jan 99
 *
 * See UWASHINGTON_4DD_README file for license
 *
 * Adapted for use in Py-ART by Jonathan J. Helmus May 2014
 */

#include "helpers.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <rsl.h> /* Sweep */

/* Test if sweeps are compatible, 0 for no, 1 for yet */
static int sweeps_compatible(Sweep *sweep1, Sweep *sweep2)
{
    /* number of rays in the sweep must be the same */
    if (sweep1->h.nrays != sweep2->h.nrays)
        return -1;
    /* gates must begin at same location */
    if (sweep1->ray[0]->h.range_bin1 != sweep2->ray[0]->h.range_bin1)
        return -1;
    /* neither gate size can be zero */
    if (sweep1->ray[0]->h.gate_size == 0)
        return -1;
    if (sweep2->ray[0]->h.gate_size == 0)
        return -1;
    return 1;
}

/*
 *  Mark velocity bins missing where the corresponding reflectivity
 *  is below lowdbz, above highdbz or missing (if rm_missing==1).
*/
int filter_by_reflectivity(
    Volume* ref_volume, Volume* vel_volume, float missingVal,
    float lowdbz, float highdbz, int rm_missing)
{
    int ray_index, sweep_index, i, j, ref_index, dbz_per_vel;
    float ref;
    div_t gate_ratio; 
    Sweep *vel_sweep, *dbz_sweep;

    if (ref_volume == NULL || vel_volume == NULL)
        return 0;
    if (vel_volume->h.nsweeps != ref_volume->h.nsweeps)
        return 0;

    for (sweep_index=0; sweep_index < (vel_volume->h.nsweeps); sweep_index++) {
        vel_sweep = vel_volume->sweep[sweep_index];
        dbz_sweep = ref_volume->sweep[sweep_index];
        if (sweeps_compatible(vel_sweep, dbz_sweep) < 0)
            return 0;
    
        /* Determine the number of reflectivity gates per velocity gate */
        gate_ratio = div(vel_sweep->ray[0]->h.gate_size,
                         dbz_sweep->ray[0]->h.gate_size);
        if (gate_ratio.rem != 0)
            return 0;     /* Gate sizes must be integer multiples */
        else
            dbz_per_vel = gate_ratio.quot;  

        /* Assign missingVal to any velocity bins where the reflectivity
         * is missing (if rm_missing==1) or outside lowdbz and highdbz. */
        for (ray_index=0; ray_index < (vel_sweep->h.nrays); ray_index++) {
            for (i=0; i < (vel_sweep->ray[0]->h.nbins); i++) {
                for (j=0; j < dbz_per_vel; j++) {
                    ref_index =i * dbz_per_vel + j;
                    if (ref_index >= (dbz_sweep->ray[0]->h.nbins))
                        return 0;     /* Index out of bounds */

                    ref = ray_val(dbz_sweep->ray[ray_index], ref_index);
                    
                    if (ref < lowdbz)
                        ray_set(vel_sweep->ray[ray_index], i, missingVal);
                    if (ref > highdbz)    
                        ray_set(vel_sweep->ray[ray_index], i, missingVal);
                    if ((rm_missing==1) && fabs(ref - missingVal) < 0.00001)
                        ray_set(vel_sweep->ray[ray_index], i, missingVal);
                }
            }
        }
    }
    return 1;
}
