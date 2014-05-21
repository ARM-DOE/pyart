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

#include <rsl.h> /* Sweep */ 

float ray_val(Ray *ray, int index);
void ray_set(Ray *ray, int index, float val); 

#endif /* DEALIAS_H */
