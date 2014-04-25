
#include <rsl.h>

/* Gate value getters */

float ray_val(Ray *ray, int index) 
{
    /* Return the value of range bin at index in ray */
    return (float) ray->h.f(ray->range[index]);
}

void ray_set(Ray *ray, int index, float val) 
{
    /* Set the gate value in ray at index to val */
    ray->range[index]=(unsigned short) ray->h.invf(val);
    return;
}
