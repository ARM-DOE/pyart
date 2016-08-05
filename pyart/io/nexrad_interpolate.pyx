"""
pyart.io.nexrad_interpolate
===========================

Interpolation of NEXRAD moments from 1000 meter to 250 meter gate spacing.

.. autosummary::
    :toctree: generated/

    _fast_interpolate_scan

"""

def _fast_interpolate_scan(
        float[:, :] data, float[:] scratch_ray, float fill_value, 
        int start, int end, int moment_ngates, int linear_interp):
    """ Interpolate a single NEXRAD moment scan from 1000 m to 250 m. """
    # This interpolation scheme is only valid for NEXRAD data where a 4:1
    # (1000 m : 250 m) interpolation is needed.
    #
    # The scheme here performs a linear interpolation between pairs of gates
    # in a ray when the both of the gates are not masked (below threshold).
    # When one of the gates is masked the interpolation changes to a nearest
    # neighbor interpolation. Nearest neighbor is also performed at the end
    # points until the new range bin would be centered beyond half of the range
    # spacing of the original range.
    #
    # Nearest neighbor interpolation is performed when linear_interp is False,
    # this is equivalent to repeating each gate four times in each ray.
    #
    # No transformation of the raw data is performed prior to interpolation, so
    # reflectivity will be interpolated in dB units, velocity in m/s, etc,
    # this may not be the best method for interpolation.
    #
    # This method was adapted from Radx
    cdef int ray_num, i, interp_ngates
    cdef float gate_val, next_val, delta

    interp_ngates = 4 * moment_ngates  # number of gates interpolated

    for ray_num in range(start, end+1):

        # repeat each gate value 4 times
        for i in range(moment_ngates):
            gate_val = data[ray_num, i]
            scratch_ray[i*4 + 0] = gate_val
            scratch_ray[i*4 + 1] = gate_val
            scratch_ray[i*4 + 2] = gate_val
            scratch_ray[i*4 + 3] = gate_val

        if linear_interp:
            # linear interpolate
            for i in range(2, interp_ngates - 4, 4):
                gate_val = scratch_ray[i]
                next_val = scratch_ray[i+4]
                if gate_val == fill_value or next_val == fill_value:
                    continue
                delta = (next_val - gate_val) / 4.
                scratch_ray[i+0] = gate_val + delta * 0.5
                scratch_ray[i+1] = gate_val + delta * 1.5
                scratch_ray[i+2] = gate_val + delta * 2.5
                scratch_ray[i+3] = gate_val + delta * 3.5

        for i in range(interp_ngates):
            data[ray_num, i] = scratch_ray[i]
