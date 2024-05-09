"""
Cython classes and functions for efficient mapping of radar gates to
a uniform grid.

"""

cimport cython
from cython.view cimport array as cvarray
from libc.math cimport asin, ceil, cos, exp, floor, sin, sqrt, tan

import numpy as np


# constants
cdef int BARNES2 = 3
cdef int NEAREST = 2
cdef int CRESSMAN = 1
cdef int BARNES = 0
cdef float PI = 3.141592653589793
cdef float R = 8494666.66666667   # 4/3 earths radius of 6371 km in meters

# This definition can be added to a .pxd file so others can defined fast
# RoI functions
cdef class RoIFunction:
    """ A class for storing radius of interest calculations. """

    cpdef float get_roi(self, float z, float y, float x):
        """ Return the radius of influence for coordinates in meters. """
        return 0


cdef class ConstantRoI(RoIFunction):
    """ Constant radius of influence class. """

    cdef float constant_roi

    def __init__(self, float constant_roi):
        """ intialize. """
        self.constant_roi = constant_roi

    cpdef float get_roi(self, float z, float y, float x):
        """ Return contstant radius of influence. """
        return self.constant_roi


cdef class DistRoI(RoIFunction):
    """ Radius of influence which expands with distance from the radar. """

    cdef float z_factor, xy_factor, min_radius
    cdef int num_offsets
    cdef float[:, :] offsets

    def __init__(self, z_factor, xy_factor, min_radius, offsets):
        """ initalize. """
        cdef int i
        self.z_factor = z_factor
        self.xy_factor = xy_factor
        self.min_radius = min_radius

        self.num_offsets = len(offsets)
        # does this array need to be explicitly de-allocated when the
        # class instance is removed?
        self.offsets = cvarray(
            shape=(self.num_offsets, 3), itemsize=sizeof(float), format='f')

        for i, (z_offset, y_offset, x_offset) in enumerate(offsets):
            self.offsets[i, 0] = z_offset
            self.offsets[i, 1] = y_offset
            self.offsets[i, 2] = x_offset

    @cython.initializedcheck(False)
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef float get_roi(self, float z, float y, float x):
        """ Return the radius of influence for coordinates in meters. """
        cdef float min_roi, roi, z_offset, y_offset, x_offset
        cdef int i

        min_roi = 999999999.0
        for i in range(self.num_offsets):
            z_offset = self.offsets[i, 0]
            y_offset = self.offsets[i, 1]
            x_offset = self.offsets[i, 2]
            roi = (self.z_factor * (z - z_offset) + self.xy_factor *
                   sqrt((x - x_offset)**2 + (y - y_offset)**2))
            if roi < self.min_radius:
                roi = self.min_radius
            if roi < min_roi:
                min_roi = roi
        return min_roi


cdef class DistBeamRoI(RoIFunction):
    """
    Radius of influence which expands with distance from multiple radars.
    """

    cdef float min_radius, beam_factor
    cdef int num_offsets
    cdef float[:, :] offsets
    cdef float[::1] h_factor

    def __init__(self, float[::1] h_factor, float nb, float bsp, float min_radius, offsets):
        """ initalize. """
        cdef int i
        self.h_factor = h_factor
        self.min_radius = min_radius
        self.beam_factor = tan(nb * bsp * PI / 180.)

        self.num_offsets = len(offsets)
        # does this array need to be explicitly de-allocated when the
        # class instance is removed?
        self.offsets = cvarray(
            shape=(self.num_offsets, 3), itemsize=sizeof(float), format='f')

        for i, (z_offset, y_offset, x_offset) in enumerate(offsets):
            self.offsets[i, 0] = z_offset
            self.offsets[i, 1] = y_offset
            self.offsets[i, 2] = x_offset

    @cython.initializedcheck(False)
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef float get_roi(self, float z, float y, float x):
        """ Return the radius of influence for coordinates in meters. """

        cdef float min_roi, roi, z_offset, y_offset, x_offset
        cdef int i

        min_roi = 999999999.0
        for i in range(self.num_offsets):
            z_offset = self.offsets[i, 0]
            y_offset = self.offsets[i, 1]
            x_offset = self.offsets[i, 2]
            roi = (sqrt((self.h_factor[0] * (z - z_offset))**2 +
                        (self.h_factor[1] * (y - y_offset))**2 +
                        (self.h_factor[2] * (x - x_offset))**2) *
                   self.beam_factor
            )
            if roi < self.min_radius:
                roi = self.min_radius
            if roi < min_roi:
                min_roi = roi
        return min_roi


cdef class GateToGridMapper:
    """
    A class for efficient mapping of radar gates to a regular grid by
    weighting all gates within a specified radius of influence by distance.

    Parameters
    ----------
    grid_shape, : tuple of ints
        Shape of the grid along the z, y, and x dimensions.
    grid_starts, grid_steps  : tuple of ints
        Starting points and step sizes in meters of the grid along the
        z, y and x dimensions.
    grid_sum, grid_wsum : 4D float32 array
        Array for collecting grid weighted values and weights for each grid
        point and field. Dimension are order z, y, x, and fields. These array
        are modified in place when mapping gates unto the grid.

    """

    cdef float x_step, y_step, z_step
    cdef float x_start, y_start, z_start
    cdef int nx, ny, nz, nfields
    cdef float[:, :, :, ::1] grid_sum
    cdef float[:, :, :, ::1] grid_wsum
    cdef double[:, :, :, :] min_dist2


    def __init__(self, tuple grid_shape, tuple grid_starts, tuple grid_steps,
                 float[:, :, :, ::1] grid_sum, float[:, :, :, ::1] grid_wsum):
        """ initialize. """

        # unpack tuples
        nz, ny, nx = grid_shape
        z_start, y_start, x_start = grid_starts
        z_step, y_step, x_step = grid_steps

        # set attributes
        self.x_step = x_step
        self.y_step = y_step
        self.z_step = z_step
        self.x_start = x_start
        self.y_start = y_start
        self.z_start = z_start
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.nfields = grid_sum.shape[3]
        self.grid_sum = grid_sum
        self.grid_wsum = grid_wsum
        self.min_dist2 = 1e30*np.ones((nz, ny, nx, self.nfields))
        return

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def find_roi_for_grid(
            self, float[:, :, ::1] roi_array, RoIFunction roi_func):
        """
        Fill in the radius of influence for each point in the grid.

        Parameters
        ----------
        roi_array : 3D float32 array
            Array which will be filled by the radius of influence for each
            point in the grid.
        roi_func : RoIFunction
            Object whose get_roi method returns the radius of influence.

        """
        cdef int ix, iy, iz
        cdef float x, y, z, roi
        for ix in range(self.nx):
            for iy in range(self.ny):
                for iz in range(self.nz):
                    x = self.x_start + self.x_step * ix
                    y = self.y_start + self.y_step * iy
                    z = self.z_start + self.z_step * iz
                    roi = roi_func.get_roi(z, y, x)
                    roi_array[iz, iy, ix] = roi
        return

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def map_gates_to_grid(
            self,
            int ngates, int nrays,
            float[:, ::1] gate_z, float[:, ::1] gate_y, float[:, ::1] gate_x,
            float[:, :, ::1] field_data,
            char[:, :, ::1] field_mask, char[:, ::1] excluded_gates,
            RoIFunction roi_func, int weighting_function,
            float[::1] dist_factor):
        """
        Map radar gates unto the regular grid.

        The grid_sum and grid_wsum arrays used to initalize the class
        are update with the mapped gate data.

        Parameters
        ----------
        ngates, nrays : int
            Number of gates and rays in the radar volume.
        gate_z, gate_y, gate_x : 2D float32 array
            Cartesian locations of the gates in meters.
        field_data : 3D float32 array
            Array containing field data for the radar, dimension are ordered
            as nrays, ngates, nfields.
        field_mask : 3D uint8 array
            Array containing masking of the field data for the radar,
            dimension are ordered as nrays, ngates, nfields.
        excluded_gates : 2D uint8 array
            Array containing gate masking information. Gates with non-zero
            values will not be included in the mapping.
        offset : tuple of floats
            Offset of the radar from the grid origin. Dimension are ordered
            as z, y, x.
            Top of atmosphere. Gates above this level are considered.
        roi_func : RoIFunction
            Object whose get_roi method returns the radius of influence.
        weighting_function : int
            Function to use for weighting gates based upon distance.
            0 for Barnes, 1 for Cressman, 2 for Nearest and 3 for Barnes 2
            neighbor weighting.
        dist_factor: 3-element float32 array
            Scaling factors for squared z,y,x difference in distance calculation.
            For example:
            A value of (0.0, 1.0, 1.0)  combined with an h_factor=(0.0, 1.0, 1.0)
            (if calling DistBeamRoI) or z_factor=0.0 (if calling DistRoI) results in
            the exclusion of the z dimension in gridding weighting and could
            serve as a potential solution for gridding a single PPI sweep.

        """

        cdef float roi
        cdef float[:] values
        cdef char[:] masks
        cdef float x, y, z

        for nray in range(nrays):
            for ngate in range(ngates):

                # continue if gate excluded
                if excluded_gates[nray, ngate]:
                    continue

                x = gate_x[nray, ngate]
                y = gate_y[nray, ngate]
                z = gate_z[nray, ngate]
                roi = roi_func.get_roi(z, y, x)
                values = field_data[nray, ngate]
                masks = field_mask[nray, ngate]

                self.map_gate(x, y, z, roi, values, masks, weighting_function,
                              dist_factor)

    @cython.initializedcheck(False)
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int map_gate(self, float x, float y, float z, float roi,
                      float[:] values, char[:] masks,
                      int weighting_function, float[:] dist_factor):
        """ Map a single gate to the grid. """

        cdef float xg, yg, zg, weight, roi2, dist2, min_dist2
        cdef int x_min, x_max, y_min, y_max, z_min, z_max
        cdef int xi, yi, zi, x_argmin, y_argmin, z_argmin

        # shift positions so that grid starts at 0
        x -= self.x_start
        y -= self.y_start
        z -= self.z_start

        x_min = find_min(x, roi, self.x_step)
        if x_min > self.nx-1:
            return 0
        x_max = find_max(x, roi, self.x_step, self.nx)
        if x_max < 0:
            return 0

        y_min = find_min(y, roi, self.y_step)
        if y_min > self.ny-1:
            return 0
        y_max = find_max(y, roi, self.y_step, self.ny)
        if y_max < 0:
            return 0

        z_min = find_min(z, roi, self.z_step)
        if z_min > self.nz-1:
            return 0
        z_max = find_max(z, roi, self.z_step, self.nz)
        if z_max < 0:
            return 0

        roi2 = roi * roi

        if weighting_function == NEAREST:
            # Get the xi, yi, zi of desired weight
            x_argmin = -1
            y_argmin = -1
            z_argmin = -1
            for xi in range(x_min, x_max+1):
                for yi in range(y_min, y_max+1):
                    for zi in range(z_min, z_max+1):
                        xg = self.x_step * xi
                        yg = self.y_step * yi
                        zg = self.z_step * zi
                        dist2 = (dist_factor[2] * (xg - x)**2 +
                                 dist_factor[1] * (yg - y)**2 +
                                 dist_factor[0] * (zg - z)**2)
                        if dist2 >= roi2:
                            continue
                        for i in range(self.nfields):
                            if dist2 < self.min_dist2[zi, yi, xi, i]:
                                self.min_dist2[zi, yi, xi, i] = dist2
                                x_argmin = xi
                                y_argmin = yi
                                z_argmin = zi
                                if masks[i]:
                                    self.grid_wsum[zi, yi, xi, i] = 0
                                    self.grid_sum[zi, yi, xi, i] = 0
                                else:
                                    self.grid_wsum[z_argmin, y_argmin, x_argmin, i] = 1
                                    self.grid_sum[z_argmin, y_argmin, x_argmin, i] = values[i]
        else:
            for xi in range(x_min, x_max+1):
                for yi in range(y_min, y_max+1):
                    for zi in range(z_min, z_max+1):
                        xg = self.x_step * xi
                        yg = self.y_step * yi
                        zg = self.z_step * zi
                        dist2 = (dist_factor[2] * (xg-x)*(xg-x) +
                                 dist_factor[1] * (yg-y)*(yg-y) +
                                 dist_factor[0] * (zg-z)*(zg-z))

                        if dist2 > roi2:
                            continue

                        if weighting_function == BARNES:
                            weight = exp(-(dist2) / (2*roi2)) + 1e-5
                        elif weighting_function == BARNES2:
                            weight = exp(-(dist2) / (roi2/4)) + 1e-5
                        else: # Cressman
                            weight = (roi2 - dist2) / (roi2 + dist2)

                        for i in range(self.nfields):
                            if masks[i]:
                                continue
                            self.grid_sum[zi, yi, xi, i] += weight * values[i]
                            self.grid_wsum[zi, yi, xi, i] += weight
        return 1


@cython.cdivision(True)
cdef int find_min(float a, float roi, float step):
    """ Find the mimumum gate index for a dimension. """
    cdef int a_min
    if step == 0:
        return 0
    a_min = <int>ceil((a - roi) / step)
    if a_min < 0:
        a_min = 0
    return a_min


@cython.cdivision(True)
cdef int find_max(float a, float roi, float step, int na):
    """ Find the maximum gate index for a dimension. """
    cdef int a_max
    if step == 0:
        return 0
    a_max = <int>floor((a + roi) / step)
    if a_max > na-1:
        a_max = na-1
    return a_max
