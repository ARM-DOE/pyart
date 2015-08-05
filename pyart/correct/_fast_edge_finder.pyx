"""
pyart.correct._fast_edge_finder
===============================

Cython routine for quickly finding edges between connected regions.

.. autosummary::
    :toctree: generated/

    _fast_edge_finder

"""

import numpy as np

cimport numpy as np
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
def _fast_edge_finder(
        int[:, ::1] labels, float[:, ::1] data, int rays_wrap_around,
        int max_gap_x, int max_gap_y, int total_nodes):
    """
    Return the gate indices and velocities of all edges between regions.
    """
    cdef int x_index, y_index, right, bottom, y_check, x_check
    cdef int label, neighbor
    cdef float vel, nvel

    collector = _EdgeCollector(total_nodes)
    right = labels.shape[0] - 1
    bottom = labels.shape[1] - 1

    for x_index in range(labels.shape[0]):
        for y_index in range(labels.shape[1]):

            label = labels[x_index, y_index]
            if label == 0:
                continue

            vel = data[x_index, y_index]

            # left
            x_check = x_index - 1
            if x_check == -1 and rays_wrap_around:
                x_check = right     # wrap around
            if x_check != -1:
                neighbor = labels[x_check, y_index]
                nvel = data[x_check, y_index]

                # if the left side gate is masked, keep looking to the left
                # until we find a valid gate or reach the maximum gap size
                if neighbor == 0:
                    for i in range(max_gap_x):
                        x_check -= 1
                        if x_check == -1:
                            if rays_wrap_around:
                                x_check = right
                            else:
                                break
                        neighbor = labels[x_check, y_index]
                        nvel = data[x_check, y_index]
                        if neighbor != 0:
                            break

                # add the edge to the collection (if valid)
                collector.add_edge(label, neighbor, vel, nvel)

            # right
            x_check = x_index + 1
            if x_check == right+1 and rays_wrap_around:
                x_check = 0     # wrap around
            if x_check != right+1:
                neighbor = labels[x_check, y_index]
                nvel = data[x_check, y_index]

                # if the right side gate is masked, keep looking to the left
                # until we find a valid gate or reach the maximum gap size
                if neighbor == 0:
                    for i in range(max_gap_x):
                        x_check += 1
                        if x_check == right+1:
                            if rays_wrap_around:
                                x_check = 0
                            else:
                                break
                        neighbor = labels[x_check, y_index]
                        nvel = data[x_check, y_index]
                        if neighbor != 0:
                            break

                # add the edge to the collection (if valid)
                collector.add_edge(label, neighbor, vel, nvel)

            # top
            y_check = y_index - 1
            if y_check != -1:
                neighbor = labels[x_index, y_check]
                nvel = data[x_index, y_check]

                # if the top side gate is masked, keep looking up
                # until we find a valid gate or reach the maximum gap size
                if neighbor == 0:
                    for i in range(max_gap_y):
                        y_check -= 1
                        if y_check == -1:
                            break
                        neighbor = labels[x_index, y_check]
                        nvel = data[x_index, y_check]
                        if neighbor != 0:
                            break

                # add the edge to the collection (if valid)
                collector.add_edge(label, neighbor, vel, nvel)

            # bottom
            y_check = y_index + 1
            if y_check != bottom + 1:
                neighbor = labels[x_index, y_check]
                nvel = data[x_index, y_check]

                # if the top side gate is masked, keep looking up
                # until we find a valid gate or reach the maximum gap size
                if neighbor == 0:
                    for i in range(max_gap_y):
                        y_check += 1
                        if y_check == bottom + 1:
                            break
                        neighbor = labels[x_index, y_check]
                        nvel = data[x_index, y_check]
                        if neighbor != 0:
                            break

                # add the edge to the collection (if valid)
                collector.add_edge(label, neighbor, vel, nvel)

    indices, velocities = collector.get_indices_and_velocities()
    return indices, velocities


# Cython implementation inspired by coo_entries in scipy/spatial/ckdtree.pyx
cdef class _EdgeCollector:
    """
    Class for collecting edges, used by _edge_sum_and_count function.
    """

    cdef np.ndarray l_index, n_index, l_velo, n_velo
    cdef np.int32_t *l_data
    cdef np.int32_t *n_data
    cdef np.float64_t *lv_data
    cdef np.float64_t *nv_data
    cdef int idx

    def __init__(self, total_nodes):
        """ initalize. """
        self.l_index = np.zeros(total_nodes * 4, dtype=np.int32)
        self.n_index = np.zeros(total_nodes * 4, dtype=np.int32)
        self.l_velo = np.zeros(total_nodes * 4, dtype=np.float64)
        self.n_velo = np.zeros(total_nodes * 4, dtype=np.float64)

        self.l_data = <np.int32_t *>np.PyArray_DATA(self.l_index)
        self.n_data = <np.int32_t *>np.PyArray_DATA(self.n_index)
        self.lv_data = <np.float64_t*>np.PyArray_DATA(self.l_velo)
        self.nv_data = <np.float64_t*>np.PyArray_DATA(self.n_velo)

        self.idx = 0

    cdef int add_edge(_EdgeCollector self, int label, int neighbor,
                      float vel, float nvel):
        """ Add an edge. """
        if neighbor == label or neighbor == 0:
            # Do not add edges between the same region (circular edges)
            # or edges to masked gates (indicated by a label of 0).
            return 0
        self.l_data[self.idx] = label
        self.n_data[self.idx] = neighbor
        self.lv_data[self.idx] = vel
        self.nv_data[self.idx] = nvel
        self.idx += 1
        return 1

    def get_indices_and_velocities(self):
        """ Return the edge indices and velocities. """
        indices = (self.l_index[:self.idx], self.n_index[:self.idx])
        velocities = (self.l_velo[:self.idx], self.n_velo[:self.idx])
        return indices, velocities
