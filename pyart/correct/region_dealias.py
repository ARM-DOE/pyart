"""
pyart.correct.region_dealias
============================

Region based dealiasing using a dynamic network reduction for region joining.

.. autosummary::
    :toctree: generated/

    dealias_region_based


"""

import numpy as np
import scipy.ndimage as ndimage
import scipy.sparse as sparse

from ..config import get_metadata

# XXX move these to a _common_dealias module
from .unwrap import _parse_fields, _parse_gatefilter, _parse_rays_wrap_around
from .unwrap import _parse_nyquist_vel


# TODO
# * segmentation parameter based number of splits
# * periodic ray boundary (respect rays_wrap_around), only in _edge_sum_and_c
# * loop over sweeps
# * mask output if needed
# * replace masked with original vel
# * optimize
# * refactor
# * document
# * unit tests
# * merge into Py-ART


def dealias_region_based(
        radar, nyquist_vel=None, gatefilter=None,
        rays_wrap_around=None, keep_original=True, vel_field=None,
        corr_vel_field=None, **kwargs):
    """
    Dealias Doppler velocities using a region based algorithm.

    Parameters
    ----------
    radar : Radar
        Radar object containing Doppler velocities to dealias.
    nyquist_velocity : float, optional
        Nyquist velocity in unit identical to those stored in the radar's
        velocity field.  None will attempt to determine this value from the
        instrument_parameters attribute.
    gatefilter : GateFilter, None or False, optional.
        A GateFilter instance which specified which gates should be
        ignored when performing de-aliasing.  A value of None, the default,
        created this filter from the radar moments using any additional
        arguments by passing them to :py:func:`moment_based_gate_filter`.
        False disables filtering including all gates in the dealiasing.
    rays_wrap_around : bool or None, optional
        True when the rays at the beginning of the sweep and end of the sweep
        should be interpreted as connected when de-aliasing (PPI scans).
        False if they edges should not be interpreted as connected (other scan
        types).  None will determine the correct value from the radar
        scan type.
    keep_original : bool, optional
        True to retain the original Doppler velocity values at gates
        where the dealiasing procedure fails or was not applied. False
        does not replacement and these gates will be masked in the corrected
        velocity field.
    vel_field : str, optional
        Field in radar to use as the Doppler velocities during dealiasing.
        None will use the default field name from the Py-ART configuration
        file.
    corr_vel_field : str, optional
        Name to use for the dealiased Doppler velocity field metadata.  None
        will use the default field name from the Py-ART configuration file.

    Returns
    -------
    corr_vel : dict
        Field dictionary containing dealiased Doppler velocities.  Dealiased
        array is stored under the 'data' key.

    """
    # parse function parameters
    vel_field, corr_vel_field = _parse_fields(vel_field, corr_vel_field)
    gatefilter = _parse_gatefilter(gatefilter, radar, **kwargs)
    rays_wrap_around = _parse_rays_wrap_around(rays_wrap_around, radar)
    nyquist_vel = _parse_nyquist_vel(nyquist_vel, radar)

    # exclude masked and invalid velocity gates
    gatefilter.exclude_masked(vel_field)
    gatefilter.exclude_invalid(vel_field)
    gfilter = gatefilter.gate_excluded

    # raw velocity data possibly with masking
    vdata = radar.fields[vel_field]['data'].copy()

    # find regions/segments in original data
    labels, nfeatures = _find_segments(vdata, nyquist_vel, gfilter)
    masked_gates, segment_sizes = _segment_sizes(labels, nfeatures)
    edge_sum, edge_count = _edge_sum_and_count(labels, nfeatures, vdata,
                                               masked_gates)

    # find unwrap number for these regions
    region_tracker = _RegionTracker(segment_sizes)
    edge_tracker = _EdgeTracker(edge_sum, edge_count, nyquist_vel)
    while True:
        if _combine_segments(region_tracker, edge_tracker):
            break

    # dealias the data using the unwrap numbers
    dealias_data = vdata.copy()
    for i in range(1, nfeatures+1):     # start from 0 to skip masked region
        nwrap = region_tracker.unwrap_number[i]
        dealias_data[labels == i] += nwrap * nyquist_vel

    # return field dictionary containing dealiased Doppler velocities
    corr_vel = get_metadata(corr_vel_field)
    corr_vel['data'] = dealias_data
    return corr_vel


def _find_segments(vel, nyquist, gfilter):
    """ Find segments """
    # label 0 indicates masked gates

    mask = ~gfilter
    inp = (vel > nyquist * 1/3.)
    inp = np.logical_and(inp, mask)
    label, nfeatures = ndimage.label(inp)

    inp = (vel < -nyquist * 1/3.)
    inp = np.logical_and(inp, mask)
    labels2, nfeatures2 = ndimage.label(inp)
    labels2[np.nonzero(labels2)] += nfeatures

    inp = np.logical_and(vel <= nyquist * 1/3., vel >= -nyquist * 1/3.)
    inp = np.logical_and(inp, mask)
    labels3, nfeatures3 = ndimage.label(inp)
    labels3[np.nonzero(labels3)] += nfeatures + nfeatures2

    return label + labels2 + labels3, nfeatures + nfeatures2 + nfeatures3


def _combine_segments(region_tracker, edge_tracker):
    """ Returns True when done. """
    # Edge parameter from edge with largest weight
    status, extra = edge_tracker.pop_edge()
    if status:
        return True
    node1, node2, weight, diff, edge_number = extra
    rdiff = np.round(diff)

    # node sizes of nodes to be merged
    node1_size = region_tracker.get_node_size(node1)
    node2_size = region_tracker.get_node_size(node2)

    # determine which nodes should be merged
    if node1_size > node2_size:
        base_node, merge_node = node1, node2
    else:
        base_node, merge_node = node2, node1
        rdiff = -rdiff

    # unwrap merge_node
    if rdiff != 0:
        region_tracker.unwrap_node(merge_node, rdiff)
        edge_tracker.unwrap_node(merge_node, rdiff)

    # merge nodes
    region_tracker.merge_nodes(base_node, merge_node)
    edge_tracker.merge_nodes(base_node, merge_node, edge_number)

    return False


def _edge_sum_and_count(labels, nfeatures, data, masked_gates):

    total_nodes = np.prod(labels.shape) - masked_gates
    l_index = np.zeros(total_nodes * 4, dtype=np.int32)
    n_index = np.zeros(total_nodes * 4, dtype=np.int32)
    e_sum = np.zeros(total_nodes * 4, dtype=np.float64)
    right, bottom = [i-1 for i in labels.shape]

    idx = 0
    for index, label in np.ndenumerate(labels):
        if label == 0:
            continue

        x, y = index
        vel = data[x, y]

        # left
        if x != 0:
            neighbor = labels[x-1, y]
            if neighbor != label and neighbor != 0:
                l_index[idx] = label
                n_index[idx] = neighbor
                e_sum[idx] = vel
                idx += 1

        # right
        if x != right:
            neighbor = labels[x+1, y]
            if neighbor != label and neighbor != 0:
                l_index[idx] = label
                n_index[idx] = neighbor
                e_sum[idx] = vel
                idx += 1

        # top
        if y != 0:
            neighbor = labels[x, y-1]
            if neighbor != label and neighbor != 0:
                l_index[idx] = label
                n_index[idx] = neighbor
                e_sum[idx] = vel
                idx += 1

        # bottom
        if y != bottom:
            neighbor = labels[x, y+1]
            if neighbor != label and neighbor != 0:
                l_index[idx] = label
                n_index[idx] = neighbor
                e_sum[idx] = vel
                idx += 1

    l_index = l_index[:idx]
    n_index = n_index[:idx]
    e_sum = e_sum[:idx]
    e_count = np.ones_like(e_sum, dtype=np.int32)

    shape = (nfeatures+1, nfeatures+1)
    edge_sum = sparse.coo_matrix((e_sum, (l_index, n_index)), shape=shape)
    edge_count = sparse.coo_matrix((e_count, (l_index, n_index)), shape=shape)

    return edge_sum.tocsr(), edge_count.tocsr()




def _segment_sizes(labels, nfeatures):
    x =  np.bincount(labels.ravel())
    return x[0], x[1:]

class _RegionTracker(object):
    """
    Tracks the location of radar volume regions contained in each node
    as the network is reduced.
    """

    def __init__(self, region_sizes):
        """ initalize. """
        # number of gates in each node
        nregions = len(region_sizes) + 1
        self.node_size = np.zeros(nregions, dtype='int32')
        self.node_size[1:] = region_sizes[:]

        # array of lists containing the regions in each node
        self.regions_in_node = np.zeros(nregions, dtype='object')
        for i in range(nregions):
            self.regions_in_node[i] = [i]

        # number of unwrappings to apply to dealias each region
        self.unwrap_number = np.zeros(nregions, dtype='int32')

    def merge_nodes(self, node_a, node_b):
        """ Merge node b into node a. """

        # move all regions from node_b to node_a
        regions_to_merge = self.regions_in_node[node_b]
        self.regions_in_node[node_a].extend(regions_to_merge)
        self.regions_in_node[node_b] = []

        # update node sizes
        self.node_size[node_a] += self.node_size[node_b]
        self.node_size[node_b] = 0
        return

    def unwrap_node(self, node, nwrap):
        """ Unwrap all gates contained a node. """
        if nwrap == 0:
            return
        # for each region in node add nwrap
        regions_to_unwrap = self.regions_in_node[node]
        self.unwrap_number[regions_to_unwrap] += nwrap
        return

    def get_node_size(self, node):
        """ Return the number of gates in a node. """
        return self.node_size[node]


class _EdgeTracker(object):

    def __init__(self, edge_sum, edge_count, nyquist):
        """ initialize """

        nedges = int(edge_count.nnz / 2)
        nnodes = edge_count.shape[0]

        # node number and different in sum for each edge
        self.node_alpha = np.zeros(nedges, dtype=np.int32)
        self.node_beta = np.zeros(nedges, dtype=np.int32)
        self.sum_diff = np.zeros(nedges, dtype=np.float32)

        # number of connections between the regions
        self.weight = np.zeros(nedges, dtype=np.int32)

        # fast finding
        self._common_finder = np.zeros(nnodes, dtype=np.bool)
        self._common_index = np.zeros(nnodes, dtype=np.int32)
        self._last_base_node = -1

        # array of linked lists pointing to each segment
        self.edges_in_node = np.zeros(nnodes, dtype='object')
        for i in range(nnodes):
            self.edges_in_node[i] = []

        # fill out data from edge_count and edge_sum
        a, b = edge_count.nonzero()
        edge = 0
        for i, j in zip(a, b):
            if i < j:
                continue

            self.node_alpha[edge] = i
            self.node_beta[edge] = j
            self.sum_diff[edge] = (edge_sum[i, j] - edge_sum[j, i]) / nyquist
            self.weight[edge] = edge_count[i, j]
            self.edges_in_node[i].append(edge)
            self.edges_in_node[j].append(edge)

            edge += 1

        # list which orders edges according to their
        # weight, highest first
        self.priority_queue = []

    def merge_nodes(self, base_node, merge_node, foo_edge):
        """ Merge nodes. """

        # remove edge between base and merge nodes
        self.weight[foo_edge] = -999
        self.edges_in_node[merge_node].remove(foo_edge)
        self.edges_in_node[base_node].remove(foo_edge)
        self._common_finder[merge_node] = False

        # find all the edges in the two segments
        edges_in_merge = list(self.edges_in_node[merge_node])

        # loop over base_node edges if last base_node was different
        if self._last_base_node != base_node:
            self._common_finder[:] = False
            edges_in_base = list(self.edges_in_node[base_node])
            for edge_num in edges_in_base:

                # reverse edge if needed so node_alpha is base_node
                if self.node_beta[edge_num] == base_node:
                    self._reverse_edge_direction(edge_num)
                assert self.node_alpha[edge_num] == base_node

                # find all neighboring nodes to base_node
                neighbor = self.node_beta[edge_num]
                self._common_finder[neighbor] = True
                self._common_index[neighbor] = edge_num

        # loop over edge nodes
        for edge_num in edges_in_merge:

            # reverse edge so that node alpha is the merge_node
            if self.node_beta[edge_num] == merge_node:
                self._reverse_edge_direction(edge_num)
            assert self.node_alpha[edge_num] == merge_node

            # update all the edges to point to the base node
            self.node_alpha[edge_num] = base_node

            # if base_node also has an edge with the neighbor combine them
            neighbor = self.node_beta[edge_num]
            if self._common_finder[neighbor]:
                base_edge_num = self._common_index[neighbor]
                self._combine_edges(base_edge_num, edge_num,
                                    base_node, merge_node, neighbor)
            # if not fill in _common_ arrays.
            else:
                self._common_finder[neighbor] = True
                self._common_index[neighbor] = edge_num

        # move all edges from merge_node to base_node
        edges = self.edges_in_node[merge_node]
        self.edges_in_node[base_node].extend(edges)
        self.edges_in_node[merge_node] = []
        self._last_base_node = int(base_node)
        return

    def _combine_edges(self, base_edge, merge_edge,
                       base_node, merge_node, neighbor_node):
        """ Combine edges into a single edge.  """
        # Merging nodes MUST be set to alpha prior to calling this function

        # combine edge weights
        self.weight[base_edge] += self.weight[merge_edge]
        self.weight[merge_edge] = -999.

        # combine sums
        self.sum_diff[base_edge] += self.sum_diff[merge_edge]

        # remove merge_edge from both node lists
        self.edges_in_node[merge_node].remove(merge_edge)
        self.edges_in_node[neighbor_node].remove(merge_edge)

        # TODO priority queue
        # remove merge_edge from edge priority queue
        # self.priority_queue.remove(merge_edge)

        # relocate base_edge in priority queue
        # self.priority_queue.sort()

    def _reverse_edge_direction(self, edge):
        """ Reverse an edges direction, change alpha and beta. """
        # swap nodes
        a = int(self.node_alpha[edge])
        b = int(self.node_beta[edge])
        self.node_alpha[edge] = b
        self.node_beta[edge] = a
        # swap sums
        self.sum_diff[edge] = -1. * self.sum_diff[edge]
        return

    def unwrap_node(self, node, nwrap):
        if nwrap == 0:
            return
        # add weight * nwrap to each edge in segment
        for edge in self.edges_in_node[node]:
            weight = self.weight[edge]
            if node == self.node_alpha[edge]:
                self.sum_diff[edge] += weight * nwrap
            else:
                assert self.node_beta[edge] == node
                self.sum_diff[edge] += -weight * nwrap
        return

    def pop_edge(self):
        """ Pop edge with largest weight.  Return node numbers and diff """

        # if len(priority_queue) == 0:
        #     return True, None

        # edge_num = self.priority_queue[0]
        edge_num = np.argmax(self.weight)
        node1 = self.node_alpha[edge_num]
        node2 = self.node_beta[edge_num]
        weight = self.weight[edge_num]
        diff = self.sum_diff[edge_num] / (float(weight))

        if weight < 0:
            return True, None
        return False, (node1, node2, weight, diff, edge_num)
