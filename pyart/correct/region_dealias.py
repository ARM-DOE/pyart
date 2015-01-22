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
# * loop over sweeps
# * mask output if needed
# * replace masked with original vel
# * nyquist free edge tracking?
# * gap jumping edge connections?
# * 3D segmentation?
# * optimize
# * refactor
# * document
# * unit tests
# * merge into Py-ART


def dealias_region_based(
        radar, segmentation_splits=3, segmentation_limits=None,
        nyquist_vel=None, gatefilter=None,
        rays_wrap_around=None, keep_original=True, vel_field=None,
        corr_vel_field=None, **kwargs):
    """
    Dealias Doppler velocities using a region based algorithm.

    Parameters
    ----------
    radar : Radar
        Radar object containing Doppler velocities to dealias.
    segmentation_splits : int, optional
        Number of segments to split the nyquist interval into when finding
        regions of similar velocity.  More splits creates a larger number of
        initial regions which takes longer to process but may result in better
        dealiasing.  The default value of 3 seems to be a good compromise
        between performance and artifact free dealiasing.  This value
        is not used if the segmentation_limits parameter is not None.
    segmentation_limits : array like or None, optional
        Velocity limits used for finding regions of similar velocity.  Should
        cover the entire nyquist interval.  None, the default value, will
        split the nyquist interval into segmentation_splits equal sized
        intervals.
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
    nyquist_interval = 2. * nyquist_vel

    # find segmentation limits
    if segmentation_limits is None:
        segmentation_limits = np.linspace(
            -nyquist_vel, nyquist_vel, segmentation_splits+1, endpoint=True)

    # exclude masked and invalid velocity gates
    gatefilter.exclude_masked(vel_field)
    gatefilter.exclude_invalid(vel_field)
    gfilter = gatefilter.gate_excluded

    # raw velocity data possibly with masking
    vdata = radar.fields[vel_field]['data'].copy()

    # find regions/segments in original data
    labels, nfeatures = _find_regions(vdata, gfilter, segmentation_limits)
    edge_sum, edge_count, segment_sizes = _edge_sum_and_count(
        labels, nfeatures, vdata, rays_wrap_around)

    # find unwrap number for these regions
    region_tracker = _RegionTracker(segment_sizes)
    edge_tracker = _EdgeTracker(edge_sum, edge_count, nyquist_interval)
    while True:
        if _combine_segments(region_tracker, edge_tracker):
            break

    # dealias the data using the unwrap numbers
    dealias_data = vdata.copy()
    for i in range(1, nfeatures+1):     # start from 0 to skip masked region
        nwrap = region_tracker.unwrap_number[i]
        if nwrap != 0:
            dealias_data[labels == i] += nwrap * nyquist_interval

    # return field dictionary containing dealiased Doppler velocities
    corr_vel = get_metadata(corr_vel_field)
    corr_vel['data'] = dealias_data
    return corr_vel


def _find_regions(vel, gfilter, limits):
    """
    Find regions of similar velocity.

    For each pair of values in the limits array (or list) find all connected
    velocity regions within these limits.

    Parameters
    ----------
    vel : 2D ndarray
        Array containing velocity data for a single sweep.
    gfilter : 2D ndarray
        Filter indicating if a particular gate should be masked.  True
        indicates the gate should be masked (excluded).
    limits : array like
        Velocity limits for segmentation.  For each pair of limits, taken from
        elements i and i+1 of the array, all connected regions with velocities
        within these limits will be found.

    Returns
    -------
    label : ndarray
        Interger array with each region labeled by a value.  The array
        ranges from 0 to nfeatures, inclusive, where a value of 0 indicates
        masked gates and non-zero indicates a region of connected gates.
    nfeatures : int
        Number of regions found.

    """
    mask = ~gfilter
    label = np.zeros(vel.shape, dtype=np.int32)
    nfeatures = 0
    for lmin, lmax in zip(limits[:-1], limits[1:]):

        # find connected regions within the limits
        inp = (lmin <= vel) & (vel < lmax) & mask
        limit_label, limit_nfeatures = ndimage.label(inp)

        # add these regions to the global regions
        limit_label[np.nonzero(limit_label)] += nfeatures
        label += limit_label
        nfeatures += limit_nfeatures

    return label, nfeatures


def _combine_segments(region_tracker, edge_tracker):
    """ Returns True when done. """
    # Edge parameter from edge with largest weight
    status, extra = edge_tracker.pop_edge()
    if status:
        return True
    node1, node2, _, diff, edge_number = extra
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


def _edge_sum_and_count(labels, nfeatures, data, rays_wrap_around):
    """ Return a sparse matrices containing the edge sums and counts. """

    bincount = np.bincount(labels.ravel())
    num_masked_gates = bincount[0]
    segment_sizes = bincount[1:]

    total_nodes = np.prod(labels.shape) - num_masked_gates
    if rays_wrap_around:
        total_nodes += labels.shape[0] * 2
    l_index = np.zeros(total_nodes * 4, dtype=np.int32)
    n_index = np.zeros(total_nodes * 4, dtype=np.int32)
    e_sum = np.zeros(total_nodes * 4, dtype=np.float64)
    right, bottom = [i-1 for i in labels.shape]

    idx = 0
    for index, label in np.ndenumerate(labels):
        if label == 0:
            continue

        x_index, y_index = index
        vel = data[x_index, y_index]

        # left
        if x_index != 0:
            neighbor = labels[x_index-1, y_index]
            if neighbor != label and neighbor != 0:
                l_index[idx] = label
                n_index[idx] = neighbor
                e_sum[idx] = vel
                idx += 1
        elif rays_wrap_around:
            neighbor = labels[right, y_index]
            if neighbor != label and neighbor != 0:
                l_index[idx] = label
                n_index[idx] = neighbor
                e_sum[idx] = vel
                idx += 1

        # right
        if x_index != right:
            neighbor = labels[x_index+1, y_index]
            if neighbor != label and neighbor != 0:
                l_index[idx] = label
                n_index[idx] = neighbor
                e_sum[idx] = vel
                idx += 1
        elif rays_wrap_around:
            neighbor = labels[0, y_index]
            if neighbor != label and neighbor != 0:
                l_index[idx] = label
                n_index[idx] = neighbor
                e_sum[idx] = vel
                idx += 1

        # top
        if y_index != 0:
            neighbor = labels[x_index, y_index-1]
            if neighbor != label and neighbor != 0:
                l_index[idx] = label
                n_index[idx] = neighbor
                e_sum[idx] = vel
                idx += 1

        # bottom
        if y_index != bottom:
            neighbor = labels[x_index, y_index+1]
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

    return edge_sum.tocsr(), edge_count.tocsr(), segment_sizes


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
    """ A class for tracking edges in a dynamic network. """

    def __init__(self, edge_sum, edge_count, nyquist_interval):
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
        edge = 0
        for i, j in zip(*edge_count.nonzero()):
            if i < j:
                continue

            self.node_alpha[edge] = i
            self.node_beta[edge] = j
            self.sum_diff[edge] = ((edge_sum[i, j] - edge_sum[j, i]) /
                                   nyquist_interval)
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
        old_alpha = int(self.node_alpha[edge])
        old_beta = int(self.node_beta[edge])
        self.node_alpha[edge] = old_beta
        self.node_beta[edge] = old_alpha
        # swap sums
        self.sum_diff[edge] = -1. * self.sum_diff[edge]
        return

    def unwrap_node(self, node, nwrap):
        """ Unwrap a node. """
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
