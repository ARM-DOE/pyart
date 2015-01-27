"""
pyart.correct.region_dealias
============================

Region based dealiasing using a dynamic network reduction for region joining.

.. autosummary::
    :toctree: generated/

    dealias_region_based
    _find_regions
    _combine_regions
    _edge_sum_and_count

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    _EdgeCollector
    _RegionTracker
    _EdgeTracker

"""

import numpy as np
import scipy.ndimage as ndimage
import scipy.sparse as sparse

from ..config import get_metadata
from ._common_dealias import _parse_fields, _parse_gatefilter
from ._common_dealias import _parse_rays_wrap_around, _parse_nyquist_vel


# Possible future improvements to the region based dealiasing algorithm:
#
# * Find connected regions in the entire volume (3D) an dealias the entire
#   volume as apposed to the current sweep-by-sweep approach.
# * Scale the Doppler velocities by the Nyquist interval so that
#   calculations are done in units of percentage of the Nyquist interval.
#   This scaling is already done in the _EdgeTracker object but could also be
#   applied before finding regions of similar velocities.
# * Currently no attmpt is made to unfold gates which are included in the
#   gatefilter.  One could add a second pass to the algorithms where gates
#   excluded by the filter are unfolded using a more elemenary approach.  In
#   this manner the filter could define only "good" gates which establish the
#   folding pattern which is then applied to all gates.
# * Either each sweep is assumed to be 'centered' or the largest region
#   is assumed to not be folded.  In some cases this may not be true and
#   all the gates in the corrected sweep should be unfolded after the routine
#   completes.  By how much the sweep should be unfolded would need to be
#   determined, perhapes by comparing against sweeps above and below the
#   current sweep, or by comparing to an atmospheric sounding.
# * Improve performance by writing portions of the _edge_sum_and_count
#   function and EdgeCollector in Cython.  This step is typically the most
#   time-consuming for real world radar volumes.
# * Improve perfornace by implementing a priority queue in the _EdgeTracker
#   object. See comments in the class for details.


def dealias_region_based(
        radar, interval_splits=3, interval_limits=None,
        skip_between_rays=2, skip_along_ray=2, centered=True,
        nyquist_vel=None, gatefilter=None, rays_wrap_around=None,
        keep_original=True, vel_field=None, corr_vel_field=None, **kwargs):
    """
    Dealias Doppler velocities using a region based algorithm.

    Performs Doppler velocity dealiasing by finding regions of similar
    velocities and unfolding and merging pairs of regions until all
    regions are unfolded.  Unfolding and merging regions is accomplished by
    modeling the problem as a dynamic network reduction.

    Parameters
    ----------
    radar : Radar
        Radar object containing Doppler velocities to dealias.
    interval_splits : int, optional
        Number of segments to split the nyquist interval into when finding
        regions of similar velocity.  More splits creates a larger number of
        initial regions which takes longer to process but may result in better
        dealiasing.  The default value of 3 seems to be a good compromise
        between performance and artifact free dealiasing.  This value
        is not used if the interval_limits parameter is not None.
    interval_limits : array like or None, optional
        Velocity limits used for finding regions of similar velocity.  Should
        cover the entire nyquist interval.  None, the default value, will
        split the Nyquist interval into interval_splits equal sized
        intervals.
    skip_between_rays, skip_along_ray : int, optional
        Maximum number of filtered gates to skip over when joining regions,
        gaps between region larger than this will not be connected.  Parameters
        specify the maximum number of filtered gates between and along a ray.
        Set these parameters to 0 to disable unfolding across filtered gates.
    centered : bool, optional
        True to apply centering to each sweep after the dealiasing algorithm
        so that the average number of unfolding is near 0.  False does not
        apply centering which may results in individual sweeps under or over
        folded by the nyquist interval.
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

    # find nyquist interval segmentation limits
    if interval_limits is None:
        interval_limits = np.linspace(
            -nyquist_vel, nyquist_vel, interval_splits+1, endpoint=True)

    # exclude masked and invalid velocity gates
    gatefilter.exclude_masked(vel_field)
    gatefilter.exclude_invalid(vel_field)
    gfilter = gatefilter.gate_excluded

    # perform dealiasing
    vdata = radar.fields[vel_field]['data'].view(np.ndarray)
    data = vdata.copy()     # dealiased velocities

    for sweep_slice in radar.iter_slice():      # loop over sweeps

        # extract sweep data
        sdata = vdata[sweep_slice].copy()   # is a copy needed here?
        scorr = data[sweep_slice]
        sfilter = gfilter[sweep_slice]

        # find regions in original data
        labels, nfeatures = _find_regions(sdata, sfilter, interval_limits)
        edge_sum, edge_count, region_sizes = _edge_sum_and_count(
            labels, nfeatures, sdata, rays_wrap_around, skip_between_rays,
            skip_along_ray)

        # find the number of folds in the regions
        region_tracker = _RegionTracker(region_sizes)
        edge_tracker = _EdgeTracker(edge_sum, edge_count, nyquist_interval)
        while True:
            if _combine_regions(region_tracker, edge_tracker):
                break

        # center sweep if requested, determine a global sweep unfold number
        # so that the average number of gate folds is zero.
        if centered:
            gates_dealiased = region_sizes.sum()
            total_folds = np.sum(
                region_sizes * region_tracker.unwrap_number[1:])
            sweep_offset = int(round(float(total_folds) / gates_dealiased))
            if sweep_offset != 0:
                region_tracker.unwrap_number -= sweep_offset

        # dealias the data using the fold numbers
        # start from label 1 to skip masked region
        for i in range(1, nfeatures+1):
            nwrap = region_tracker.unwrap_number[i]
            if nwrap != 0:
                scorr[labels == i] += nwrap * nyquist_interval

    # mask filtered gates
    if np.any(gfilter):
        data = np.ma.array(data, mask=gfilter)

    # restore original values where dealiasing not applied
    if keep_original:
        data[gfilter] = radar.fields[vel_field]['data'][gfilter]

    # return field dictionary containing dealiased Doppler velocities
    corr_vel = get_metadata(corr_vel_field)
    corr_vel['data'] = data
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
        Velocity limits for region finding.  For each pair of limits, taken
        from elements i and i+1 of the array, all connected regions with
        velocities within these limits will be found.

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


# This function takes considerable time and would be a good canidate for
# rewriting in Cython for better performance.
def _edge_sum_and_count(labels, nfeatures, data,
                        rays_wrap_around, max_gap_x, max_gap_y):
    """
    Return sparse matrices containing the edge sums and counts between regions
    and a array of the number of gates in each region.
    """

    bincount = np.bincount(labels.ravel())
    num_masked_gates = bincount[0]
    region_sizes = bincount[1:]

    total_nodes = np.prod(labels.shape) - num_masked_gates
    if rays_wrap_around:
        total_nodes += labels.shape[0] * 2
    collector = _EdgeCollector(total_nodes)
    right, bottom = [i-1 for i in labels.shape]

    for index, label in np.ndenumerate(labels):
        if label == 0:
            continue

        x_index, y_index = index
        vel = data[x_index, y_index]

        # left
        x_check = x_index - 1
        if x_check == -1 and rays_wrap_around:
            x_check = right     # wrap around
        if x_check != -1:
            neighbor = labels[x_check, y_index]

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
                    if neighbor != 0:
                        break

            # add the edge to the collection (if valid)
            collector.add_edge(label, neighbor, vel)

        # right
        x_check = x_index + 1
        if x_check == right+1 and rays_wrap_around:
            x_check = 0     # wrap around
        if x_check != right+1:
            neighbor = labels[x_check, y_index]

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
                    if neighbor != 0:
                        break

            # add the edge to the collection (if valid)
            collector.add_edge(label, neighbor, vel)

        # top
        y_check = y_index - 1
        if y_check != -1:
            neighbor = labels[x_index, y_check]

            # if the top side gate is masked, keep looking up
            # until we find a valid gate or reach the maximum gap size
            if neighbor == 0:
                for i in range(max_gap_y):
                    y_check -= 1
                    if y_check == -1:
                        break
                    neighbor = labels[x_index, y_check]
                    if neighbor != 0:
                        break

            # add the edge to the collection (if valid)
            collector.add_edge(label, neighbor, vel)

        # bottom
        y_check = y_index + 1
        if y_check != bottom + 1:
            neighbor = labels[x_index, y_check]

            # if the top side gate is masked, keep looking up
            # until we find a valid gate or reach the maximum gap size
            if neighbor == 0:
                for i in range(max_gap_y):
                    y_check += 1
                    if y_check == bottom + 1:
                        break
                    neighbor = labels[x_index, y_check]
                    if neighbor != 0:
                        break

            # add the edge to the collection (if valid)
            collector.add_edge(label, neighbor, vel)

    edge_sum, edge_count = collector.make_edge_matrices(nfeatures)
    return edge_sum, edge_count, region_sizes


class _EdgeCollector(object):
    """
    Class for collecting edges, used by _edge_sum_and_count function.
    """

    def __init__(self, total_nodes):
        """ initalize. """
        self.l_index = np.zeros(total_nodes * 4, dtype=np.int32)
        self.n_index = np.zeros(total_nodes * 4, dtype=np.int32)
        self.e_sum = np.zeros(total_nodes * 4, dtype=np.float64)
        self.idx = 0

    def add_edge(self, label, neighbor, vel):
        """ Add an edge. """
        if neighbor == label or neighbor == 0:
            # Do not add edges between the same region (circular edges)
            # or edges to masked gates (indicated by a label of 0).
            return
        self.l_index[self.idx] = label
        self.n_index[self.idx] = neighbor
        self.e_sum[self.idx] = vel
        self.idx += 1
        return

    def make_edge_matrices(self, nfeatures):
        """ Return sparse matrices for the edge sums and counts. """
        matrix_indices = (self.l_index[:self.idx], self.n_index[:self.idx])
        e_sum = self.e_sum[:self.idx]
        e_count = np.ones_like(e_sum, dtype=np.int32)
        shape = (nfeatures+1, nfeatures+1)
        edge_sum = sparse.coo_matrix((e_sum, matrix_indices), shape=shape)
        edge_count = sparse.coo_matrix((e_count, matrix_indices), shape=shape)
        return edge_sum.tocsr(), edge_count.tocsr()


def _combine_regions(region_tracker, edge_tracker):
    """ Returns True when done. """
    # Edge parameters from edge with largest weight
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

        # array of linked lists pointing to each node
        self.edges_in_node = np.zeros(nnodes, dtype='object')
        for i in range(nnodes):
            self.edges_in_node[i] = []

        # fill out data from edge_count and edge_sum
        edge = 0
        for i, j in zip(*edge_count.nonzero()):
            if i < j:
                continue
            assert edge_count[i, j] == edge_count[j, i]
            self.node_alpha[edge] = i
            self.node_beta[edge] = j
            self.sum_diff[edge] = ((edge_sum[i, j] - edge_sum[j, i]) /
                                   nyquist_interval)
            self.weight[edge] = edge_count[i, j]
            self.edges_in_node[i].append(edge)
            self.edges_in_node[j].append(edge)

            edge += 1

        # list which orders edges according to their weight, highest first
        # TODO
        self.priority_queue = []

    def merge_nodes(self, base_node, merge_node, foo_edge):
        """ Merge nodes. """

        # remove edge between base and merge nodes
        self.weight[foo_edge] = -999
        self.edges_in_node[merge_node].remove(foo_edge)
        self.edges_in_node[base_node].remove(foo_edge)
        self._common_finder[merge_node] = False

        # find all the edges in the two nodes
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
                                    merge_node, neighbor)
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
                       merge_node, neighbor_node):
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
        # add weight * nwrap to each edge in node
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
