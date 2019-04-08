"""
pyart.correct.region_dealias
============================

Region based dealiasing using a dynamic network reduction for region joining.

.. autosummary::
    :toctree: generated/

    dealias_region_based
    _find_regions
    _find_sweep_interval_splits
    _combine_regions
    _edge_sum_and_count

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    _RegionTracker
    _EdgeTracker

"""

import warnings

import numpy as np
import scipy.ndimage as ndimage

from scipy.optimize import fmin_l_bfgs_b

from ..config import get_metadata, get_fillvalue
from ._common_dealias import _parse_fields, _parse_gatefilter, _set_limits
from ._common_dealias import _parse_rays_wrap_around, _parse_nyquist_vel
from ._fast_edge_finder import _fast_edge_finder


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
# * Improve perfornace by implementing a priority queue in the _EdgeTracker
#   object. See comments in the class for details.


def dealias_region_based(
        radar, ref_vel_field=None, interval_splits=3, interval_limits=None,
        skip_between_rays=100, skip_along_ray=100, centered=True,
        nyquist_vel=None, check_nyquist_uniform=True, gatefilter=False,
        rays_wrap_around=None, keep_original=False, set_limits=True,
        vel_field=None, corr_vel_field=None, **kwargs):
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
    ref_vel_field : str or None, optional
         Field in radar containing a reference velocity field used to anchor
         the unfolded velocities once the algorithm completes. Typically this
         field is created by simulating the radial velocities from wind data
         from an atmospheric sonding using
         :py:func:`pyart.util.simulated_vel_from_profile`.
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
    nyquist_velocity : array like or float, optional
        Nyquist velocity in unit identical to those stored in the radar's
        velocity field, either for each sweep or a single value which will be
        used for all sweeps.  None will attempt to determine this value from
        the Radar object.
    check_nyquist_uniform : bool, optional
        True to check if the Nyquist velocities are uniform for all rays
        within a sweep, False will skip this check.  This parameter is ignored
        when the nyquist_velocity parameter is not None.
    gatefilter : GateFilter, None or False, optional.
        A GateFilter instance which specified which gates should be
        ignored when performing de-aliasing.  A value of None created this
        filter from the radar moments using any additional arguments by
        passing them to :py:func:`moment_based_gate_filter`. False, the
        default, disables filtering including all gates in the dealiasing.
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
    set_limits : bool, optional
        True to set valid_min and valid_max elements in the returned
        dictionary.  False will not set these dictionary elements.
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
    nyquist_vel = _parse_nyquist_vel(nyquist_vel, radar, check_nyquist_uniform)

    # parse ref_vel_field
    if ref_vel_field is None:
        ref_vdata = None
    else:
        ref_vdata = radar.fields[ref_vel_field]['data']

    # exclude masked and invalid velocity gates
    gatefilter.exclude_masked(vel_field)
    gatefilter.exclude_invalid(vel_field)
    gfilter = gatefilter.gate_excluded

    # perform dealiasing
    vdata = radar.fields[vel_field]['data'].view(np.ndarray)
    data = vdata.copy()     # dealiased velocities

    # loop over sweeps
    for nsweep, sweep_slice in enumerate(radar.iter_slice()):

        # extract sweep data
        sdata = vdata[sweep_slice].copy()   # is a copy needed here?
        scorr = data[sweep_slice]
        sfilter = gfilter[sweep_slice]

        # find nyquist velocity and interval segmentation limits
        nyquist_interval = nyquist_vel[nsweep] * 2.
        if interval_limits is None:
            nvel = nyquist_vel[nsweep]
            valid_sdata = sdata[~sfilter]
            s_interval_limits = _find_sweep_interval_splits(
                nvel, interval_splits, valid_sdata, nsweep)
        else:
            s_interval_limits = interval_limits

        # find regions in original data
        labels, nfeatures = _find_regions(sdata, sfilter, s_interval_limits)
        # skip sweep if all gates are masked or only a single region
        if nfeatures < 2:
            continue
        bincount = np.bincount(labels.ravel())
        num_masked_gates = bincount[0]
        region_sizes = bincount[1:]

        # find all edges between regions
        indices, edge_count, velos = _edge_sum_and_count(
            labels, num_masked_gates, sdata, rays_wrap_around,
            skip_between_rays, skip_along_ray)

        # no unfolding required if no edges exist between regions
        if len(edge_count) == 0:
            continue

        # find the number of folds in the regions
        region_tracker = _RegionTracker(region_sizes)
        edge_tracker = _EdgeTracker(indices, edge_count, velos,
                                    nyquist_interval, nfeatures+1)
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
        nwrap = np.take(region_tracker.unwrap_number, labels)
        scorr += nwrap * nyquist_interval

        # anchor unfolded velocities against reference velocity
        if ref_vdata is not None:
            sref = ref_vdata[sweep_slice]
            gfold = (sref-scorr).mean()/nyquist_interval
            gfold = round(gfold)

            # Anchor specific regions against reference velocity
            # Do this by constraining cost function due to difference
            # from reference velocity and to 2D continuity
            new_interval_limits = np.linspace(scorr.min(), scorr.max(), 10)
            labels_corr, nfeatures_corr = _find_regions(
                scorr, sfilter, new_interval_limits)

            bounds_list = [(x, y) for (x, y) in zip(-6*np.ones(nfeatures_corr),
                           5*np.ones(nfeatures_corr))]
            scorr_means = np.zeros(nfeatures_corr)
            sref_means = np.zeros(nfeatures_corr)
            for reg in range(1, nfeatures_corr+1):
                scorr_means[reg-1] = np.ma.mean(scorr[labels_corr == reg])
                sref_means[reg-1] = np.ma.mean(sref[labels_corr == reg])

            def cost_function(x):
                return _cost_function(x, scorr_means, sref_means,
                                      nyquist_interval, nfeatures_corr)

            def gradient(x):
                return _gradient(x, scorr_means, sref_means,
                                 nyquist_interval, nfeatures_corr)

            nyq_adjustments = fmin_l_bfgs_b(
                cost_function, gfold*np.ones((nfeatures_corr)), disp=True,
                fprime=gradient, bounds=bounds_list, maxiter=200,
                )

            i = 0
            for reg in range(1, nfeatures_corr):
                scorr[labels == reg] += (nyquist_interval *
                                         np.round(nyq_adjustments[0][i]))
                i = i + 1

    # fill_value from the velocity dictionary if present
    fill_value = radar.fields[vel_field].get('_FillValue', get_fillvalue())

    # mask filtered gates
    if np.any(gfilter):
        data = np.ma.array(data, mask=gfilter, fill_value=fill_value)

    # restore original values where dealiasing not applied
    if keep_original:
        data[gfilter] = radar.fields[vel_field]['data'][gfilter]

    # return field dictionary containing dealiased Doppler velocities
    corr_vel = get_metadata(corr_vel_field)
    corr_vel['data'] = data
    corr_vel['_FillValue'] = fill_value

    if set_limits:
        # set valid_min and valid_max in corr_vel
        _set_limits(data, nyquist_vel, corr_vel)

    return corr_vel


def _find_sweep_interval_splits(nyquist, interval_splits, velocities, nsweep):
    """ Return the interval limits for a given sweep. """
    # The Nyquist interval is split into interval_splits  equal sized areas.
    # If velocities outside the Nyquist are present the number and
    # limits of the interval splits must be adjusted so that theses
    # velocities are included in one of the splits.
    add_start = add_end = 0
    interval = (2. * nyquist) / (interval_splits)
    # no change from default if all gates filtered
    if len(velocities) != 0:
        max_vel = velocities.max()
        min_vel = velocities.min()
        if max_vel > nyquist or min_vel < -nyquist:
            msg = ("Velocities outside of the Nyquist interval found in "
                   "sweep %i." % (nsweep))
            warnings.warn(msg, UserWarning)
            # additional intervals must be added to capture the velocities
            # outside the nyquist limits
            add_start = int(np.ceil((max_vel - nyquist) / (interval)))
            add_end = int(np.ceil(-(min_vel + nyquist) / (interval)))

    start = -nyquist - add_start * interval
    end = nyquist + add_end * interval
    num = interval_splits + 1 + add_start + add_end
    return np.linspace(start, end, num, endpoint=True)


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


def _edge_sum_and_count(labels, num_masked_gates, data,
                        rays_wrap_around, max_gap_x, max_gap_y):
    """
    Find all edges between labels regions.

    Returns the indices, count and velocities of all edges.
    """
    total_nodes = labels.shape[0] * labels.shape[1] - num_masked_gates
    if rays_wrap_around:
        total_nodes += labels.shape[0] * 2

    indices, velocities = _fast_edge_finder(
        labels.astype('int32'), data.astype('float32'),
        rays_wrap_around, max_gap_x, max_gap_y, total_nodes)
    index1, index2 = indices
    vel1, vel2 = velocities
    count = np.ones_like(vel1, dtype=np.int32)

    # return early if not edges were found
    if len(vel1) == 0:
        return ([], []), [], ([], [])

    # find the unique edges, procedure based on method in
    # scipy.sparse.coo_matrix.sum_duplicates
    # except we have three data arrays, vel1, vel2, and count
    order = np.lexsort((index1, index2))
    index1 = index1[order]
    index2 = index2[order]
    vel1 = vel1[order]
    vel2 = vel2[order]
    count = count[order]

    unique_mask = ((index1[1:] != index1[:-1]) |
                   (index2[1:] != index2[:-1]))
    unique_mask = np.append(True, unique_mask)
    index1 = index1[unique_mask]
    index2 = index2[unique_mask]

    unique_inds, = np.nonzero(unique_mask)
    vel1 = np.add.reduceat(vel1, unique_inds, dtype=vel1.dtype)
    vel2 = np.add.reduceat(vel2, unique_inds, dtype=vel2.dtype)
    count = np.add.reduceat(count, unique_inds, dtype=count.dtype)

    return (index1, index2), count, (vel1, vel2)


def _combine_regions(region_tracker, edge_tracker):
    """ Returns True when done. """
    # Edge parameters from edge with largest weight
    status, extra = edge_tracker.pop_edge()
    if status:
        return True
    node1, node2, weight, diff, edge_number = extra
    rdiff = int(np.round(diff))

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


# Minimize cost function that is sum of difference between regions and
# sounding
def _cost_function(nyq_vector, vels_slice_means,
                   svels_slice_means, v_nyq_vel, nfeatures):
    """ Cost function for minimization in region based algorithm """
    cost = 0
    i = 0

    for reg in range(nfeatures):
        add_value = 0
        # Deviance from sounding
        add_value = 1*(vels_slice_means[reg] +
                       np.round(nyq_vector[i])*v_nyq_vel -
                       svels_slice_means[reg])**2

        # Region continuity
        vels_without_cur = np.delete(vels_slice_means, reg)
        diffs = np.square(vels_slice_means[reg]-vels_without_cur)
        if(len(diffs) > 0):
            the_min = np.argmin(diffs)
            vel_wo_cur = vels_without_cur[the_min]
        else:
            vel_wo_cur = vels_slice_means[reg]
        add_value2 = np.square(vels_slice_means[reg] +
                               np.round(nyq_vector[i])*v_nyq_vel -
                               vel_wo_cur)
        if np.isfinite(add_value2):
            add_value += .1*add_value2

        if np.isfinite(add_value):
            cost += add_value
        i = i + 1

    return cost


def _gradient(nyq_vector, vels_slice_means, svels_slice_means,
              v_nyq_vel, nfeatures):
    """ Gradient of cost function for minimization
        in region based algorithm """
    gradient_vector = np.zeros(len(nyq_vector))
    i = 0
    for reg in range(nfeatures):
        add_value = (vels_slice_means[reg] +
                     np.round(nyq_vector[i])*v_nyq_vel -
                     svels_slice_means[reg])
        if np.isfinite(add_value):
            gradient_vector[i] = 2*add_value*v_nyq_vel

        # Regional continuity
        vels_without_cur = np.delete(vels_slice_means, reg)
        diffs = np.square(vels_slice_means[reg]-vels_without_cur)
        if(len(diffs) > 0):
            the_min = np.argmin(diffs)
            vel_wo_cur = vels_without_cur[the_min]
        else:
            vel_wo_cur = vels_slice_means[reg]

        add_value2 = (vels_slice_means[reg] +
                      np.round(nyq_vector[i])*v_nyq_vel -
                      vel_wo_cur)

        if np.isfinite(add_value2):
            gradient_vector[i] += 2*.1*add_value2*v_nyq_vel

        i = i + 1

    return gradient_vector


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

    def __init__(self, indices, edge_count, velocities, nyquist_interval,
                 nnodes):
        """ initialize """

        nedges = int(len(indices[0]) / 2)

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

        # fill out data from the provides indicies, edge counts and velocities
        edge = 0
        idx1, idx2 = indices
        vel1, vel2 = velocities
        for i, j, count, vel, nvel in zip(idx1, idx2, edge_count, vel1, vel2):
            if i < j:
                continue
            self.node_alpha[edge] = i
            self.node_beta[edge] = j
            self.sum_diff[edge] = ((vel - nvel) / nyquist_interval)
            self.weight[edge] = count
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
