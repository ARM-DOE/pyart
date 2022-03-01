"""
Utilities for finding gates with equivalent locations between radars for easy
comparison.

"""

from copy import deepcopy
import numpy as np

from scipy.spatial import KDTree

from ..core import Radar, geographic_to_cartesian


class GateMapper(object):
    """
    The GateMapper class will, given one radar's gate, find the gate in another
    radar's volume that is closest in location to the specified gate.
    GateMapper will use a kd-tree in order to generate a mapping function that
    provides the sweep and ray index in the other radar that is closest in
    physical location to the first radar's gate. This functionality provides
    easy mapping of equivalent locations between radar objects with simple
    indexing. In addition, returning a mapped radar object is also supported.

    Examples
    --------
    >>> grid_mapper = pyart.map.GridMapper(src, dest)
    >>> # Get the destination radar's equivalent of (2, 2) in the source radar's
    >>> # coordinates
    >>> dest_index = grid_mapper[2, 2]
    >>> radar_mapped = grid_mapper.mapped_radar(['reflectivity'])

    Parameters
    ----------
    src_radar : pyart.core.Radar
        The source radar data.
    dest_radar : pyart.core.Radar
        The destination radar to map the source radar data onto.
    tol : float
        The difference in meters between the source and destination gate allowed
        for an adequate match.

    """
    def __init__(self, src_radar: Radar, dest_radar: Radar, tol=500.):
        gate_x = dest_radar.gate_x['data']
        gate_y = dest_radar.gate_y['data']
        gate_z = dest_radar.gate_altitude['data']
        proj_dict = {'proj': 'pyart_aeqd', 'lon_0': dest_radar.longitude["data"],
                     'lat_0': dest_radar.latitude["data"]}
        self.src_radar_x, self.src_radar_y = geographic_to_cartesian(
            src_radar.longitude["data"], src_radar.latitude["data"],
            proj_dict)
        data = np.stack(
            [gate_x.flatten(), gate_y.flatten(), gate_z.flatten()], axis=1)
        self._kdtree = KDTree(data)
        self.tolerance = tol
        self.dest_radar = dest_radar
        self.src_radar = src_radar
        self._index_map = np.stack(
            [np.nan * np.ones_like(
                self.src_radar.gate_x["data"]),
                np.nan * np.ones_like(self.src_radar.gate_x["data"])],
                axis=2)

        new_points = np.stack(
            [self.src_radar.gate_x["data"].flatten() + self.src_radar_x,
             self.src_radar.gate_y["data"].flatten() + self.src_radar_y,
             self.src_radar.gate_altitude["data"].flatten()], axis=1)
        dists, inds = self._kdtree.query(new_points, distance_upper_bound=tol)
        inds = np.where(
            np.abs(dists) < self.tolerance, inds[:], -32767).astype(int)
        inds = np.reshape(inds, self.src_radar.gate_x["data"].shape)
        self._index_map[:, :, 0] = (
            inds / self.dest_radar.gate_x["data"].shape[1]).astype(int)
        self._index_map[:, :, 1] = (
            inds - self.dest_radar.gate_x["data"].shape[1] *
            (inds / self.dest_radar.gate_x["data"].shape[1]).astype(int))

    def __getitem__(self, key: tuple):
        """ Return the equivalent index in the destination radar's
        coordinates. """
        x = (int(self._index_map[key[0], key[1], 0]),
                 int(self._index_map[key[0], key[1], 1]))
        if x[0] > 0:
            return x
        else:
            return None, None

    def mapped_radar(self, field_list):
        """
        This returns a version of the destination radar with the fields in
        field_list from the source radar mapped into the destination radar's
        coordinate system.

        Parameters
        ----------
        field_list: list of str or str
            The list of fields to map.

        Returns
        -------
        mapped_radar:
            The destination radar with the fields from the source radar mapped
            into the destination radar's coordinate system.
        """
        mapped_radar = deepcopy(self.dest_radar)
        if isinstance(field_list, str):
            field_list = [field_list]

        for field in field_list:
            mapped_radar.fields[field]['data'] = np.ma.masked_where(
                True, mapped_radar.fields[field]['data'])

        for i in range(self.src_radar.nrays):
            for j in range(self.src_radar.ngates):
                index = self[i, j]
                if index[0] is not None:
                    for field in field_list:
                        mapped_radar.fields[
                            field]['data'][index] = self.src_radar.fields[
                                field]['data'][i, j]

        return mapped_radar
