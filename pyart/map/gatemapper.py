import numpy as np

from copy import deepcopy
from scipy.spatial import KDTree
from ..core import Radar, geographic_to_cartesian

class GateMapper(object):
    """

    The GateMapper object will, given one radar's gate, find the gate in another radar's volume
    that is closest in location to the specified gate. GateMapper will generate this mapping function
    that provides the sweep and ray index in the other radar that is closest in physical location to
    the first radar's gate.

    """
    def __init__(self, original_radar: Radar, new_radar: Radar, tol=500.):
        gate_x = new_radar.gate_x['data']
        gate_y = new_radar.gate_y['data']
        gate_z = new_radar.gate_altitude['data']
        proj_dict = {'proj': 'pyart_aeqd', 'lon_0': new_radar.longitude["data"],
                     'lat_0': new_radar.latitude["data"]}
        self.original_radar_x, self.original_radar_y = geographic_to_cartesian(
            new_radar.longitude["data"], new_radar.latitude["data"],
            proj_dict)
        data = np.stack([gate_x.flatten(), gate_y.flatten(), gate_z.flatten()], axis=1)
        self._kdtree = KDTree(data)
        self.tolerance = tol
        self.new_radar = new_radar
        self.original_radar = original_radar
        self._index_map = np.stack([np.nan * np.ones_like(self.original_radar.gate_x["data"]),
                                   np.nan * np.ones_like(self.original_radar.gate_x["data"])],
                                  axis=2)

        new_points = np.stack([self.original_radar.gate_x["data"].flatten(),
                               self.original_radar.gate_y["data"].flatten(),
                               self.original_radar.gate_altitude["data"].flatten()], axis=1)
        dists, inds = self._kdtree.query(new_points)
        inds = np.where(np.abs(dists) < self.tolerance, inds[:], -32767).astype(int)
        inds = np.reshape(inds, self.original_radar.gate_x["data"].shape)
        self._index_map[:, :, 0] = (inds / self.new_radar.gate_x["data"].shape[1]).astype(int)
        self._index_map[:, :, 1] = (inds - self.new_radar.gate_x["data"].shape[1] *
                                   (inds / self.new_radar.gate_x["data"].shape[1]).astype(int))

    def __getitem__(self, key: tuple):
        return (int(self._index_map[key[0], key[1], 0]), int(self._index_map[key[0], key[1], 1]))

    def mapped_radar(self, field_list):
        mapped_radar = deepcopy(self.new_radar)
        for field in field_list:
            mapped_radar.fields[field]['data'] = np.ma.masked_where(True, mapped_radar.fields[field]['data'])

        for i in range(self.original_radar.nrays):
            for j in range(self.original_radar.ngates):
                index = self[i, j]
                if index[0] > 0:
                    for field in field_list:
                        mapped_radar.fields[field]['data'][index] = self.original_radar.fields[field]['data'][i, j]

        return mapped_radar