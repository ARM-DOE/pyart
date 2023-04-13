"""
Utilities for finding gates with equivalent locations between radars for easy comparison.

"""

from copy import deepcopy

import numpy as np
from scipy.spatial import KDTree

import pyart.filters

from ..core import Radar, geographic_to_cartesian


class GateMapper:
    """
    The GateMapper class will, given one radar's gate, find the gate in another radar's volume
    that is closest in location to the specified gate. GateMapper will use a kd-tree in order to generate a
    mapping function that provides the sweep and ray index in the other radar that is closest in physical
    location to the first radar's gate. This functionality provides easy mapping of equivalent locations between
    radar objects with simple indexing. In addition, returning a mapped radar object is also supported.

    Attributes
    ----------
    src_radar_x, src_radar_y: float
        The source radar's x and y location in the the destination radar's Cartesian coordinates.
    distance_tolerance: float
        The distance tolerance in meters for each gate in meters.
    time_tolerance: float
        The time tolerance in meters for each gate in seconds.
    gatefilter_src: pyart.filters.GateFilter
        The gatefilter to apply to the source radar data when mapping to the destination.
    src_radar: pyart.core.Radar
        The source radar data.
    dest_radar: pyart.core.Radar
        The destination radar to map the source radar data onto.
    src_radar_time: float np.array
        The array of times of each gate in the source radar.
    dest_radar_time: float np.array
        The array of times of each gate in the destination radar.

    Examples
    --------
    >>> gate_mapper = pyart.map.GateMapper(src, dest)
    >>> # Get the destination radar's equivalent of (2, 2) in the source radar's coordinates
    >>> dest_index = gate_mapper[2, 2]
    >>> radar_mapped = gate_mapper.mapped_radar(['reflectivity'])

    Parameters
    ----------
    src_radar: pyart.core.Radar
        The source radar data.
    dest_radar: pyart.core.Radar
        The destination radar to map the source radar data onto.
    gatefilter_src: pyart.filters.GateFilter, or None
        The gatefilter to apply to the source radar data before mapping
    distance_tolerance: float
        The difference in meters between the source and destination gate allowed for an adequate match.
    time_tolerance: float
        The difference in time between the source and destination radar rays.
    """

    def __init__(
        self,
        src_radar: Radar,
        dest_radar: Radar,
        distance_tolerance=500.0,
        gatefilter_src=None,
        time_tolerance=60.0,
    ):
        gate_x = dest_radar.gate_x["data"]
        gate_y = dest_radar.gate_y["data"]
        gate_z = dest_radar.gate_altitude["data"]
        proj_dict = {
            "proj": "pyart_aeqd",
            "lon_0": dest_radar.longitude["data"],
            "lat_0": dest_radar.latitude["data"],
        }
        self.src_radar_x, self.src_radar_y = geographic_to_cartesian(
            src_radar.longitude["data"], src_radar.latitude["data"], proj_dict
        )
        data = np.stack([gate_x.flatten(), gate_y.flatten(), gate_z.flatten()], axis=1)
        self._kdtree = KDTree(data)
        self.dest_radar = dest_radar
        self.src_radar = src_radar
        self.src_radar_time = np.tile(
            src_radar.time["data"], (1, src_radar.ngates)
        ).flatten()
        self.dest_radar_time = np.tile(
            dest_radar.time["data"], (1, dest_radar.ngates)
        ).flatten()
        if gatefilter_src is None:
            self.gatefilter_src = pyart.filters.GateFilter(self.src_radar)
            self.gatefilter_src.exclude_none()
        else:
            self.gatefilter_src = gatefilter_src
        self._time_tolerance = time_tolerance
        self._distance_tolerance = distance_tolerance
        self.distance_tolerance = distance_tolerance

    def __getitem__(self, key: tuple):
        """Return the equivalent index in the destination radar's coordinates"""
        x = (
            int(self._index_map[key[0], key[1], 0]),
            int(self._index_map[key[0], key[1], 1]),
        )
        if x[0] > 0:
            return x
        else:
            return None, None

    @property
    def distance_tolerance(self):
        """
        Getter for distance_tolerance property of GateMapper class.

        Returns
        -------
        tolerance: float
            The current distance tolerance of the GateMapper in meters.
        """
        return self._distance_tolerance

    @property
    def time_tolerance(self):
        """
        Getter for time_tolerance property of GateMapper class.

        Returns
        -------
        time_tolerance: float
            The current time tolerance of the GateMapper in meters.
        """
        return self._time_tolerance

    @time_tolerance.setter
    def time_tolerance(self, time_tolerance):
        """
        Setter for the time_tolerance property of GateMapper class. This will
        also reset the mapping function inside GateMapper so that no additional
        calls are needed to reset the tolerance.

        Parameters
        ----------
        time_tolerance: float
            The new time tolerance of the GateMapper in seconds.
        """
        self._time_tolerance = time_tolerance
        self.tolerance = self._distance_tolerance

    @distance_tolerance.setter
    def distance_tolerance(self, distance_tolerance):
        """
        Setter for the distance_tolerance property of GateMapper class. This will
        also reset the mapping function inside GateMapper so that no additional
        calls are needed to reset the tolerance.

        Parameters
        ----------
        distance_tolerance: float
            The new distance tolerance of the GateMapper in meters.
        """
        self._index_map = np.stack(
            [
                np.nan * np.ones_like(self.src_radar.gate_x["data"]),
                np.nan * np.ones_like(self.src_radar.gate_x["data"]),
            ],
            axis=2,
        )
        new_points = np.stack(
            [
                self.src_radar.gate_x["data"].flatten() + self.src_radar_x,
                self.src_radar.gate_y["data"].flatten() + self.src_radar_y,
                self.src_radar.gate_altitude["data"].flatten(),
            ],
            axis=1,
        )
        dists, inds = self._kdtree.query(
            new_points, distance_upper_bound=distance_tolerance
        )
        inds[inds == len(self.dest_radar_time)] = (
            inds[inds == len(self.dest_radar_time)] - 1
        )
        times = np.abs(self.src_radar_time - self.dest_radar_time[inds])
        inds = np.where(
            np.logical_and(
                times < self._time_tolerance, np.abs(dists) < distance_tolerance
            ),
            inds[:],
            -32767,
        ).astype(int)
        inds = np.reshape(inds, self.src_radar.gate_x["data"].shape)

        self._index_map[:, :, 0] = (
            inds / self.dest_radar.gate_x["data"].shape[1]
        ).astype(int)
        self._index_map[:, :, 1] = inds - self.dest_radar.gate_x["data"].shape[1] * (
            inds / self.dest_radar.gate_x["data"].shape[1]
        ).astype(int)
        self._distance_tolerance = distance_tolerance

    def mapped_radar(self, field_list):
        """
        This returns a version of the destination radar with the fields in field_list from the source radar
        mapped into the destination radar's coordinate system.

        Parameters
        ----------
        field_list: list of str or str
            The list of fields to map.

        Returns
        -------
        mapped_radar:
            The destination radar with the fields from the source radar mapped into the destination radar's
            coordinate system.
        """
        mapped_radar = deepcopy(self.dest_radar)
        if isinstance(field_list, str):
            field_list = [field_list]

        src_fields = {}
        for field in field_list:
            if field in list(mapped_radar.fields.keys()):
                mapped_radar.fields[field]["data"] = np.ma.masked_where(
                    True, mapped_radar.fields[field]["data"]
                )
            else:
                mapped_radar.fields[field] = deepcopy(self.src_radar.fields[field])
                mapped_radar.fields[field]["data"] = np.ma.masked_where(
                    True, np.ma.zeros((mapped_radar.nrays, mapped_radar.ngates))
                )
            src_fields[field] = np.ma.masked_where(
                self.gatefilter_src.gate_excluded, self.src_radar.fields[field]["data"]
            )

        for i in range(self.src_radar.nrays):
            for j in range(self.src_radar.ngates):
                index = self[i, j]
                if index[0] is not None:
                    for field in field_list:
                        mapped_radar.fields[field]["data"][index] = src_fields[field][
                            i, j
                        ]
        del src_fields
        return mapped_radar
