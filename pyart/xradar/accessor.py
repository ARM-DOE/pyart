"""
Utilities for interfacing between xradar and Py-ART

"""


import copy

import numpy as np
import pandas as pd
from datatree import DataTree, formatting, formatting_html
from datatree.treenode import NodePath
from xarray import DataArray, Dataset, concat
from xarray.core import utils
from xradar.util import get_sweep_keys

from ..config import get_metadata
from ..core.transforms import (
    antenna_vectors_to_cartesian,
    cartesian_to_geographic,
    cartesian_vectors_to_geographic,
)
from ..lazydict import LazyLoadDict


class Xgrid:
    def __init__(self, grid_ds):
        """
        Wraps a Cf-compliant xarray Dataset into a PyART Grid Object.
        Note that the times must not be decoded by xr.open_dataset when loading the file.

        Parameters
        ----------
        grid_ds: xarray Dataset
            The xarray Dataset to convert to a Py-ART grid.
        """
        if "units" not in list(grid_ds["time"].attrs.keys()):
            raise RuntimeError(
                "decode_times must be set to false when opening grid file!"
            )
        self.ds = grid_ds
        self.time = dict(data=np.atleast_1d(self.ds["time"].values))
        self.time.update(self.ds["time"].attrs)
        self.fields = {}
        self._find_fields()
        self.origin_altitude = dict(
            data=np.atleast_1d(self.ds["origin_altitude"].values)
        )
        self.origin_altitude.update(self.ds["origin_altitude"].attrs)
        self.origin_latitude = dict(
            data=np.atleast_1d(self.ds["origin_latitude"].values)
        )
        self.origin_latitude.update(self.ds["origin_latitude"].attrs)
        self.origin_longitude = dict(
            data=np.atleast_1d(self.ds["origin_longitude"].values)
        )
        self.origin_longitude.update(self.ds["origin_longitude"].attrs)
        self.z = dict(data=np.atleast_1d(self.ds["z"].values))
        self.z.update(self.ds["z"].attrs)
        self.y = dict(data=np.atleast_1d(self.ds["y"].values))
        self.y.update(self.ds["y"].attrs)
        self.x = dict(data=np.atleast_1d(self.ds["x"].values))
        self.x.update(self.ds["x"].attrs)
        self.nradar = len(self.ds["nradar"].values)
        self.radar_altitude = dict(data=np.atleast_1d(self.ds["radar_altitude"].values))
        self.radar_altitude.update(self.ds["radar_altitude"].attrs)
        self.radar_longitude = dict(
            data=np.atleast_1d(self.ds["radar_longitude"].values)
        )
        self.radar_longitude.update(self.ds["radar_longitude"].attrs)
        self.radar_latitude = dict(data=np.atleast_1d(self.ds["radar_latitude"].values))
        self.radar_latitude.update(self.ds["radar_latitude"].attrs)
        self.radar_time = dict(data=np.atleast_1d(self.ds["radar_time"].values))
        self.radar_time.update(self.ds["radar_time"].attrs)
        self.radar_name = dict(data=self.ds["radar_name"].values.astype("<U12"))
        self.radar_name.update(self.ds["radar_name"].attrs)
        self.projection = self.ds["projection"].attrs
        if "_include_lon_0_lat_0" in list(self.projection.keys()):
            if self.projection["_include_lon_0_lat_0"].lower() == "true":
                self.projection["_include_lon_0_lat_0"] = True
            else:
                self.projection["_include_lon_0_lat_0"] = False

        self.init_point_altitude()
        self.init_point_longitude_latitude()
        self.init_point_x_y_z()

    def _find_fields(self):
        for key in list(self.ds.variables.keys()):
            if self.ds[key].dims == ("time", "z", "y", "x"):
                self.fields[key] = {}
                self.fields[key]["data"] = self.ds[key].values.squeeze()
                self.fields[key].update(self.ds[key].attrs)

    def get_projparams(self):
        projparams = self.projection.copy()
        if projparams.pop("_include_lon_0_lat_0", False):
            projparams["lon_0"] = self.origin_longitude["data"][0]
            projparams["lat_0"] = self.origin_latitude["data"][0]
        return projparams

    @property
    def metadata(self):
        return self.ds.attrs

    @property
    def ny(self):
        return self.ds.sizes["y"]

    @property
    def nx(self):
        return self.ds.sizes["x"]

    @property
    def nz(self):
        return self.ds.sizes["z"]

    # Attribute init/reset methods
    def init_point_x_y_z(self):
        """Initialize or reset the point_{x, y, z} attributes."""
        self.point_x = LazyLoadDict(get_metadata("point_x"))
        self.point_x.set_lazy("data", _point_data_factory(self, "x"))

        self.point_y = LazyLoadDict(get_metadata("point_y"))
        self.point_y.set_lazy("data", _point_data_factory(self, "y"))

        self.point_z = LazyLoadDict(get_metadata("point_z"))
        self.point_z.set_lazy("data", _point_data_factory(self, "z"))

    def init_point_longitude_latitude(self):
        """
        Initialize or reset the point_{longitude, latitudes} attributes.
        """
        point_longitude = LazyLoadDict(get_metadata("point_longitude"))
        point_longitude.set_lazy("data", _point_lon_lat_data_factory(self, 0))
        self.point_longitude = point_longitude

        point_latitude = LazyLoadDict(get_metadata("point_latitude"))
        point_latitude.set_lazy("data", _point_lon_lat_data_factory(self, 1))
        self.point_latitude = point_latitude

    def init_point_altitude(self):
        """Initialize the point_altitude attribute."""
        point_altitude = LazyLoadDict(get_metadata("point_altitude"))
        point_altitude.set_lazy("data", _point_altitude_data_factory(self))
        self.point_altitude = point_altitude

    def get_point_longitude_latitude(self, level=0, edges=False):
        """
        Return arrays of longitude and latitude for a given grid height level.

        Parameters
        ----------
        level : int, optional
            Grid height level at which to determine latitudes and longitudes.
            This is not currently used as all height level have the same
            layout.
        edges : bool, optional
            True to calculate the latitude and longitudes of the edges by
            interpolating between Cartesian coordinates points and
            extrapolating at the boundaries. False to calculate the locations
            at the centers.

        Returns
        -------
        longitude, latitude : 2D array
            Arrays containing the latitude and longitudes, in degrees, of the
            grid points or edges between grid points for the given height.

        """
        x = self.x["data"]
        y = self.y["data"]
        projparams = self.get_projparams()
        return cartesian_vectors_to_geographic(x, y, projparams, edges=edges)

    def add_field(self, field_name, dic, replace_existing=False):
        """
        Add a field to the object.

        Parameters
        ----------
        field_name : str
            Name of the field to add to the dictionary of fields.
        dic : dict
            Dictionary contain field data and metadata.
        replace_existing : bool, optional
            True to replace the existing field with key field_name if it
            exists, loosing any existing data. False will raise a ValueError
            when the field already exists.

        """
        # check that the field dictionary to add is valid
        if field_name in self.fields and replace_existing is False:
            err = f"A field with name: {field_name} already exists"
            raise ValueError(err)
        if "data" not in dic:
            raise KeyError("dic must contain a 'data' key")
        if dic["data"].shape != (self.nz, self.ny, self.nx):
            t = (self.nz, self.ny, self.nx)
            err = "'data' has invalid shape, should be (%i, %i)" % t
            raise ValueError(err)
        self.fields[field_name] = dic

    def to_xarray(self):
        """
        Convert the Grid object to an xarray format.

        Attributes
        ----------
        time : dict
            Time of the grid.
        fields : dict of dicts
            Moments from radars or other variables.
        longitude, latitude : dict, 2D
            Arrays of latitude and longitude for the grid height level.
        x, y, z : dict, 1D
            Distance from the grid origin for each Cartesian coordinate axis
            in a one dimensional array.

        """

        lon, lat = self.get_point_longitude_latitude()
        z = self.z["data"]
        y = self.y["data"]
        x = self.x["data"]

        time = self.ds.time.values

        ds = Dataset()
        for field in list(self.fields.keys()):
            field_data = self.fields[field]["data"]
            data = DataArray(
                np.expand_dims(field_data, 0),
                dims=("time", "z", "y", "x"),
                coords={
                    "time": (["time"], time),
                    "z": (["z"], z),
                    "lat": (["y", "x"], lat),
                    "lon": (["y", "x"], lon),
                    "y": (["y"], y),
                    "x": (["x"], x),
                },
            )
            for meta in list(self.fields[field].keys()):
                if meta != "data":
                    data.attrs.update({meta: self.fields[field][meta]})

            ds[field] = data
            ds.lon.attrs = [
                ("long_name", "longitude of grid cell center"),
                ("units", "degree_E"),
                ("standard_name", "Longitude"),
            ]
            ds.lat.attrs = [
                ("long_name", "latitude of grid cell center"),
                ("units", "degree_N"),
                ("standard_name", "Latitude"),
            ]

            ds.z.attrs = get_metadata("z")
            ds.y.attrs = get_metadata("y")
            ds.x.attrs = get_metadata("x")

            ds.z.encoding["_FillValue"] = None
            ds.lat.encoding["_FillValue"] = None
            ds.lon.encoding["_FillValue"] = None
            ds.close()
        return ds


class Xradar:
    def __init__(self, xradar, default_sweep="sweep_0", scan_type=None):
        self.xradar = xradar
        self.scan_type = scan_type or "ppi"
        self.sweep_group_names = get_sweep_keys(self.xradar)
        self.nsweeps = len(self.sweep_group_names)
        self.combined_sweeps = self._combine_sweeps()
        self.fields = self._find_fields(self.combined_sweeps)
        self.time = dict(
            data=(self.combined_sweeps.time - self.combined_sweeps.time.min()).astype(
                "int64"
            )
            / 1e9,
            units=f"seconds since {pd.to_datetime(self.combined_sweeps.time.min().values).strftime('%Y-%m-%d %H:%M:%S.0')}",
            calendar="gregorian",
        )
        self.range = dict(data=self.combined_sweeps.range.values)
        self.azimuth = dict(data=self.combined_sweeps.azimuth.values)
        self.elevation = dict(data=self.combined_sweeps.elevation.values)
        self.fixed_angle = dict(data=self.combined_sweeps.sweep_fixed_angle.values)
        self.antenna_transition = None
        self.latitude = dict(
            data=np.expand_dims(self.xradar["latitude"].values, axis=0)
        )
        self.longitude = dict(
            data=np.expand_dims(self.xradar["longitude"].values, axis=0)
        )
        self.altitude = dict(
            data=np.expand_dims(self.xradar["altitude"].values, axis=0)
        )
        self.sweep_end_ray_index = dict(
            data=self.combined_sweeps.ngates.groupby("sweep_number").max().values
        )
        self.sweep_start_ray_index = dict(
            data=self.combined_sweeps.ngates.groupby("sweep_number").min().values
        )
        self.metadata = dict(**self.xradar.attrs)
        self.ngates = len(self.range["data"])
        self.nrays = len(self.azimuth["data"])
        self.instrument_parameters = self.find_instrument_parameters()
        self.init_gate_x_y_z()
        self.init_gate_alt()

    def __repr__(self):
        return formatting.datatree_repr(self.xradar)

    def _repr_html_(self):
        return formatting_html.datatree_repr(self.xradar)

    def __getitem__(self: DataTree, key):
        """
        Access child nodes, variables, or coordinates stored anywhere in this tree.

        Returned object will be either a DataTree or DataArray object depending on whether the key given points to a
        child or variable.

        Parameters
        ----------
        key : str
            Name of variable / child within this node, or unix-like path to variable / child within another node.

        Returns
        -------
        Union[DataTree, DataArray]
        """

        # Either:
        if utils.is_dict_like(key):
            # dict-like indexing
            raise NotImplementedError("Should this index over whole tree?")
        elif isinstance(key, str):
            # path-like: a name of a node/variable, or path to a node/variable
            path = NodePath(key)
            return self.xradar._get_item(path)
        elif utils.is_list_like(key):
            # iterable of variable names
            raise NotImplementedError(
                "Selecting via tags is deprecated, and selecting multiple items should be "
                "implemented via .subset"
            )
        else:
            raise ValueError(f"Invalid format for key: {key}")

        # Iterators

    def find_instrument_parameters(self):
        # By default, check the radar_parameters first
        if "radar_parameters" in list(self.xradar.children):
            radar_param_dict = self.xradar["radar_parameters"].ds.to_dict(data="array")
            instrument_parameters = radar_param_dict["data_vars"]
            instrument_parameters.update(radar_param_dict["attrs"])

        else:
            instrument_parameters = {}

        # Check to see if the root dataset has this info
        if len(self.xradar.ds) > 0:
            root_param_dict = self.xradar.ds.to_dict(data="array")
            instrument_parameters.update(root_param_dict["data_vars"])
            instrument_parameters.update(root_param_dict["attrs"])

        if len(instrument_parameters.keys()) > 0:
            for field in instrument_parameters.keys():
                field_dict = instrument_parameters[field]
                if isinstance(field_dict, dict):
                    if "attrs" in field_dict:
                        for param in field_dict["attrs"]:
                            field_dict[param] = field_dict["attrs"][param]
                        del field_dict["attrs"]

                    if "dims" in field_dict:
                        del field_dict["dims"]
                instrument_parameters[field] = field_dict

        return instrument_parameters

    def iter_start(self):
        """Return an iterator over the sweep start indices."""
        return (s for s in self.sweep_start_ray_index["data"])

    def iter_end(self):
        """Return an iterator over the sweep end indices."""
        return (s for s in self.sweep_end_ray_index["data"])

    def iter_start_end(self):
        """Return an iterator over the sweep start and end indices."""
        return ((s, e) for s, e in zip(self.iter_start(), self.iter_end()))

    def iter_slice(self):
        """Return an iterator which returns sweep slice objects."""
        return (slice(s, e + 1) for s, e in self.iter_start_end())

    def iter_field(self, field_name):
        """Return an iterator which returns sweep field data."""
        self.check_field_exists(field_name)
        return (self.fields[field_name]["data"][s] for s in self.iter_slice())

    def iter_azimuth(self):
        """Return an iterator which returns sweep azimuth data."""
        return (self.azimuth["data"][s] for s in self.iter_slice())

    def iter_elevation(self):
        """Return an iterator which returns sweep elevation data."""
        return (self.elevation["data"][s] for s in self.iter_slice())

    def add_field(self, field_name, dic, replace_existing=False):
        """
        Add a field to the object.

        Parameters
        ----------
        field_name : str
            Name of the field to add to the dictionary of fields.
        dic : dict
            Dictionary contain field data and metadata.
        replace_existing : bool, optional
            True to replace the existing field with key field_name if it
            exists, loosing any existing data. False will raise a ValueError
            when the field already exists.

        """
        # check that the field dictionary to add is valid
        if field_name in self.fields and replace_existing is False:
            err = f"A field with name: {field_name} already exists"
            raise ValueError(err)
        if "data" not in dic:
            raise KeyError("dic must contain a 'data' key")
        if dic["data"].shape != (self.nrays, self.ngates):
            t = (self.nrays, self.ngates)
            err = "'data' has invalid shape, should be (%i, %i)" % t
            raise ValueError(err)
        # add the field
        self.fields[field_name] = dic
        for sweep in range(self.nsweeps):
            sweep_ds = (
                self.xradar[f"sweep_{sweep}"].to_dataset().drop_duplicates("azimuth")
            )
            sweep_ds[field_name] = (
                ("azimuth", "range"),
                self.fields[field_name]["data"][self.get_slice(sweep)],
            )
            attrs = dic.copy()
            del attrs["data"]
            sweep_ds[field_name].attrs = attrs
            self.xradar[f"sweep_{sweep}"].ds = sweep_ds
        return

    def get_field(self, sweep, field_name, copy=False):
        """
        Return the field data for a given sweep.

        When used with :py:func:`get_gate_x_y_z` this method can be used to
        obtain the data needed for plotting a radar field with the correct
        spatial context.

        Parameters
        ----------
        sweep : int
            Sweep number to retrieve data for, 0 based.
        field_name : str
            Name of the field from which data should be retrieved.
        copy : bool, optional
            True to return a copy of the data. False, the default, returns
            a view of the data (when possible), changing this data will
            change the data in the underlying Radar object.

        Returns
        -------
        data : array
            Array containing data for the requested sweep and field.
        """
        self.check_field_exists(field_name)
        s = self.get_slice(sweep)
        data = self.fields[field_name]["data"][s]
        if copy:
            return data.copy()
        else:
            return data

    def check_field_exists(self, field_name):
        """
        Check that a field exists in the fields dictionary.

        If the field does not exist raise a KeyError.

        Parameters
        ----------
        field_name : str
            Name of field to check.

        """
        if field_name not in self.fields:
            raise KeyError("Field not available: " + field_name)
        return

    def get_gate_x_y_z(self, sweep, edges=False, filter_transitions=False):
        """
        Return the x, y and z gate locations in meters for a given sweep.

        With the default parameter this method returns the same data as
        contained in the gate_x, gate_y and gate_z attributes but this method
        performs the gate location calculations only for the specified sweep
        and therefore is more efficient than accessing this data through these
        attribute.

        When used with :py:func:`get_field` this method can be used to obtain
        the data needed for plotting a radar field with the correct spatial
        context.

        Parameters
        ----------
        sweep : int
            Sweep number to retrieve gate locations from, 0 based.
        edges : bool, optional
            True to return the locations of the gate edges calculated by
            interpolating between the range, azimuths and elevations.
            False (the default) will return the locations of the gate centers
            with no interpolation.
        filter_transitions : bool, optional
            True to remove rays where the antenna was in transition between
            sweeps. False will include these rays. No rays will be removed
            if the antenna_transition attribute is not available (set to None).

        Returns
        -------
        x, y, z : 2D array
            Array containing the x, y and z, distances from the radar in
            meters for the center (or edges) for all gates in the sweep.

        """
        # Check to see if the data needs to be georeferenced
        if "x" not in self.xradar[f"sweep_{sweep}"].coords:
            self.combined_sweeps = self.combined_sweeps.xradar.georeference()

        data = self.combined_sweeps.sel(sweep_number=sweep)
        return data["x"].values, data["y"].values, data["z"].values

    def init_gate_x_y_z(self):
        """Initialize or reset the gate_{x, y, z} attributes."""

        ranges = self.range["data"]
        azimuths = self.azimuth["data"]
        elevations = self.elevation["data"]
        cartesian_coords = antenna_vectors_to_cartesian(
            ranges, azimuths, elevations, edges=False
        )

        if not hasattr(self, "gate_x"):
            self.gate_x = dict()

        if not hasattr(self, "gate_y"):
            self.gate_y = dict()

        if not hasattr(self, "gate_z"):
            self.gate_z = dict()

        self.gate_x = dict(data=cartesian_coords[0])
        self.gate_y = dict(data=cartesian_coords[1])
        self.gate_z = dict(data=cartesian_coords[2])

    def init_gate_alt(self):
        if not hasattr(self, "gate_altitude"):
            self.gate_altitude = dict()

        try:
            self.gate_altitude = dict(data=self.altitude["data"] + self.gate_z["data"])
        except ValueError:
            self.gate_altitude = dict(
                data=np.mean(self.altitude["data"]) + self.gate_z["data"]
            )

    def _combine_sweeps(self):
        # Loop through and extract the different datasets
        ds_list = []
        for sweep in self.sweep_group_names:
            ds_list.append(self.xradar[sweep].ds.drop_duplicates("azimuth"))

        # Merge based on the sweep number
        merged = concat(ds_list, dim="sweep_number")

        # Stack the sweep number and azimuth together
        stacked = merged.stack(gates=["sweep_number", "azimuth"]).transpose()

        # Drop the missing gates
        cleaned = stacked.where(stacked.time == stacked.time.dropna("gates"))

        # Add in number of gates variable
        cleaned["ngates"] = ("gates", np.arange(len(cleaned.gates)))

        # Return the non-missing times, ensuring valid data is returned
        return cleaned

    def add_filter(self, gatefilter, replace_existing=False, include_fields=None):
        """
        Updates the radar object with an applied gatefilter provided
        by the user that masks values in fields within the radar object.

        Parameters
        ----------
        gatefilter : GateFilter
            GateFilter instance. This filter will exclude equal to
            the conditions provided in the gatefilter and mask values
            in fields specified or all fields if include_fields is None.
        replace_existing : bool, optional
            If True, replaces the fields in the radar object with
            copies of those fields with the applied gatefilter.
            False will return new fields with the appended 'filtered_'
            prefix.
        include_fields : list, optional
            List of fields to have filtered applied to. If none, all
            fields will have applied filter.

        """
        # If include_fields is None, sets list to all fields to include.
        if include_fields is None:
            include_fields = [*self.fields.keys()]

        try:
            # Replace current fields with masked versions with applied gatefilter.
            if replace_existing:
                for field in include_fields:
                    self.fields[field]["data"] = np.ma.masked_where(
                        gatefilter.gate_excluded, self.fields[field]["data"]
                    )
            # Add new fields with prefix 'filtered_'
            else:
                for field in include_fields:
                    field_dict = copy.deepcopy(self.fields[field])
                    field_dict["data"] = np.ma.masked_where(
                        gatefilter.gate_excluded, field_dict["data"]
                    )
                    self.add_field(
                        "filtered_" + field, field_dict, replace_existing=True
                    )

        # If fields don't match up throw an error.
        except KeyError:
            raise KeyError(
                field + " not found in the original radar object, "
                "please check that names in the include_fields list "
                "match those in the radar object."
            )
        return

    def get_nyquist_vel(self, sweep, check_uniform=True):
        """
        Return the Nyquist velocity in meters per second for a given sweep.

        Raises a LookupError if the Nyquist velocity is not available, an
        Exception is raised if the velocities are not uniform in the sweep
        unless check_uniform is set to False.

        Parameters
        ----------
        sweep : int
            Sweep number to retrieve data for, 0 based.
        check_uniform : bool
            True to check to perform a check on the Nyquist velocities that
            they are uniform in the sweep, False will skip this check and
            return the velocity of the first ray in the sweep.

        Returns
        -------
        nyquist_velocity : float
            Array containing the Nyquist velocity in m/s for a given sweep.

        """
        s = self.get_slice(sweep)
        try:
            nyq_vel = self.instrument_parameters["nyquist_velocity"]["data"][s]
        except TypeError:
            raise LookupError("Nyquist velocity unavailable")
        if check_uniform:
            if np.any(nyq_vel != nyq_vel[0]):
                raise Exception("Nyquist velocities are not uniform in sweep")
        return float(nyq_vel[0])

    def get_start(self, sweep):
        """Return the starting ray index for a given sweep."""
        return int(self.combined_sweeps.ngates.sel(sweep_number=sweep).min())

    def get_end(self, sweep):
        """Return the ending ray for a given sweep."""
        return self.sweep_end_ray_index["data"][sweep]

    def get_start_end(self, sweep):
        """Return the starting and ending ray for a given sweep."""
        return self.get_start(sweep), self.get_end(sweep)

    def get_slice(self, sweep):
        """Return a slice for selecting rays for a given sweep."""
        start, end = self.get_start_end(sweep)
        return slice(start, end + 1)

    def _find_fields(self, ds):
        fields = {}
        for field in self.combined_sweeps.variables:
            if self.combined_sweeps[field].dims == ("gates", "range"):
                fields[field] = {
                    "data": self.combined_sweeps[field].values,
                    **self.combined_sweeps[field].attrs,
                }
        return fields

    def get_azimuth(self, sweep, copy=False):
        """
        Return an array of azimuth angles for a given sweep.

        Parameters
        ----------
        sweep : int
            Sweep number to retrieve data for, 0 based.
        copy : bool, optional
            True to return a copy of the azimuths. False, the default, returns
            a view of the azimuths (when possible), changing this data will
            change the data in the underlying Radar object.

        Returns
        -------
        azimuths : array
            Array containing the azimuth angles for a given sweep.

        """
        s = self.get_slice(sweep)
        azimuths = self.azimuth["data"][s]
        if copy:
            return azimuths.copy()
        else:
            return azimuths


def _point_data_factory(grid, coordinate):
    """Return a function which returns the locations of all points."""

    def _point_data():
        """The function which returns the locations of all points."""
        reg_x = grid.x["data"]
        reg_y = grid.y["data"]
        reg_z = grid.z["data"]
        if coordinate == "x":
            return np.tile(reg_x, (len(reg_z), len(reg_y), 1)).swapaxes(2, 2)
        elif coordinate == "y":
            return np.tile(reg_y, (len(reg_z), len(reg_x), 1)).swapaxes(1, 2)
        else:
            assert coordinate == "z"
            return np.tile(reg_z, (len(reg_x), len(reg_y), 1)).swapaxes(0, 2)

    return _point_data


def _point_lon_lat_data_factory(grid, coordinate):
    """Return a function which returns the geographic locations of points."""

    def _point_lon_lat_data():
        """The function which returns the geographic point locations."""
        x = grid.point_x["data"]
        y = grid.point_y["data"]
        projparams = grid.get_projparams()
        geographic_coords = cartesian_to_geographic(x, y, projparams)
        # Set point_latitude['data'] when point_longitude['data'] is evaluated
        # and vice-versa.  This ensures that both attributes contain data from
        # the same map projection and that the map projection only needs to be
        # evaluated once.
        if coordinate == 0:
            grid.point_latitude["data"] = geographic_coords[1]
        else:
            grid.point_longitude["data"] = geographic_coords[0]
        return geographic_coords[coordinate]

    return _point_lon_lat_data


def _point_altitude_data_factory(grid):
    """Return a function which returns the point altitudes."""

    def _point_altitude_data():
        """The function which returns the point altitudes."""
        return grid.origin_altitude["data"][0] + grid.point_z["data"]

    return _point_altitude_data
