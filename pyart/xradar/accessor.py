"""
Utilities for interfacing between xradar and Py-ART

"""


class Xradar:
    def __init__(self, xradar):
        self.xradar = xradar

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
        data = self.xradar[f"sweep_{sweep}"][field_name].values

        if copy:
            return data.copy()
        else:
            return data

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
        if "x" not in self.xradar[f"sweep_{0}"].coords:
            self.xradar = self.xradar.xradar.georeference()

        data = self.xradar[f"sweep_{sweep}"].xradar.georeference()
        return data["x"].values, data["y"].values, data["x"].values
