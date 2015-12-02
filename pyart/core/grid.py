"""
pyart.core.grid
===============

An class for holding gridded Radar data.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    Grid

"""


class Grid(object):
    """
    A class for storing rectilinear gridded radar data in Cartesian coordinate.

    Parameters
    ----------
    fields : dict
        Dictionary of field dictionaries.
    metadata : dict
        Dictionary of metadata.
    axes : dict
        Dictionary of axes dictionaries.

    Attributes
    ----------
    fields: dict of dicts
        Moments from radars or other variables.
    metadata: dict
        Metadata describing the grid.
    time : dict
        Time of the grid.
    origin_longitude, origin_latitude, origin_altitude : dict
        Geographic coordinate of the origin of the grid.
    regular_x, regular_y, regular_z : dict
        Regular locations of grid points from the origin in the three
        Cartesian coordinates.
    axes : dict
        Dictionary of axes dictionaries.
        This attribute is Depreciated, it will be removed in the next Py-ART
        release.

    """
    def __init__(self, time, fields, metadata,
                 origin_latitude, origin_longitude, origin_altitude,
                 regular_x, regular_y, regular_z):
        """ Initalize object. """

        self.time = time
        self.fields = fields
        self.metadata = metadata
        self.origin_latitude = origin_latitude
        self.origin_longitude = origin_longitude
        self.origin_altitude = origin_altitude
        self.regular_x = regular_x
        self.regular_y = regular_y
        self.regular_z = regular_z

        # Depreciated axes attribute
        axes = {'time': time,
                'time_start': time,  # incorrect metadata
                'time_end': time,    # incorrect metadata
                'z_disp': regular_z,
                'y_disp': regular_y,
                'x_disp': regular_x,
                'alt': origin_altitude,
                'lat': origin_latitude,
                'lon': origin_longitude}
        self.axes = axes

        return

    def write(self, filename, format='NETCDF4', arm_time_variables=False):
        """
        Write the the Grid object to a NetCDF file.

        Parameters
        ----------
        filename : str
            Filename to save to.
        format : str, optional
            NetCDF format, one of 'NETCDF4', 'NETCDF4_CLASSIC',
            'NETCDF3_CLASSIC' or 'NETCDF3_64BIT'.
        arm_time_variables : bool
            True to write the ARM standard time variables base_time and
            time_offset. False will not write these variables.

        """
        # delayed import to avoid circular import
        from ..io.grid_io import write_grid

        write_grid(filename, self, format=format,
                   arm_time_variables=arm_time_variables)

    def add_field(self, field_name, field_dict, replace_existing=False):
        """
        Add a field to the object.

        Parameters
        ----------
        field_name : str
            Name of the field to the fields dictionary.
        field_dict : dict
            Dictionary containing field data and metadata.
        replace_existing : bool, optional
            True to replace the existing field with key field_name if it
            exists, overwriting the existing data. If False, a ValueError is
            raised if field_name already exists.

        """

        # grid fields should have dimensions of (nz, ny, nx)
        nz, ny, nx = self.fields[self.fields.keys()[0]]['data'].shape

        # checks to make sure input field dictionary is valid
        if 'data' not in field_dict:
            raise KeyError('Field dictionary must contain a "data" key')
        if field_name in self.fields and replace_existing is False:
            raise ValueError('A field named %s already exists' % (field_name))
        if field_dict['data'].shape != (nz, ny, nx):
            raise ValueError('Field has invalid shape')

        self.fields[field_name] = field_dict

        return
