"""
pyart.core.grid
===============

An class for holding gridded Radar data.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    Grid

.. autosummary::
    :toctree: generated/

    _point_data_factory
    _point_lon_lat_data_factory
    _point_altitude_data_factory

"""

import numpy as np

from ..config import get_metadata
from ..lazydict import LazyLoadDict
from .transforms import cartesian_to_geographic


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
    point_x, point_y, point_z : LazyLoadDict
        The Cartesian locations of all grid points from the origin in the
        three Cartesian coordinates.  The three dimensional data arrays
        contained these attributes are calculated from the regular_x,
        regular_y, and regular_z attributes.  If these attributes are changed
        use :py:func:`init_point_x_y_z` to reset the attributes.
    point_longitude, point_latitude : LazyLoadDict
        Geographic location of each grid point. The projection parameter(s)
        defined in the `projection` attribute are used to perform an inverse
        map projection from the Cartesian grid point locations relative to
        the grid origin. If these attributes are changed use
        :py:func:`init_point_longitude_latitude` to reset the attributes.
    point_altitude : LazyLoadDict
        The altitude of each grid point as calculated from the altitude of the
        grid origin and the Cartesian z location of each grid point.  If this
        attribute is changed use :py:func:`init_gate_altitude` to reset the
        attribute.
    projection : dic or str
        Projection parameters defining the map projection used to transform
        from Cartesian to geographic coordinates.  The default dictionary sets
        the 'proj' key to 'pyart_aeqd' indicating that the native Py-ART
        azimuthal equidistant projection is used. This can be modified to
        specify a valid pyproj.Proj projparams dictionary or string.
        The special key '_include_lon_0_lat_0' is removed when interpreting
        this dictionary. If this key is present and set to True, which is
        required when proj='pyart_aeqd', then the radar longitude and
        latitude will be added to the dictionary as 'lon_0' and 'lat_0'.
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
        self.projection = {'proj': 'pyart_aeqd', '_include_lon_0_lat_0': True}

        # initialize attributes with Lazy load dictionaries
        self.init_point_x_y_z()
        self.init_point_longitude_latitude()
        self.init_point_altitude()

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

    # Attribute init/reset methods
    def init_point_x_y_z(self):
        """ Initialize or reset the point_{x, y, z} attributes. """
        self.point_x = LazyLoadDict(get_metadata('point_x'))
        self.point_x.set_lazy('data', _point_data_factory(self, 'x'))

        self.point_y = LazyLoadDict(get_metadata('point_y'))
        self.point_y.set_lazy('data', _point_data_factory(self, 'y'))

        self.point_z = LazyLoadDict(get_metadata('point_z'))
        self.point_z.set_lazy('data', _point_data_factory(self, 'z'))

    def init_point_longitude_latitude(self):
        """
        Initialize or reset the point_{longitude, latitudes} attributes.
        """
        point_longitude = LazyLoadDict(get_metadata('point_longitude'))
        point_longitude.set_lazy('data', _point_lon_lat_data_factory(self, 0))
        self.point_longitude = point_longitude

        point_latitude = LazyLoadDict(get_metadata('point_latitude'))
        point_latitude.set_lazy('data', _point_lon_lat_data_factory(self, 1))
        self.point_latitude = point_latitude

    def init_point_altitude(self):
        """ Initialize the point_altitude attribute. """
        point_altitude = LazyLoadDict(get_metadata('point_altitude'))
        point_altitude.set_lazy('data', _point_altitude_data_factory(self))
        self.point_altitude = point_altitude

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


def _point_data_factory(grid, coordinate):
    """ Return a function which returns the locations of all points.  """
    def _point_data():
        """ The function which returns the locations of all points. """
        reg_x = grid.regular_x['data']
        reg_y = grid.regular_y['data']
        reg_z = grid.regular_z['data']
        if coordinate == 'x':
            return np.tile(reg_x, (len(reg_z), len(reg_y), 1)).swapaxes(2, 2)
        elif coordinate == 'y':
            return np.tile(reg_y, (len(reg_z), len(reg_x), 1)).swapaxes(1, 2)
        else:
            assert coordinate == 'z'
            return np.tile(reg_z, (len(reg_x), len(reg_y), 1)).swapaxes(0, 2)
    return _point_data


def _point_lon_lat_data_factory(grid, coordinate):
    """ Return a function which returns the geographic locations of points. """
    def _point_lon_lat_data():
        """ The function which returns the geographic point locations. """
        x = grid.point_x['data']
        y = grid.point_y['data']
        projparams = grid.projection.copy()
        if projparams.pop('_include_lon_0_lat_0', False):
            projparams['lon_0'] = grid.origin_longitude['data'][0]
            projparams['lat_0'] = grid.origin_latitude['data'][0]
        geographic_coords = cartesian_to_geographic(x, y, projparams)
        # set the other geographic coordinate
        if coordinate == 0:
            grid.point_latitude['data'] = geographic_coords[1]
        else:
            grid.point_longitude['data'] = geographic_coords[0]
        return geographic_coords[coordinate]
    return _point_lon_lat_data


def _point_altitude_data_factory(grid):
    """ Return a function which returns the point altitudes. """
    def _point_altitude_data():
        """ The function which returns the point altitudes. """
        return grid.origin_altitude['data'][0] + grid.point_z['data']
    return _point_altitude_data
