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
    An object for holding gridded Radar data.

    Parameters
    ----------
    fields : dict
        Dictionary of field dictionaries.
    axes : dict
        Dictionary of axes dictionaries.
    metadata : dict
        Dictionary of metadata.

    Attributes
    ----------
    fields: dict
        Dictionary of field dictionaries.
    axes: dict
        Dictionary of axes dictionaries.
    metadata: dict
        Dictionary of metadata.

    """
    def __init__(self, fields, axes, metadata):
        """ Initalize object. """
        self.fields = fields
        self.metadata = metadata
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

    def add_field(self, field_name, field_dict):
        """ Add field to Grid object.

        Parameters
        ----------
        field_name : str
            Name of field to add.
        field_dict : dict
            Dictionary containing field data and metadata.
        """

        nz, ny, nx = self.fields[self.fields.keys[0]]['data'].shape

        if 'data' not in field_dict:
            raise KeyError('Field dictionary must contain a "data" key')
        if field_name in self.fields:
            raise ValueError('A field named %s already exists' % (field_name))
        if field_dict['data'].shape != (nz, ny, nx):
            raise ValueError('Field has invalid shape')

        self.fields[field_name] = field_dict
        return
