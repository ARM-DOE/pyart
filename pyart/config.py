"""
pyart.config
============

Py-ART configuration.

.. autosummary::
    :toctree: generated/

    load_config
    get_metadata
    get_fillvalue
    get_field_name
    get_field_colormap
    get_field_limits
    get_field_mapping
    FileMetadata

"""

import os
import sys
import traceback
import warnings


# the path to the default configuration file
_dirname = os.path.dirname(__file__)
_DEFAULT_CONFIG_FILE = os.path.join(_dirname, 'default_config.py')


def load_config(filename=None):
    """
    Load a Py-ART configuration from a config file.

    The default values for a number of Py-ART parameters and metadata is
    controlled by a single Python configuration file.  An self-descriping
    example of this file can be found in the Py-ART source directory named
    **default_config.py**.  These defaults can modified by setting the
    environmental variable `PYART_CONFIG` to point to a new configuration
    file. If this variable is not set then the settings contained in
    the **default_config.py** file are used.

    The code the configuration file is executed as-is with full permission,
    this may present a security issue, do not load un-trusted configuration
    files.

    The recommended method for changing these defaults is for users to
    copy this file into their home directory, rename it to .pyart_config.py,
    make any changes, and adjust their login scripts to set the PYART_CONFIG
    environmental variable to point to .pyart_config.py in their home
    directory.

    Py-ART's configuration can also be modified within a script or shell
    session using this function, the modification will last until a the end
    of the script/session or until a new configuration is loaded.

    Parameters
    ----------
    filename : str
        Filename of configuration file.  If None the default configuration
        file is loaded from the Py-ART source code directory.

    """
    if filename is None:
        filename = _DEFAULT_CONFIG_FILE

    # these are private since they should not be accessed by users or other
    # modules, use the get_ functions.
    global _DEFAULT_METADATA
    global _FILE_SPECIFIC_METADATA
    global _FIELD_MAPPINGS
    global _FILL_VALUE
    global _DEFAULT_FIELD_NAMES
    global _DEFAULT_FIELD_COLORMAP
    global _DEFAULT_FIELD_LIMITS

    if sys.version_info[:2] >= (3, 4):
        from importlib.machinery import SourceFileLoader
        cfile = SourceFileLoader('metadata_config', filename).load_module()
    else:
        import imp
        cfile = imp.load_source('metadata_config', filename)
    _DEFAULT_METADATA = cfile.DEFAULT_METADATA
    _FILE_SPECIFIC_METADATA = cfile.FILE_SPECIFIC_METADATA
    _FIELD_MAPPINGS = cfile.FIELD_MAPPINGS
    _FILL_VALUE = cfile.FILL_VALUE
    _DEFAULT_FIELD_NAMES = cfile.DEFAULT_FIELD_NAMES
    _DEFAULT_FIELD_COLORMAP = cfile.DEFAULT_FIELD_COLORMAP
    _DEFAULT_FIELD_LIMITS = cfile.DEFAULT_FIELD_LIMITS
    return

# load the configuration from the enviromental parameter if it is set
# if the load fails issue a warning and load the default config.
_config_file = os.environ.get('PYART_CONFIG')
if _config_file is None:
    load_config(_DEFAULT_CONFIG_FILE)
else:
    try:
        load_config(_config_file)
    except:
        msg = ("\nLoading configuration from PYART_CONFIG enviromental "
               "variable failed:"
               "\n--- START IGNORED TRACEBACK --- \n" +
               traceback.format_exc() +
               "\n --- END IGNORED TRACEBACK ---"
               "\nLoading default configuration")
        warnings.warn(msg)
        load_config(_DEFAULT_CONFIG_FILE)


def get_metadata(p):
    """
    Return a dictionary of metadata for a given parameter, p.

    An empty dictionary will be returned in no metadata dictionary exists for
    parameter p.
    """
    if p in _DEFAULT_METADATA:
        return _DEFAULT_METADATA[p].copy()
    else:
        return {}


def get_fillvalue():
    """
    Return the current fill value.
    """
    return _FILL_VALUE


def get_field_name(field):
    """
    Return the field name from the configuration file for a given field.
    """
    return str(_DEFAULT_FIELD_NAMES[field])


def get_field_colormap(field):
    """
    Return the colormap name from the configuration file for a field name.
    """
    if field in _DEFAULT_FIELD_COLORMAP:
        return _DEFAULT_FIELD_COLORMAP[field]
    else:
        import matplotlib.cm
        return matplotlib.cm.get_cmap().name


def get_field_limits(field, container=None, selection=0):
    """
    Return the data limits from the configuration file for a given field,
    radar and sweep

    Parameters
    ----------
    field: str
        Field name.
    container: Radar, Grid or None
        This is an optional parameter that will be use to get informations
        related to the field, like for instace nyquist velocity.
    selection: int
        selection of the data in the container, case container is a Radar this
        is the sweep to be considered

    Returns
    -------
    vmin, vmax: 2-tuplet of float
        Minimun and Maximun teorical value for field, if field is not
        in the configuration file returns (None, None).
    """
    if field in _DEFAULT_FIELD_LIMITS:
        limits = _DEFAULT_FIELD_LIMITS[field]
        if callable(limits):
            limits = limits(container, selection)
        return limits
    else:
        return None, None


def get_field_mapping(filetype):
    """
    Return a copy of the default field mapping for a given file type.

    Parameters
    ----------
    filetype : str
        Filetype to return field mappings for.

    Returns
    -------
    field_mappings : dict
        Dictionary mapping field names from one type to another.
    """
    return _FIELD_MAPPINGS[filetype].copy()


class FileMetadata():
    """
    A class for accessing metadata needed when reading files.

    Parameters
    ----------
    filetype: str
        Type of file being read.
    field_names : dict
        Dictionary mapping file field names to radar fields names.
    additional_metadata : dict of dicts
        Additional metadata to use during read.
    file_field_names : bool
        True to keep the field names in the file.
    exclude_fields : list of strings
        Fields to exclude during readings.
    include_fields : list of strings
        Fields to include during readings.
    """

    def __init__(self, filetype, field_names=None, additional_metadata=None,
                 file_field_names=False, exclude_fields=None,
                 include_fields=None):
        """
        Initialize.
        """

        # parse filetype parameter
        if filetype in _FILE_SPECIFIC_METADATA:
            self._file_specific_metadata = _FILE_SPECIFIC_METADATA[filetype]
        else:
            self._file_specific_metadata = {}

        # parse additional_metadata
        if additional_metadata is None:
            self._additional_metadata = {}
        else:
            self._additional_metadata = additional_metadata

        # parse field_names and file_field_names
        if file_field_names:
            self._field_names = None
        elif field_names is None and filetype in _FIELD_MAPPINGS:
            # if filetype is missing there is no mapping
            self._field_names = _FIELD_MAPPINGS[filetype]
        else:
            self._field_names = field_names

        # parse exclude_fields
        if exclude_fields is None:
            self._exclude_fields = []
        else:
            self._exclude_fields = exclude_fields
 
        if include_fields is None:
            self._include_fields = None
        else:
            self._include_fields = include_fields 

    def get_metadata(self, p):
        """
        Retrieve metadata for a parameter `p`.

        Parameters
        ----------
        p : str
            Parameter to retrieve metadata for.

        Returns
        -------
        dic : dict
            Dictionary of metadata for the parameter.
        """

        # additional_metadata is queued first
        if p in self._additional_metadata:
            return self._additional_metadata[p].copy()

        # then the file specific metadata
        elif p in self._file_specific_metadata:
            return self._file_specific_metadata[p].copy()

        # and finally the default metadata
        elif p in _DEFAULT_METADATA:
            return _DEFAULT_METADATA[p].copy()

        # return a empty dict if the parameter is in none of the above
        else:
            return {}

    def __call__(self, p):
        """
        Retrieve metadata for parameter `p`.
        """
        return self.get_metadata(p)

    def get_field_name(self, file_field_name):
        """
        Return the name radar field for a given file field name

        Parameters
        ----------
        file_field_name : str
            Field name in file being read.

        Returns
        -------
        field_name : str or None
            Field name in radar object fields dictionary, None indicated
            that the field should not be included.

        """
        if self._field_names is None:
            field_name = file_field_name
        elif file_field_name in self._field_names:
            field_name = self._field_names[file_field_name]
        else:
            return None     # field is not mapped

        if field_name in self._exclude_fields:
            return None     # field is excluded
        elif self._include_fields is not None:
            if not field_name in self._include_fields:
                return None
            else:
                return field_name
        else:
            return field_name
