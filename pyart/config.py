"""
pyart.config
============

Py-ART configuration functions.

.. autosummary::
    :toctree: generated/

    load_config
    FileMetadata

"""

import os
import imp


def load_config(filename):
    """
    Load a configuration from a Py-ART config file.

    Parameters
    ----------
    filename : str
        Filename of Py-ART configuration file.  An example file
        can be found in the pyart source directory named `default_config.py`.
        The code in this file is executed as is which may be a security
        issue, do not load un-trusted configuration files.

    """
    global DEFAULT_METADATA
    global FILE_SPECIFIC_METADATA
    global FIELD_MAPPINGS

    config = imp.load_source('metadata_config', filename)
    DEFAULT_METADATA = config.DEFAULT_METADATA
    FILE_SPECIFIC_METADATA = config.FILE_SPECIFIC_METADATA
    FIELD_MAPPINGS = config.FIELD_MAPPINGS

_config_file = os.environ.get('PYART_CONFIG')
if _config_file is None:
    dirname = os.path.dirname(__file__)
    _config_file = os.path.join(dirname, 'default_config.py')
load_config(_config_file)


class _FileMetadata():

    def __init__(self, filetype, field_names=None, additional_metadata=None,
                 file_field_names=False, exclude_fields=None):
        """
        """

        # parse filetype parameter
        if filetype in FILE_SPECIFIC_METADATA:
            self.file_specific_metadata = FILE_SPECIFIC_METADATA[filetype]
        else:
            self.file_specific_metadata = {}

        # parse additional_metadata
        if additional_metadata is None:
            self.additional_metadata = {}
        else:
            self.additional_metadata = additional_metadata

        # parse field_names and file_field_names
        if file_field_names:
            self.field_names = None
        elif field_names is None:
            self.field_names = FIELD_MAPPINGS[filetype]
        else:
            self.field_names = field_names

        # parse exclude_fields
        if exclude_fields is None:
            self.exclude_fields = []
        else:
            self.exclude_fields = exclude_fields

    def get_metadata(self, p):
        """
        """

        # additional_metadata is queued first
        if p in self.additional_metadata:
            return self.additional_metadata[p].copy()

        # then the file specific metadata
        elif p in self.file_specific_metadata:
            return self.file_specific_metadata[p].copy()

        # and finally the default metadata
        elif p in DEFAULT_METADATA:
            return DEFAULT_METADATA[p].copy()

        # return a empty dict if the parameter is in none of the above
        else:
            return {}

    def __call__(self, p):
        """
        """
        return self.get_metadata(p)

    def get_field_name(self, file_field_name):
        """
        """

        if self.field_names is None:
            field_name = file_field_names
        elif file_field_name in self.field_names:
            field_name = self.field_names[file_field_name]
        else:
            return None     # field is not mapped

        if field_name in self.exclude_fields:
            return None     # field is excluded
        else:
            return field_name
