"""
pyart.exceptions
================

Custom Py-ART exceptions.

.. autosummary::
    :toctree: generated/

    MissingOptionalDependency

"""


class MissingOptionalDependency(Exception):
    """ Exception raised when a optional depency is needed by not found. """
    pass
