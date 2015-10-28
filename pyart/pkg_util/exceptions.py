"""
pyart.pkg_util.exceptions
=========================

Various exceptions.

.. autosummary::
    :toctree: generated/

    MissingOptionalDepedency

"""


class MissingOptionalDepedency(Exception):
    """ Exception raised when a optional depency is needed by not found. """
    pass
