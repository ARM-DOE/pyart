"""
pyart.pkg_util.decorators
=========================

Various decorators.

.. autosummary::
    :toctree: generated/

    requires

"""

from ..exceptions import MissingOptionalDependency


_REQUIRES_MESSAGE = (
    "The use of '%s' requires that the package %s be installed.")


def requires(requirement_name, requirement_check):
    """
    Decorate a function which requires an optional dependency.

    Use this function as a decorator to a function or method which requires
    that a optional depency be installed.  The function can then be imported
    into the package namespace.  When called the decorated function or method
    will raise a MissingOptionalDependency error with a helpful message to
    providing information about what package must be installed to use the
    particular function or method in question

    Parameters
    ----------
    requirement_name : str
        Name of package which is required for use of the decorated function or
        method.
    requirement_check : bool
        Boolean indicating if the requirement has been met.  Typically this is
        a variable which determines if required package is installed at run
        time.

    """
    def wrap_with_requirement(func):
        """ Wrap a function/method with a check for a requirement. """
        def wrapper(*args, **kwargs):
            """ Actual wrapper around the function/method. """
            if not requirement_check:
                raise MissingOptionalDependency(
                    _REQUIRES_MESSAGE % (func.__name__, requirement_name))
            return func(*args, **kwargs)
        return wrapper
    return wrap_with_requirement
