""" Unit tests for the decorators module. """

from numpy.testing import assert_raises

from pyart.pkg_util.decorators import requires
from pyart.pkg_util.exceptions import MissingOptionalDepedency


HAS_PACKAGE1 = True     # package 1 is installed
HAS_PACKAGE2 = False    # package 2 is not installed


@requires('Package 1', HAS_PACKAGE1)
def function_that_requires_package1():
    return 1


@requires('Package 2', HAS_PACKAGE2)
def function_that_requires_package2():
    return 2


class ExampleClass(object):

    value = 9

    @requires('Package 1', HAS_PACKAGE1)
    def method1(self):
        return self.value

    @requires('Package 2', HAS_PACKAGE2)
    def method2(self):
        return self.value


def test_requires_met():

    assert function_that_requires_package1() == 1

    example = ExampleClass()
    assert example.method1() == 9


def test_requires_unmet():

    # accessing the function/method with unmet requirement should not raise an
    # exception, only calling them
    assert hasattr(function_that_requires_package2, '__call__')
    example = ExampleClass()
    assert hasattr(example.method2, '__call__')

    assert_raises(MissingOptionalDepedency, function_that_requires_package2)
    assert_raises(MissingOptionalDepedency, example.method2)
