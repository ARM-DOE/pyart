""" Test NetCDF to NetCDF conversion through a Radar object """

import tempfile
import os

import netCDF4
from numpy.testing import assert_array_equal
import pyart


NETCDF3_FILE = 'cfrad.20080604_002217_000_SPOL_v36_SUR_netcdf3.nc'
NETCDF4_FILE = 'cfrad.20080604_002217_000_SPOL_v36_SUR_netcdf4.nc'


def test_netcdf4_to_netcdf():
    # CF/Radial example file -> Radar object -> netCDF file
    # Check that all data in resuling netCDF file matches that in original
    tmpfile = tempfile.mkstemp(suffix='.nc', dir='.')[1]
    radar = pyart.io.read_netcdf(NETCDF4_FILE)
    pyart.io.write_netcdf(tmpfile, radar)
    ref = netCDF4.Dataset(NETCDF4_FILE)
    dset = netCDF4.Dataset(tmpfile)
    check_dataset_to_ref(dset, ref)
    os.remove(tmpfile)


def test_netcdf3_to_netcdf():
    # CF/Radial example file -> Radar object -> netCDF file
    # Check that all data in resuling netCDF file matches that in original
    tmpfile = tempfile.mkstemp(suffix='.nc', dir='.')[1]
    radar = pyart.io.read_netcdf(NETCDF3_FILE)
    pyart.io.write_netcdf(tmpfile, radar)
    ref = netCDF4.Dataset(NETCDF3_FILE)
    dset = netCDF4.Dataset(tmpfile)
    check_dataset_to_ref(dset, ref)
    os.remove(tmpfile)


def check_dataset_to_ref(dset, ref):
    """ Check that all data in Dataset is contained in the ref Dataset. """

    # file format (these are not expected to match)
    #assert dset.file_format == ref.file_format

    # global attributes
    dset_attributes = dset.ncattrs()
    ref_attributes = ref.ncattrs()

    for attr in dset_attributes:
        print "Global attribute:", attr
        assert attr in ref_attributes
        attribute_equal(dset, ref, attr, allow_str_case_diff=True)

    # cmptypes, expected to be empty
    print "Checking cmptypes..."
    assert dset.cmptypes == ref.cmptypes

    # groups, expected to be empty
    print "Checking groups..."
    assert dset.groups == ref.groups

    # dimensions (these are not expected to match)
    #print "Checking dimensions"
    #for dim in dset.dimensions.keys():
    #    print "Dimension", dim
    #    assert dim in ref.dimensions.keys()

    # variables
    print "Checking variables..."
    dset_vars = dset.variables
    ref_vars = ref.variables
    for v in dset_vars.keys():
        # XXX these variables do not have the correct size/value
        if v in ['time_coverage_start', 'time_coverage_end', 'volume_number']:
            continue

        print "Variable", v
        check_variable_to_ref(dset_vars[v], ref_vars[v])


def check_variable_to_ref(var, ref_var):
    """ Check that the data/metadata in var matches the ref. variable. """
    # check variable attributes
    for attr in var.ncattrs():
        print "Checking attribute", attr
        if attr == '_FillValue':
            continue
        if attr == 'calendar':  # we add a calendar attribute to time var
            continue

        assert attr in ref_var.ncattrs()
        attribute_equal(var, ref_var, attr)

    # instance variables
    #assert var.dimensions == ref_var.dimensions  # This does not match for
                                                  # unicode/strings
    #assert var.dtype == ref_var.dtype  # These do not match, the resulting
                                        # dataset dtype match the sliced
                                        # dataset dtype
    assert var.ndim == ref_var.ndim
    assert var.shape == ref_var.shape

    if 'least_significant_digit' in dir(var):
        assert 'least_significant_digit' in dir(ref_var)
        assert var.least_significant_digit == ref_var.least_significant_digit

    # properties
    assert var.size == ref_var.size
    assert var.maskandscale == ref_var.maskandscale

    # data and the mask
    data = var[:]
    ref_data = ref_var[:]
    assert ('mask' in dir(data)) == ('mask' in dir(ref_data))
    if 'mask' in dir(data):
        assert_array_equal(data, ref_data)
        assert_array_equal(data.mask, ref_data.mask)
    else:
        assert_array_equal(data, ref_data)


def attribute_equal(class1, class2, key, allow_str_case_diff=True):
    """ Check that an attribute of two classes are equal. """
    a1 = getattr(class1, key)
    a2 = getattr(class2, key)

    assert type(a1) == type(a2)

    if type(a1) is unicode and allow_str_case_diff:
        assert a1.upper() == a2.upper()
    else:
        assert a1 == a2
    return
