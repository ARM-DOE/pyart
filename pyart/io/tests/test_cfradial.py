""" Unit Tests for Py-ART's io/cfradial.py module. """

from __future__ import print_function

import warnings

import numpy as np
from numpy.ma.core import MaskedArray
from numpy.testing import assert_array_equal, assert_almost_equal
import netCDF4
import pytest

import pyart

#################################################
# read_cfradial tests (verify radar attributes) #
#################################################

# read in the sample file and create a a Radar object
radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)


# time attribute
def test_time():
    assert 'comment' in radar.time.keys()
    assert 'long_name' in radar.time.keys()
    assert 'standard_name' in radar.time.keys()
    assert 'units' in radar.time.keys()
    assert 'calendar' in radar.time.keys()
    assert 'data' in radar.time.keys()
    assert radar.time['units'] == 'seconds since 2011-05-20T10:54:08Z'
    assert radar.time['data'].shape == (40, )
    assert_almost_equal(radar.time['data'][38], 7, 0)


# range attribute
def test_range():
    assert 'long_name' in radar.range
    assert 'standard_name' in radar.range
    assert 'meters_to_center_of_first_gate' in radar.range
    assert 'meters_between_gates' in radar.range
    assert 'units' in radar.range
    assert 'data' in radar.range
    assert 'spacing_is_constant' in radar.range
    assert radar.range['data'].shape == (42, )
    assert_almost_equal(radar.range['data'][0], 0, 0)


# fields attribute is tested later


# metadata attribute
def test_metadata():
    assert 'Conventions' in radar.metadata
    assert 'comment' in radar.metadata
    assert 'history' in radar.metadata
    assert 'institution' in radar.metadata
    assert 'instrument_name' in radar.metadata
    assert 'instrument_type' in radar.metadata
    assert 'platform_type' in radar.metadata
    assert 'primary_axis' in radar.metadata
    assert 'references' in radar.metadata
    assert 'source' in radar.metadata
    assert 'title' in radar.metadata
    assert 'version' in radar.metadata
    assert 'volume_number' in radar.metadata


# scan_type attribute
def test_scan_type():
    assert radar.scan_type == 'ppi'


# latitude attribute
def test_latitude():
    assert 'data' in radar.latitude
    assert 'standard_name' in radar.latitude
    assert 'units' in radar.latitude
    assert radar.latitude['data'].shape == (1, )
    assert_almost_equal(radar.latitude['data'], 36.0, 0)


# longitude attribute
def test_longitude():
    assert 'data' in radar.longitude
    assert 'standard_name' in radar.longitude
    assert 'units' in radar.longitude
    assert radar.longitude['data'].shape == (1, )
    assert_almost_equal(radar.longitude['data'], -98, 0)


# altitude attribute
def test_altitude():
    assert 'data' in radar.altitude
    assert 'standard_name' in radar.altitude
    assert 'units' in radar.altitude
    assert 'positive' in radar.altitude
    assert radar.altitude['data'].shape == (1, )
    assert_almost_equal(radar.altitude['data'], 214, 0)


# altitude_agl attribute
def test_altitude_agl():
    assert radar.altitude_agl is None


# sweep_number attribute
def test_sweep_number():
    assert 'standard_name' in radar.sweep_number
    assert np.all(radar.sweep_number['data'] == range(1))


# sweep_mode attribute
def test_sweep_mode():
    assert 'standard_name' in radar.sweep_mode
    assert radar.sweep_mode['data'].shape == (1, 32)
    assert radar.sweep_mode['data'].dtype.char == 'S'
    str_array = netCDF4.chartostring(radar.sweep_mode['data'])
    try:
        assert np.all(str_array == ['azimuth_surveillance'])
    except AssertionError:
        assert np.all(str_array == [b'azimuth_surveillance'])


# fixed_angle attribute
def test_fixed_angle():
    assert 'standard_name' in radar.fixed_angle
    assert 'units' in radar.fixed_angle
    assert radar.fixed_angle['data'].shape == (1, )
    assert_almost_equal(radar.fixed_angle['data'][0], 0.50, 2)


# sweep_start_ray_index attribute
def test_sweep_start_ray_index():
    assert 'long_name' in radar.sweep_start_ray_index
    assert radar.sweep_start_ray_index['data'].shape == (1, )
    assert_almost_equal(radar.sweep_start_ray_index['data'][0], 0, 0)


# sweep_end_ray_index attribute
def test_sweep_end_ray_index():
    assert 'long_name' in radar.sweep_end_ray_index
    assert radar.sweep_end_ray_index['data'].shape == (1, )
    assert_almost_equal(radar.sweep_end_ray_index['data'][0], 39, 0)


# target_scan_rate attribute
def test_target_scan_rate():
    assert radar.target_scan_rate is None


# azimuth attribute
def test_azimuth():
    assert 'standard_name' in radar.azimuth
    assert 'units' in radar.azimuth
    assert_almost_equal(radar.azimuth['data'][0], 360, 0)
    assert_almost_equal(radar.azimuth['data'][10], 90., 0)


# elevation attribute
def test_elevation():
    assert 'standard_name' in radar.elevation
    assert 'units' in radar.elevation
    assert radar.elevation['data'].shape == (40, )
    assert_almost_equal(radar.elevation['data'][0], 0.48, 2)


# scan_rate attribute
def test_scan_rate():
    assert radar.scan_rate is None


# antenna_transition attribute
def test_antenna_transition():
    assert radar.antenna_transition is None


# instrument_parameters attribute
@pytest.mark.parametrize(
    "keys", ['prt_mode', 'prt', 'nyquist_velocity', 'unambiguous_range'])
def test_instument_parameters(keys):
    # instrument_parameter sub-convention
    description = 'instrument_parameters: %s' % keys
    check_instrument_parameter.description = description
    check_instrument_parameter(keys)


def check_instrument_parameter(param):
    assert param in radar.instrument_parameters
    param_dic = radar.instrument_parameters[param]
    assert param_dic['meta_group'] == 'instrument_parameters'


# ngates attribute
def test_ngates():
    assert radar.ngates == 42


# nrays attribute
def test_nrays():
    assert radar.nrays == 40


# nsweeps attribute
def test_nsweeps():
    assert radar.nsweeps == 1


####################
# fields attribute #
####################


@pytest.mark.parametrize("field", ['reflectivity_horizontal', ])
def test_field_dics(field):
    description = "field : %s, dictionary" % field
    check_field_dic.description = description
    check_field_dic(field)


def check_field_dic(field):
    """ Check that the required keys are present in a field dictionary. """
    assert 'standard_name' in radar.fields[field]
    assert 'units' in radar.fields[field]
    assert '_FillValue' in radar.fields[field]
    assert 'coordinates' in radar.fields[field]


@pytest.mark.parametrize("field", ['reflectivity_horizontal', ])
def test_field_shapes(field):
    description = "field : %s, shape" % field
    check_field_shape.description = description
    check_field_shape(field)


def check_field_shape(field):
    assert radar.fields[field]['data'].shape == (40, 42)


fields = {'reflectivity_horizontal': MaskedArray, }
@pytest.mark.parametrize(
    "field, field_type", fields.items(), ids=list(fields.keys()))
def test_field_types(field, field_type):
    description = "field : %s, type" % field
    check_field_type.description = description
    check_field_type(field, field_type)


def check_field_type(field, field_type):
    assert type(radar.fields[field]['data']) is field_type


fields = {'reflectivity_horizontal': -6.0, }
@pytest.mark.parametrize(
    "field, field_value", fields.items(), ids=list(fields.keys()))
def test_field_first_points(field, field_value):
    # these values can be found using:
    # [round(radar.fields[f]['data'][0,0]) for f in radar.fields]
    description = "field : %s, first point" % field
    check_field_first_point.description = description
    check_field_first_point(field, field_value)


def check_field_first_point(field, value):
    assert_almost_equal(radar.fields[field]['data'][0, 0], value, 0)

########################################################################
# write_cfradial tests (verify data in written netCDF matches original #
########################################################################


def test_write_ppi():
    # CF/Radial example file -> Radar object -> netCDF file
    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_ppi.nc'
        radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
        pyart.io.write_cfradial(tmpfile, radar)
        ref = netCDF4.Dataset(pyart.testing.CFRADIAL_PPI_FILE)
        dset = netCDF4.Dataset(tmpfile)
        check_dataset_to_ref(dset, ref)
        dset.close()


def test_write_ppi_unicode():
    # GitHub issue #381, unicode dtypes cause write to fail
    # CF/Radial example file -> Radar object -> netCDF file
    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_ppi_unicode.nc'
        radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
        # force the sweep mode array to be of unicode
        radar.sweep_mode['data'] = np.array(
            ['azimuth_surveillance'], dtype='U')
        pyart.io.write_cfradial(tmpfile, radar)
        ref = netCDF4.Dataset(pyart.testing.CFRADIAL_PPI_FILE)
        dset = netCDF4.Dataset(tmpfile)
        check_dataset_to_ref(dset, ref)
        dset.close()


def test_write_ppi_U1():
    # GitHub issue #381, unicode dtypes cause write to fail
    # CF/Radial example file -> Radar object -> netCDF file
    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_ppi_unicode.nc'
        radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
        # force the sweep mode array to be a masked U1 array
        radar.sweep_mode['data'] = radar.sweep_mode['data'].astype('U')
        pyart.io.write_cfradial(tmpfile, radar)
        ref = netCDF4.Dataset(pyart.testing.CFRADIAL_PPI_FILE)
        dset = netCDF4.Dataset(tmpfile)
        check_dataset_to_ref(dset, ref)
        dset.close()


def test_write_ppi_unknown_instrument_parameter_element():
    # CF/Radial example file -> Radar object -> netCDF file
    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_ppi_unknonw_ip.nc'
        radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
        radar.instrument_parameters['foobar'] = {'data': np.zeros(40)}
        pytest.warns(UserWarning, pyart.io.write_cfradial, tmpfile, radar)
        ref = netCDF4.Dataset(pyart.testing.CFRADIAL_PPI_FILE)
        dset = netCDF4.Dataset(tmpfile)
        check_dataset_to_ref(dset, ref)
        dset.close()


def test_write_rhi():
    # CF/Radial example file -> Radar object -> netCDF file
    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_rhi.nc'
        radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_RHI_FILE)
        pyart.io.write_cfradial(tmpfile, radar)
        ref = netCDF4.Dataset(pyart.testing.CFRADIAL_RHI_FILE)
        dset = netCDF4.Dataset(tmpfile)
        check_dataset_to_ref(dset, ref)
        dset.close()


def test_write_ppi_arm_time_vars():
    # CF/Radial example file -> Radar object -> netCDF file
    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_arm_time_vars.nc'
        radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_PPI_FILE)
        pyart.io.write_cfradial(tmpfile, radar, arm_time_variables=True)
        dset = netCDF4.Dataset(tmpfile)
        assert 'base_time' in dset.variables
        assert 'time_offset' in dset.variables

        base_time = dset.variables['base_time']
        assert base_time[:] == 1305888856
        # The month is abbreviated according to the current locale, for
        # en_US the complete string is 20-May-2011,10:54:16 GMT.
        # In order to support other locales test only the non-locale dependent
        # portions of the string.
        assert base_time.string.startswith('20-')
        assert base_time.string.endswith('2011,10:54:16 GMT')

        time_offset = dset.variables['time_offset']
        assert time_offset.units == 'seconds since 2011-05-20 10:54:08'
        assert_almost_equal(time_offset[10], 11, 0)
        dset.close()


def check_dataset_to_ref(dset, ref):
    """ Check that all data in Dataset is contained in the ref Dataset. """

    # file format (these are not expected to match)
    assert dset.file_format == ref.file_format

    # global attributes
    dset_attributes = dset.ncattrs()
    ref_attributes = ref.ncattrs()

    for attr in dset_attributes:
        print("Global attribute:", attr)
        assert attr in ref_attributes
        attribute_equal(dset, ref, attr, allow_str_case_diff=True)

    # cmptypes, expected to be empty
    print("Checking cmptypes...")
    assert dset.cmptypes == ref.cmptypes

    # groups, expected to be empty
    print("Checking groups...")
    assert dset.groups == ref.groups

    # dimensions (these are not expected to match)
    print("Checking dimensions")
    for dim in dset.dimensions.keys():
        print("Dimension", dim)
        assert dim in ref.dimensions.keys()

    # variables
    print("Checking variables...")
    dset_vars = dset.variables
    ref_vars = ref.variables
    for v in dset_vars.keys():
        if v in ['platform_type', 'instrument_type', 'primary_axis']:
            continue    # default value created by pyart.
        print("Variable", v)
        check_variable_to_ref(dset_vars[v], ref_vars[v])


def check_variable_to_ref(var, ref_var):
    """ Check that the data/metadata in var matches the ref. variable. """
    # check variable attributes
    for attr in var.ncattrs():
        print("Checking attribute", attr)
        assert attr in ref_var.ncattrs()
        attribute_equal(var, ref_var, attr)

    # instance variables
    assert var.dimensions == ref_var.dimensions
    assert var.dtype == ref_var.dtype
    assert var.ndim == ref_var.ndim
    assert var.shape == ref_var.shape

    if 'least_significant_digit' in dir(var):
        assert 'least_significant_digit' in dir(ref_var)
        assert var.least_significant_digit == ref_var.least_significant_digit

    # properties
    assert var.size == ref_var.size
    # netCDF4 version < 1.1.1 use the maskandscale attribute
    if hasattr(ref_var, 'maskandscale'):
        assert var.maskandscale == ref_var.maskandscale
    # netCDF4 version >= 1.1.1 use a seperate mask and scale attributes
    if hasattr(ref_var, 'mask'):
        assert var.mask == ref_var.mask
    if hasattr(ref_var, 'scale'):
        assert var.scale == ref_var.scale

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

    if type(a1) is str and allow_str_case_diff:
        assert a1.upper() == a2.upper()
    else:
        assert a1 == a2
    return


def test_auto_history_and_conventions():
    # history and Conventions metadata should be created on write if
    # they do not exist in the original radar object
    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_auto_history.nc'
        radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_RHI_FILE)
        radar.metadata.pop('Conventions')
        radar.metadata.pop('history')
        pyart.io.write_cfradial(tmpfile, radar)
        radar2 = pyart.io.read_cfradial(tmpfile)
        assert 'Conventions' in radar2.metadata
        assert 'history' in radar2.metadata


def test_rhi_sweep_end_ray_index():
    radar = pyart.io.read_cfradial(pyart.testing.CFRADIAL_RHI_FILE)
    assert 'long_name' in radar.sweep_end_ray_index
    assert radar.sweep_end_ray_index['data'].shape == (1, )
    assert_almost_equal(radar.sweep_end_ray_index['data'][0], 39, 0)


def test_delay_field_loading():
    radar = pyart.io.read_cfradial(
        pyart.testing.CFRADIAL_PPI_FILE, delay_field_loading=True)
    assert isinstance(radar.fields['reflectivity_horizontal'],
                      pyart.lazydict.LazyLoadDict)
    data = radar.fields['reflectivity_horizontal']['data']
    assert isinstance(data, MaskedArray)
    assert data.shape == (40, 42)
    assert_almost_equal(data[0, 0], -6.0, 0)


def test_create_ncvar_different_dtype():
    # test _Write_as_dtype key handling in _create_ncvar
    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_nvcar_different_dtype.nc'
        dset = netCDF4.Dataset(tmpfile, mode='w')
        dset.createDimension('x', 256)
        dic = {
            'data': np.arange(256, dtype='float32') * 0.5 + 100.,
            '_Write_as_dtype': 'u1',
            '_FillValue': 100.,
        }
        assert 'add_offset' not in dic
        assert 'scale_factor' not in dic

        pyart.io.cfradial._create_ncvar(dic, dset, 'foo', ('x'))

        assert 'add_offset' in dic
        assert 'scale_factor' in dic
        assert '_FillValue' in dic

        foo = dset.variables['foo']
        assert foo.dtype == np.dtype('uint8')
        assert_almost_equal(foo.scale_factor, 0.5)
        assert_almost_equal(foo.add_offset, 100)
        assert foo[0] is np.ma.masked
        assert_almost_equal(foo[1], 100.5)
        assert_almost_equal(foo[-1], 227.5)

        dset.close()


def test_calculate_scale_and_offset():
    data = np.arange(128, dtype='float32')
    scale, offset, fill = pyart.io.cfradial._calculate_scale_and_offset(
        {'data': data}, np.dtype('u1'))

    assert_almost_equal(scale, 0.5)
    assert_almost_equal(offset, -0.5)
    assert_almost_equal(fill, 0)

    # check that nans and infs are ignored
    data[5] = np.nan
    data[6] = np.inf
    scale, offset, fill = pyart.io.cfradial._calculate_scale_and_offset(
        {'data': data}, np.dtype('u1'))
    assert_almost_equal(scale, 0.5)
    assert_almost_equal(offset, -0.5)
    assert_almost_equal(fill, 0)

    # maximum < minimum raises ValueError
    pytest.raises(ValueError, pyart.io.cfradial._calculate_scale_and_offset,
                  {'data': data}, np.dtype('u1'), 100., 99)

    # maximum == minimum issues warning
    with warnings.catch_warnings(record=True) as w:
        assert len(w) == 0
        pyart.io.cfradial._calculate_scale_and_offset(
            {'data': data}, np.dtype('u1'), 100, 100)
        assert len(w) == 1
