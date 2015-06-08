""" Unit Tests for Py-ART's io/mdv_grid.py module. """

from __future__ import print_function

import os
import tempfile
import datetime

import numpy as np
from numpy.testing import assert_raises, assert_warns, assert_almost_equal

import pyart


class Mdv_grid_Tests(object):
    """
    Class for declaring unit tests for the io/mdv_grid.py module which
    require a temporary file which will be removed at the end of the test.
    """

    def setUp(self):
        self.tmpfile = tempfile.mkstemp(suffix='.mdv', dir='.')[1]

    def tearDown(self):
        if os.path.isfile(self.tmpfile):
            os.remove(self.tmpfile)

    def test_write_read_target(self):
        # write and read in target grid
        original_grid = pyart.testing.make_target_grid()
        pyart.io.write_grid_mdv(self.tmpfile, original_grid)
        grid = pyart.io.read_grid_mdv(self.tmpfile, file_field_names=True)

        # check the contents of the grid which was read
        self.check_target_reflectivity_field(grid)
        keys_to_check = grid.axes.keys()
        for key in keys_to_check:
            self.check_axes_dic(grid.axes, original_grid.axes, key)

    @staticmethod
    def check_target_reflectivity_field(grid):
        # check the reflectivity field data
        assert 'reflectivity' in grid.fields
        fdata = grid.fields['reflectivity']['data']

        assert fdata.shape == (2, 400, 320)
        assert fdata.dtype == np.float32

        assert np.abs(fdata[0, 0, 0] - 0.) <= 0.01
        assert np.abs(fdata[0, 60, 50] - 10.) <= 0.01
        assert np.abs(fdata[0, 110, 90] - 20.) <= 0.01
        assert np.abs(fdata[0, 160, 130] - 30.) <= 0.01

        assert np.abs(fdata[1, 0, 0] - 5.) <= 0.01
        assert np.abs(fdata[1, 60, 50] - 15.) <= 0.01
        assert np.abs(fdata[1, 110, 90] - 25.) <= 0.01
        assert np.abs(fdata[1, 160, 130] - 35.) <= 0.01

    @staticmethod
    def check_axes_dic(axes, axes2, key):
        print("Checking:", key)
        # check that the data and units keys are similar
        assert np.allclose(axes[key]['data'], axes2[key]['data'])
        assert axes[key]['data'].shape == axes2[key]['data'].shape
        assert axes[key]['data'].dtype == axes2[key]['data'].dtype
        if 'units' in axes[key]:
            assert axes[key]['units'] == axes2[key]['units']

    def test_write_read_target_delay(self):
        # write and read in the target grid with delayed field loading
        original_grid = pyart.testing.make_target_grid()
        pyart.io.write_grid_mdv(self.tmpfile, original_grid)
        grid = pyart.io.read_grid_mdv(self.tmpfile, file_field_names=True,
                                      delay_field_loading=True)

        # check the contents of the grid which was read
        self.check_target_reflectivity_field(grid)
        keys_to_check = grid.axes.keys()
        for key in keys_to_check:
            self.check_axes_dic(grid.axes, original_grid.axes, key)

    def test_write_read_target_modified(self):
        # write and read in target grid with modifications
        original_grid = pyart.testing.make_target_grid()
        del original_grid.axes['time_start']
        del original_grid.axes['time_end']
        original_grid.metadata['radar_0_lon'] = -98.1
        original_grid.metadata['radar_0_lat'] = 36.74
        original_grid.metadata['radar_0_alt'] = 300.
        original_grid.metadata['instrument_name'] = 'testtesttest'
        original_grid.fields['reflectivity']['standard_name'] = 'fieldname'
        pyart.io.write_grid_mdv(self.tmpfile, original_grid)
        grid = pyart.io.read_grid_mdv(self.tmpfile, file_field_names=True)

        # check the contents of the grid which was read
        self.check_target_reflectivity_field(grid)
        keys_to_check = ['time', 'x_disp', 'y_disp', 'z_disp', 'lat', 'lon',
                         'alt']
        for key in keys_to_check:
            self.check_axes_dic(grid.axes, original_grid.axes, key)
        assert grid.metadata['instrument_name'] == 'testtesttest'

    def test_read_write_mask(self):
        # Read and write a grid with field which is partially masked
        original_grid = pyart.testing.make_empty_grid(
            (20, 20, 20), ((0, 1), (0, 1), (0, 1)))
        fdata = np.ma.ones((20, 20, 20), dtype=np.float32)
        fdata[0, 0, 1] = np.ma.masked
        original_grid.fields['field_one'] = {'data': fdata}
        pyart.io.write_grid_mdv(self.tmpfile, original_grid)

        # check the new grid object
        grid = pyart.io.read_grid_mdv(self.tmpfile, file_field_names=True)
        assert 'field_one' in grid.fields
        fdata = grid.fields['field_one']['data']
        assert fdata.dtype == np.float32
        assert fdata.shape == (20, 20, 20)
        assert np.ma.isMA(fdata)
        assert np.ma.is_masked(fdata[0, 0, 1])
        assert not np.ma.is_masked(fdata[0, 0, 0])

    def test_raise_warning_to_many_z(self):
        # Grids greater than 122 level should raise a warning and
        # write out the first 122 levels.
        grid = pyart.testing.make_empty_grid(
            (123, 2, 2), ((0, 1), (0, 1), (0, 1)))
        grid.fields['reflectivity'] = {'data': np.zeros((123, 2, 2))}
        assert_warns(
            UserWarning, pyart.io.write_grid_mdv, self.tmpfile, grid)
        rgrid = pyart.io.read_grid_mdv(self.tmpfile)
        assert len(grid.axes['z_disp']['data']) == 123
        assert len(rgrid.axes['z_disp']['data']) == 122

    def test_write_types(self):
        # Write various data types
        original_grid = pyart.testing.make_empty_grid(
            (2, 2, 2), ((0, 1), (0, 1), (0, 1)))

        uint8_data = np.ones((2, 2, 2), dtype=np.float32) * 100. + 5
        original_grid.fields['field_one'] = {
            'data': uint8_data,
            'scale_factor': 10.,
            'add_offset': 5.,
            '_Write_as_dtype': 'uint8',
            '_FillValue': 128, }

        uint16_data = np.ones((2, 2, 2), dtype=np.float32) * 22.
        original_grid.fields['field_two'] = {
            'data': uint16_data,
            '_Write_as_dtype': 'uint16',
            'missing_value': 256, }

        float32_data = np.ones((2, 2, 2), dtype=np.float32) * 98. + 2.
        original_grid.fields['field_three'] = {
            'data': float32_data,
            'scale_factor': 10.,
            'add_offset': 2.,
            }
        pyart.io.write_grid_mdv(self.tmpfile, original_grid)

        # verify that the data get read in correctly
        grid = pyart.io.read_grid_mdv(self.tmpfile, file_field_names=True)
        assert np.all(grid.fields['field_one']['data'] == 105)
        assert np.all(grid.fields['field_two']['data'] == 22)
        assert np.all(grid.fields['field_three']['data'] == 100)

        # vertify that the data is encoded to the correct type in the file
        mdvfile = pyart.io.mdv_common.MdvFile(self.tmpfile)
        field_name = [h['field_name'] for h in mdvfile.field_headers]
        field_one = field_name.index('field_one')
        field_two = field_name.index('field_two')
        field_three = field_name.index('field_three')
        encoding = [h['encoding_type'] for h in mdvfile.field_headers]
        assert encoding[field_one] == pyart.io.mdv_common.ENCODING_INT8
        assert encoding[field_two] == pyart.io.mdv_common.ENCODING_INT16
        assert encoding[field_three] == pyart.io.mdv_common.ENCODING_FLOAT32
        mdvfile.close()

    def test_raise_error_on_invalid_dtype(self):
        #  Unsupported dtyped field should raise a TypeError
        grid = pyart.testing.make_empty_grid(
            (2, 2, 2), ((0, 1), (0, 1), (0, 1)))
        grid.fields['reflectivity'] = {
            'data': np.zeros((2, 2, 2), dtype=np.float32),
            '_Write_as_dtype': 'float16'}
        assert_raises(TypeError, pyart.io.write_grid_mdv, self.tmpfile, grid)


def test_time_dic_to_datetime():
    dt = pyart.io.mdv_grid._time_dic_to_datetime(
        {'data': [20], 'units': 'seconds since 1970-01-01 00:00:00'})
    assert dt == datetime.datetime(1970, 1, 1, 0, 0, 20)

    dt = pyart.io.mdv_grid._time_dic_to_datetime(
        {'data': [20], 'units': 'seconds since 1970-01-01 00:00:00',
         'calendar': 'standard'})
    assert dt == datetime.datetime(1970, 1, 1, 0, 0, 20)


def test_mdv_degree_grid():
    grid = pyart.io.read_grid_mdv(
        pyart.testing.MDV_GRID_FILE, file_field_names=True)

    assert 'refl' in grid.fields.keys()
    fdata = grid.fields['refl']['data']
    assert fdata.shape == (1, 1837, 3661)
    assert np.ma.is_masked(fdata[0, 0, 0])
    assert_almost_equal(fdata[0, 130, 2536], 20.0, 1)

    assert grid.axes['x_disp']['units'] == 'degree_E'
    assert_almost_equal(grid.axes['x_disp']['data'][0], -129.99, 2)
    assert grid.axes['y_disp']['units'] == 'degree_N'
    assert_almost_equal(grid.axes['y_disp']['data'][0], 20.01, 2)
