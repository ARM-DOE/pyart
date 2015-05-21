""" Unit Tests for Py-ART's io/mdv_grid.py module. """

import os
import tempfile

import numpy as np
from numpy.testing import assert_raises, assert_warns

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
        print "Checking:", key
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
        original_grid.fields['field_one'] = {
            'data': np.zeros((2, 2, 2), dtype=np.uint8),
            'scale_factor': 10.,
            'add_offset': 5.,
            '_FillValue': 128, }
        original_grid.fields['field_two'] = {
            'data': np.zeros((2, 2, 2), dtype=np.uint16),
            'missing_value': 256, }
        original_grid.fields['field_three'] = {
            'data': np.zeros((2, 2, 2), dtype=np.float32)}
        pyart.io.write_grid_mdv(self.tmpfile, original_grid)
        grid = pyart.io.read_grid_mdv(self.tmpfile, file_field_names=True)
        # XXX how should other dtypes and scale/offset be treated, where
        # should the conversion be done?
        assert np.all(grid.fields['field_one']['data'] == 5)
        assert np.all(grid.fields['field_two']['data'] == 0)
        assert np.all(grid.fields['field_three']['data'] == 0)

    def test_write_types2(self):
        # Write various data types and verify
        original_grid = pyart.testing.make_empty_grid(
            (20, 20, 20), ((0, 1), (0, 1), (0, 1)))
        original_grid.fields['field_one'] = {
            'data': np.ones((20, 20, 20), dtype=np.uint8)*1,
            'scale_factor': 10.,
            'add_offset': 5.,
            '_FillValue': 128, }
        pyart.io.write_grid_mdv(self.tmpfile, original_grid)
        grid = pyart.io.read_grid_mdv(self.tmpfile, file_field_names=True)
        # XXX how should other dtypes and scale/offset be treated, where
        # should the conversion be done?
        assert np.all(grid.fields['field_one']['data'] == 5)

    def test_raise_errors(self):
        # XXX how should other dtypes and scale/offset be treated, where
        # should the conversion be done?
        # Odd dtyped field data should raise a TypeError
        grid = pyart.testing.make_empty_grid(
            (2, 2, 2), ((0, 1), (0, 1), (0, 1)))
        grid.fields['reflectivity'] = {
            'data': np.zeros((2, 2, 2), dtype=np.float16)}
        assert_raises(
            TypeError, pyart.io.write_grid_mdv, self.tmpfile, grid)


def test_time_dic_to_unixtime():
    r = pyart.io.mdv_grid.time_dict_to_unixtime(
        {'data': [20], 'units': 'seconds since 1970-01-01 00:00:00'})
    assert abs(r - 20.) <= 0.01
