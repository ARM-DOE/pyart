""" Unit Tests for Py-ART's io/mdv_grid.py module. """

import datetime
from io import BytesIO

import numpy as np
from numpy.testing import assert_raises, assert_warns, assert_almost_equal

import pyart


class Mdv_grid_Tests(object):
    """
    Class for declaring unit tests for the io/mdv_grid.py module which
    require a temporary file which will be removed at the end of the test.
    """

    def test_write_read_target(self):
        # write and read in target grid
        original_grid = pyart.testing.make_target_grid()
        tmpfile = BytesIO()
        pyart.io.write_grid_mdv(tmpfile, original_grid)
        tmpfile.seek(0)
        grid = pyart.io.read_grid_mdv(tmpfile, file_field_names=True)

        # check the contents of the grid which was read
        self.check_target_reflectivity_field(grid)
        attrs_to_check = [
            'time', 'origin_latitude', 'origin_longitude', 'origin_altitude',
            'x', 'y', 'z']
        for attr in attrs_to_check:
            self.check_attr_dics(grid, original_grid, attr)

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
    def check_attr_dics(grid1, grid2, attr_name):
        print("Checking:", attr_name)
        dic1 = getattr(grid1, attr_name)
        dic2 = getattr(grid2, attr_name)
        # check that the data and units keys are similar
        assert np.allclose(dic1['data'], dic2['data'])
        assert dic1['data'].shape == dic2['data'].shape
        assert dic1['data'].dtype == dic2['data'].dtype
        if 'units' in dic1:
            assert dic1['units'] == dic2['units']

    def test_write_read_target_delay(self):
        # write and read in the target grid with delayed field loading
        original_grid = pyart.testing.make_target_grid()
        tmpfile = BytesIO()
        pyart.io.write_grid_mdv(tmpfile, original_grid)
        tmpfile.seek(0)
        grid = pyart.io.read_grid_mdv(tmpfile, file_field_names=True,
                                      delay_field_loading=True)

        # check the contents of the grid which was read
        self.check_target_reflectivity_field(grid)
        attrs_to_check = [
            'time', 'origin_latitude', 'origin_longitude', 'origin_altitude',
            'x', 'y', 'z']
        for attr in attrs_to_check:
            self.check_attr_dics(grid, original_grid, attr)

    def test_write_read_target_modified(self):
        # write and read in target grid with modifications
        original_grid = pyart.testing.make_target_grid()
        original_grid.metadata['radar_0_lon'] = -98.1
        original_grid.metadata['radar_0_lat'] = 36.74
        original_grid.metadata['radar_0_alt'] = 300.
        original_grid.metadata['instrument_name'] = 'testtesttest'
        original_grid.fields['reflectivity']['standard_name'] = 'fieldname'
        tmpfile = BytesIO()
        pyart.io.write_grid_mdv(tmpfile, original_grid)
        tmpfile.seek(0)
        grid = pyart.io.read_grid_mdv(tmpfile, file_field_names=True)

        # check the contents of the grid which was read
        self.check_target_reflectivity_field(grid)
        attrs_to_check = [
            'time', 'origin_latitude', 'origin_longitude', 'origin_altitude',
            'x', 'y', 'z']
        for attr in attrs_to_check:
            self.check_attr_dics(grid, original_grid, attr)
        assert grid.metadata['instrument_name'] == 'testtesttest'

    def test_read_write_mask(self):
        # Read and write a grid with field which is partially masked
        original_grid = pyart.testing.make_empty_grid(
            (20, 20, 20), ((0, 1), (0, 1), (0, 1)))
        fdata = np.ma.ones((20, 20, 20), dtype=np.float32)
        fdata[0, 0, 1] = np.ma.masked
        original_grid.fields['field_one'] = {'data': fdata}
        tmpfile = BytesIO()
        pyart.io.write_grid_mdv(tmpfile, original_grid)
        tmpfile.seek(0)

        # check the new grid object
        grid = pyart.io.read_grid_mdv(tmpfile, file_field_names=True)
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
        tmpfile = BytesIO()
        assert_warns(
            UserWarning, pyart.io.write_grid_mdv, tmpfile, grid)
        tmpfile.seek(0)
        rgrid = pyart.io.read_grid_mdv(tmpfile)
        assert len(grid.z['data']) == 123
        assert len(rgrid.z['data']) == 122

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
        tmpfile = BytesIO()
        pyart.io.write_grid_mdv(tmpfile, original_grid)
        tmpfile.seek(0)

        # verify that the data get read in correctly
        # use delay_field_loading to prevent closing of tmpfile
        grid = pyart.io.read_grid_mdv(
            tmpfile, file_field_names=True, delay_field_loading=True)
        assert np.all(grid.fields['field_one']['data'] == 105)
        assert np.all(grid.fields['field_two']['data'] == 22)
        assert np.all(grid.fields['field_three']['data'] == 100)
        tmpfile.seek(0)

        # vertify that the data is encoded to the correct type in the file
        mdvfile = pyart.io.mdv_common.MdvFile(tmpfile)
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
        tmpfile = BytesIO()
        assert_raises(TypeError, pyart.io.write_grid_mdv, tmpfile, grid)


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

    assert grid.x['units'] == 'degree_E'
    assert_almost_equal(grid.x['data'][0], -129.99, 2)
    assert grid.y['units'] == 'degree_N'
    assert_almost_equal(grid.y['data'][0], 20.01, 2)
