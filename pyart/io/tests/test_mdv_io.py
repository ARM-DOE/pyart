""" Unit Tests for Py-ART's io/mdv_io.py module. """

import os
import tempfile
import warnings
from datetime import datetime

import numpy as np
from numpy.testing import assert_raises, assert_warns

import pyart


class Mdv_io_Tests(object):
    """
    Class for declaring unit tests for the io/mdv_io.py module.

    self.tmpfile is a temporary file which will be removed at the end
    of the test if the file exists.
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

    def test_mdv_file_read_write_radar(self):
        with warnings.catch_warnings(record=True) as w:
            mdvfile_orig = pyart.io.mdv.MdvFile(pyart.testing.MDV_PPI_FILE)
            mdvfile_orig.read_all_fields()

            mdvfile_orig.write(self.tmpfile)
            mdvfile_orig.close()
            # check that a UserWarning was issued since zlib compression used
            assert len(w) == 1
            assert issubclass(w[-1].category, UserWarning)

            mdvfile = pyart.io.mdv.MdvFile(self.tmpfile)
            self.check_mdvfile_ppi(mdvfile)

    def test_read_write_file_objects(self):
        # read from and write two a file object
        f = open(pyart.testing.MDV_PPI_FILE)
        mdvfile = pyart.io.mdv.MdvFile(f)
        mdvfile.read_all_fields()
        self.check_mdvfile_ppi(mdvfile)

        # write out the file using an file handler
        f2 = open(self.tmpfile, 'w')
        mdvfile.write(f2)
        f.close()
        f2.close()

        # re-read and check object
        mdvfile = pyart.io.mdv.MdvFile(self.tmpfile)
        self.check_mdvfile_ppi(mdvfile)

    @staticmethod
    def check_mdvfile_ppi(mdvfile):
        # check some parameters
        assert mdvfile.master_header['time_end'] == 1305889595
        assert mdvfile.master_header['chunk_hdr_offset'] == 2464
        assert len(mdvfile.field_headers) == 1
        assert mdvfile.field_headers[0]['bad_data_value'] == 0.0
        assert len(mdvfile.vlevel_headers) == 1
        assert mdvfile.vlevel_headers[0]['record_len1'] == 1016
        assert mdvfile.chunk_headers[0]['chunk_id'] == 3
        assert mdvfile.calib_info['day'] == 15
        assert mdvfile.radar_info['scan_mode'] == 8
        az_deg, range_km, el_deg = mdvfile._calc_geometry()
        assert np.all(az_deg == np.arange(360))
        assert len(range_km) == 110
        assert round(range_km[10], 2) == 1.32
        assert len(el_deg) == 1
        assert el_deg[0] == 0.75
        assert mdvfile.times['time_begin'] == datetime(2011, 5, 20, 11, 1)
        assert len(mdvfile.elevations) == 0

    def test_mdv_file_read_write_radar_rhi(self):
        # write and read the RHI file which contains a elevation chunk
        mdvfile = pyart.io.mdv.MdvFile(pyart.testing.MDV_RHI_FILE)
        mdvfile.read_all_fields()
        mdvfile.write(self.tmpfile)
        mdvfile.close()
        mdvfile2 = pyart.io.mdv.MdvFile(self.tmpfile)

        # verify that the elevations are similar
        assert len(mdvfile2.elevations) == len(mdvfile.elevations)
        assert mdvfile2.elevations[0] == mdvfile.elevations[0]


def test_time_dic_to_unixtime():
    r = pyart.io.mdv_io.time_dict_to_unixtime(
        {'data': [20], 'units': 'seconds since 1970-01-01 00:00:00'})
    assert abs(r - 20.) <= 0.01


def test_empty_mdvfile():
    mdvfile = pyart.io.mdv.MdvFile(None)
    # XXX why are the following defined if they are never called?
    # Should they raise a NotImplemented or other Error rather than returning
    # incorrect data?
    mdvfile._get_chunk_header()
    mdvfile._get_radar_info()
    mdvfile._get_calib()
    mdvfile._get_compression_info()
    mdvfile._get_levels_info(10)


def test_mdv_methods():
    mdvfile = pyart.io.mdv.MdvFile(pyart.testing.MDV_PPI_FILE)
    # XXX why are this available
    mdvfile._time_dict_into_header()


def test_radar_exclude_fields():
    # skip fields
    radar = pyart.io.read_mdv(
        pyart.testing.MDV_PPI_FILE, exclude_fields=['reflectivity'])
    assert 'reflectivity' not in radar.fields


def test_rhi_cart():
    mdvfile = pyart.io.mdv.MdvFile(pyart.testing.MDV_RHI_FILE)
    assert mdvfile.projection == 'rhi'
    carts = mdvfile._make_carts_dict()
    assert carts['x'].shape == (1, 283, 125)
    assert carts['y'].shape == (1, 283, 125)
    assert carts['z'].shape == (1, 283, 125)
    assert round(carts['x'][0, 1, 2], 2) == -52.63
    assert round(carts['y'][0, 1, 2], 2) == -332.31
    assert round(carts['z'][0, 1, 2], 2) == 121.47
