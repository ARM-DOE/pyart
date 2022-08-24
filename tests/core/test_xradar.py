from datatree import DataTree
import pyart
import xarray as xr

# Load the sample radar data
test_radar = pyart.io.read(pyart.testing.get_test_data('110635.mdv'))

def test_create_dataset_from_sweep(test_radar=test_radar):
    xradar_single_sweep = pyart.core.xradar.create_dataset_from_sweep(test_radar)
    assert isinstance(xradar_single_sweep, xr.Dataset)

def test_convert_to_xradar(test_radar=test_radar):
    xradar_dataset = pyart.core.xradar.convert_to_xradar(test_radar)
    assert isinstance(xradar_dataset, DataTree)
