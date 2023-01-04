import pathlib

from pyart.testing.example_data import DATASETS, get_test_data, locate


def test_registry():
    files = DATASETS.registry_files
    assert len(files) > 0


def test_locate():
    p = locate()
    print(p)
    assert "datasets" in p
    assert pathlib.Path(p).exists


def test_get_test_data():
    test_file = "034142.mdv"
    file_path = get_test_data(test_file)
    assert pathlib.Path(file_path).exists
