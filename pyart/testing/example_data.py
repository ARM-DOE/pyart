import importlib.resources

import pooch

DATASETS = pooch.create(
    path=pooch.os_cache("pyart-datasets"),
    base_url="https://adc.arm.gov/pyart/example_data/",
    env="PYART_DATASETS_DIR",
)


with open(importlib.resources.files("pyart.testing") / "registry.txt") as registry_file:
    DATASETS.load_registry(registry_file)


def locate():
    """The absolute path to the sample data storage location on disk.
    This is where the data are saved on your computer. The location is
    dependent on the operating system. The folder locations are defined by the
    ``appdirs``  package (see the `appdirs documentation
    <https://github.com/ActiveState/appdirs>`__).
    The location can be overwritten by the ``PYTHIA_DATASETS_DIR`` environment
    variable to the desired destination.
    Returns
    -------
    path : str
        The local data storage location.
    """
    return str(DATASETS.abspath)


def get_test_data(file):
    """
    Accesses the desired test data storage
    Parameters
    ----------
    file = str
        The name of the desired file
    Returns
    -------
    path_to_file = str
        Local path to the desired file
    """
    return DATASETS.fetch(file)
