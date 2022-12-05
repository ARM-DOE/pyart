""" Unit tests for Py-ART's util/datetime_utils.py module. """

import cftime
import netCDF4

import pyart

radar = pyart.testing.make_empty_ppi_radar(10, 36, 1)


def test_datetime_from_radar():
    test_date = cftime.DatetimeGregorian(1989, 1, 1, 0, 0, 1, 0)

    date = pyart.util.datetime_from_radar(radar, epoch=False)
    assert date == test_date

    date_epoch = pyart.util.datetime_from_radar(radar, epoch=True)
    assert date_epoch == test_date


def test_datetimes_from_radar():
    test_date_first = cftime.DatetimeGregorian(1989, 1, 1, 0, 0, 1, 0)
    test_date_last = cftime.DatetimeGregorian(1989, 1, 1, 0, 0, 36, 0)

    dates = pyart.util.datetimes_from_radar(radar, epoch=False)
    assert dates[0] == test_date_first
    assert dates[-1] == test_date_last

    dates_epoch = pyart.util.datetimes_from_radar(radar, epoch=True)
    assert dates_epoch[0] == test_date_first
    assert dates_epoch[-1] == test_date_last


ds = netCDF4.Dataset(pyart.testing.CFRADIAL_PPI_FILE)


def test_datetime_from_dataset():
    test_date = cftime.DatetimeGregorian(2011, 5, 20, 10, 54, 16, 0)

    date = pyart.util.datetime_from_dataset(ds, epoch=False)
    assert date == test_date

    date_epoch = pyart.util.datetime_from_dataset(ds, epoch=True)
    assert date_epoch == test_date


def test_datetimes_from_dataset():
    test_date_first = cftime.DatetimeGregorian(2011, 5, 20, 10, 54, 16, 0)
    test_date_last = cftime.DatetimeGregorian(2011, 5, 20, 10, 54, 15, 0)

    dates = pyart.util.datetimes_from_dataset(ds, epoch=False)
    assert dates[0] == test_date_first
    assert dates[-1] == test_date_last

    dates_epoch = pyart.util.datetimes_from_dataset(ds, epoch=True)
    assert dates_epoch[0] == test_date_first
    assert dates_epoch[-1] == test_date_last


grid = pyart.testing.make_empty_grid(
    (2, 5, 5), ((0, 1000), (-2500, 2500), (-2500, 2500))
)


def test_datetime_from_grid():
    test_date = cftime.DatetimeGregorian(2000, 1, 1, 0, 0, 0, 0)

    date = pyart.util.datetime_from_grid(grid, epoch=False)
    assert date == test_date

    date_epoch = pyart.util.datetime_from_grid(grid, epoch=True)
    assert date_epoch == test_date
