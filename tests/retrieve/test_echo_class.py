""" Unit Tests for Py-ART's retrieve/echo_class.py module. """

import numpy as np
import pytest
from numpy.testing import assert_allclose

import pyart


def test_steiner_conv_strat_default():
    grid = pyart.testing.make_storm_grid()
    eclass_default = pyart.retrieve.steiner_conv_strat(grid)
    assert np.all(
        eclass_default["data"][25]
        == np.array(
            [
                0,
                0,
                0,
                0,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                0,
                0,
                0,
                0,
            ]
        )
    )


@pytest.mark.parametrize("area_relation", ["small", "medium", "large", "sgp"])
def test_steiner_conv_strat_modify_area(area_relation):
    grid = pyart.testing.make_storm_grid()
    eclass = pyart.retrieve.steiner_conv_strat(grid, area_relation=area_relation)
    assert eclass["data"].min() == 0
    assert eclass["data"].max() == 2


def test_conv_strat__yuter_default():
    grid = pyart.testing.make_storm_grid()
    dict = pyart.retrieve.conv_strat_yuter(grid, bkg_rad_km=50)

    assert "convsf" in dict.keys()
    assert "convsf_under" in dict.keys()
    assert "convsf_over" in dict.keys()
    assert np.all(
        dict["convsf"]["data"][25]
        == np.array(
            [
                0,
                0,
                0,
                0,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                2,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                0,
                0,
                0,
                0,
            ]
        )
    )


def test_conv_strat_yuter_noest():
    grid = pyart.testing.make_storm_grid()
    dict = pyart.retrieve.conv_strat_yuter(grid, bkg_rad_km=50, estimate_flag=False)

    assert "convsf" in dict.keys()
    assert "convsf_under" not in dict.keys()
    assert "convsf_over" not in dict.keys()


def test_hydroclass_semisupervised():
    radar = pyart.io.read(pyart.testing.NEXRAD_ARCHIVE_MSG31_FILE)
    temp_dict = pyart.config.get_metadata("temperature")
    temp = np.empty(radar.fields["reflectivity"]["data"].shape)
    temp[:] = -10.0
    temp_dict["data"] = temp
    radar.add_field("temperature", temp_dict, replace_existing=True)
    radar.fields["specific_differential_phase"] = radar.fields.pop("differential_phase")
    radar_small = radar.extract_sweeps([0])

    mass_centers = pyart.retrieve.echo_class._get_mass_centers(10e9)

    hydro_nofreq = pyart.retrieve.hydroclass_semisupervised(radar_small)
    assert "units" in hydro_nofreq.keys()
    assert "standard_name" in hydro_nofreq.keys()
    assert "long_name" in hydro_nofreq.keys()
    assert "coordinates" in hydro_nofreq.keys()

    assert hydro_nofreq["data"].dtype == "int64"
    assert_allclose(hydro_nofreq["data"][0][0:5], [6, 6, 6, 6, 6])
    assert_allclose(hydro_nofreq["data"][0][-5:], [2, 2, 2, 2, 2])

    radar_small.instrument_parameters["frequency"] = {"data": np.array([5e9])}
    hydro_freq = pyart.retrieve.hydroclass_semisupervised(radar_small)
    assert_allclose(hydro_freq["data"][0][0:5], [6, 6, 6, 6, 6])
    assert_allclose(hydro_freq["data"][0][-5:], [2, 2, 2, 2, 2])

    hydro = pyart.retrieve.hydroclass_semisupervised(
        radar_small, mass_centers=mass_centers
    )
    assert_allclose(hydro["data"][0][0:5], [6, 6, 6, 6, 6])
    assert_allclose(hydro["data"][0][-5:], [2, 2, 2, 2, 2])


def test_data_limits_table():
    dlimits_dict = pyart.retrieve.echo_class._data_limits_table()
    test_dict = {
        "Zh": (60.0, -10.0),
        "ZDR": (5.0, -5.0),
        "KDP": (7.0, -10.0),
        "RhoHV": (-5.23, -50.0),
    }

    assert isinstance(dlimits_dict, dict)
    assert dlimits_dict == test_dict


def test_get_freq():
    freq_band_s = pyart.retrieve.get_freq_band(3e9)
    assert freq_band_s == "S"

    freq_band_c = pyart.retrieve.get_freq_band(6e9)
    assert freq_band_c == "C"

    freq_band_x = pyart.retrieve.get_freq_band(10e9)
    assert freq_band_x == "X"

    freq_band_bad = pyart.retrieve.get_freq_band(10)
    assert freq_band_bad is None


def test_mass_centers_dict():
    mass_centers_dict = pyart.retrieve.echo_class._mass_centers_table()

    assert type(mass_centers_dict)
    assert "C" in mass_centers_dict.keys()
    assert "X" in mass_centers_dict.keys()

    assert_allclose(
        mass_centers_dict["C"][0],
        [1.35829e01, 4.06300e-01, 4.97000e-02, 9.86800e-01, 1.33030e03],
        atol=1e-7,
    )
    assert_allclose(
        mass_centers_dict["C"][8],
        [5.06186e01, -6.49000e-02, 9.46000e-02, 9.90400e-01, 1.17990e03],
        atol=1e-7,
    )

    assert_allclose(
        mass_centers_dict["X"][0],
        [1.9077e01, 4.1390e-01, 9.9000e-03, 9.8410e-01, 1.0617e03],
        atol=1e-7,
    )

    assert_allclose(
        mass_centers_dict["X"][8],
        [4.42216e01, -3.41900e-01, 6.87000e-02, 9.68300e-01, 1.27270e03],
        atol=1e-7,
    )


def test_mass_centers():
    mass_centers_good_freq = pyart.retrieve.echo_class._get_mass_centers(10e9)
    assert_allclose(
        mass_centers_good_freq[0],
        [1.90770e01, 4.13900e-01, 9.90000e-03, 9.84100e-01, 1.06170e03],
        atol=1e-7,
    )

    low_bad_freq = pyart.retrieve.echo_class._get_mass_centers(3e9)
    assert_allclose(
        low_bad_freq[0],
        [1.35829e01, 4.06300e-01, 4.97000e-02, 9.86800e-01, 1.33030e03],
        atol=1e-7,
    )

    high_bad_freq = pyart.retrieve.echo_class._get_mass_centers(13e9)
    assert_allclose(
        high_bad_freq[0],
        [1.90770e01, 4.13900e-01, 9.90000e-03, 9.84100e-01, 1.06170e03],
        atol=1e-7,
    )


def test_assign_to_class():
    rhohv = np.array(([0.2, 0.2, 0.2, 0.2, 0.2], [0.2, 0.2, 0.2, 0.2, 0.2]))
    zdr = np.array(([-7.8, -7.8, -7.8, -7.8, -7.8], [-7.8, -7.8, -7.8, -7.8, -7.8]))
    kdp = np.array(
        ([180.5, 180.5, 180.5, 180.5, 180.5], [180.5, 180.5, 180.5, 180.5, 180.5])
    )
    relh = np.array(([10.0, 10.0, 10.0, 10.0, 10.0], [10.0, 10.0, 10.0, 10.0, 10.0]))
    zh = np.array(
        ([-32.0, -32.0, -32.0, -32.0, -32.0], [-32.0, -32.0, -32.0, -32.0, -32.0])
    )

    mass_centers = pyart.retrieve.echo_class._get_mass_centers(10e9)
    hydroclass, min_dist = pyart.retrieve.echo_class._assign_to_class(
        zh, zdr, kdp, rhohv, relh, mass_centers
    )

    assert_allclose(hydroclass[0], [7, 7, 7, 7, 7], atol=1e-7)
    assert_allclose(
        min_dist[0],
        [227.0343910, 227.0343910, 227.0343910, 227.0343910, 227.0343910],
        atol=1e-7,
    )


def test_standardize():
    kdp = np.array(
        ([180.5, 180.5, 180.5, 180.5, 180.5], [180.5, 180.5, 180.5, 180.5, 180.5])
    )
    field_std_kdp = pyart.retrieve.echo_class._standardize(
        kdp, "KDP", mx=180.0, mn=181.0
    )
    assert_allclose(field_std_kdp[0], [-1.0, -1.0, -1.0, -1.0, -1.0], atol=1e-6)

    relh = np.array(([10.0, 10.0, 10.0, 10.0, 10.0], [10.0, 10.0, 10.0, 10.0, 10.0]))
    field_std_relh = pyart.retrieve.echo_class._standardize(
        relh, "relH", mx=9.0, mn=5.0
    )
    assert_allclose(
        field_std_relh[0], [0.024994, 0.024994, 0.024994, 0.024994, 0.024994], atol=1e-6
    )

    rhohv = np.array(([0.2, 0.2, 0.2, 0.2, 0.2], [0.2, 0.2, 0.2, 0.2, 0.2]))
    field_std_rhohv = pyart.retrieve.echo_class._standardize(
        rhohv, "RhoHV", mx=0.1, mn=1.0
    )
    assert_allclose(field_std_rhohv[0], [-1.0, -1.0, -1.0, -1.0, -1.0], atol=1e-6)

    field_std_no_limits = pyart.retrieve.echo_class._standardize(
        rhohv, "RhoHV", mx=None, mn=None
    )
    assert_allclose(field_std_no_limits[0], [1.0, 1.0, 1.0, 1.0, 1.0], atol=1e-6)

    pytest.raises(ValueError, pyart.retrieve.echo_class._standardize, rhohv, "foo")
