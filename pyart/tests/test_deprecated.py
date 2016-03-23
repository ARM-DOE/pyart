""" Unit tests for deprecated Py-ART function names. """

import types

import pyart


def test_io_common_radar_coords_to_cart_exists():
    assert isinstance(pyart.io.common.radar_coords_to_cart, types.FunctionType)


def test_graph_common_radar_coords_to_cart_exists():
    assert isinstance(
        pyart.graph.common.radar_coords_to_cart, types.FunctionType)


def test_graph_common_sweep_coords_to_cart_exists():
    assert isinstance(
        pyart.graph.common.sweep_coords_to_cart, types.FunctionType)
