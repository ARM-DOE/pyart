""" Unit Tests for Py-ART's graph/cm.py module. """


import matplotlib

from pyart.graph import cm_colorblind


def test_colormaps_exist():
    assert isinstance(cm_colorblind.HomeyerRainbow, matplotlib.colors.Colormap)
    assert isinstance(cm_colorblind.HomeyerRainbow, matplotlib.colors.Colormap)


def test_colormaps_registered():
    cmap = matplotlib.cm.get_cmap('pyart_HomeyerRainbow')
    assert isinstance(cmap, matplotlib.colors.Colormap)

    cmap = matplotlib.cm.get_cmap('pyart_HomeyerRainbow_r')
    assert isinstance(cmap, matplotlib.colors.Colormap)
