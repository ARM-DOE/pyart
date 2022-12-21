""" Unit Tests for Py-ART's graph/cm_colorblind.py module. """


import matplotlib

from pyart.graph import cm_colorblind


def test_colormaps_exist():
    assert isinstance(cm_colorblind.HomeyerRainbow, matplotlib.colors.Colormap)
    assert isinstance(cm_colorblind.HomeyerRainbow_r, matplotlib.colors.Colormap)
    assert isinstance(cm_colorblind.balance, matplotlib.colors.Colormap)
    assert isinstance(cm_colorblind.balance_r, matplotlib.colors.Colormap)
    assert isinstance(cm_colorblind.ChaseSpectral, matplotlib.colors.Colormap)
    assert isinstance(cm_colorblind.ChaseSpectral_r, matplotlib.colors.Colormap)
    assert isinstance(cm_colorblind.SpectralExtended, matplotlib.colors.Colormap)
    assert isinstance(cm_colorblind.SpectralExtended_r, matplotlib.colors.Colormap)


def test_colormaps_registered():
    cmap = matplotlib.colormaps.get_cmap("pyart_HomeyerRainbow")
    assert isinstance(cmap, matplotlib.colors.Colormap)

    cmap = matplotlib.colormaps.get_cmap("pyart_HomeyerRainbow_r")
    assert isinstance(cmap, matplotlib.colors.Colormap)

    cmap = matplotlib.colormaps.get_cmap("pyart_balance")
    assert isinstance(cmap, matplotlib.colors.Colormap)

    cmap = matplotlib.colormaps.get_cmap("pyart_balance_r")
    assert isinstance(cmap, matplotlib.colors.Colormap)

    cmap = matplotlib.colormaps.get_cmap("pyart_ChaseSpectral")
    assert isinstance(cmap, matplotlib.colors.Colormap)

    cmap = matplotlib.colormaps.get_cmap("pyart_ChaseSpectral_r")
    assert isinstance(cmap, matplotlib.colors.Colormap)

    cmap = matplotlib.colormaps.get_cmap("pyart_SpectralExtended")
    assert isinstance(cmap, matplotlib.colors.Colormap)

    cmap = matplotlib.colormaps.get_cmap("pyart_SpectralExtended_r")
    assert isinstance(cmap, matplotlib.colors.Colormap)
