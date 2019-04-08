"""
pyart.io.write_grid_geotiff
===========================

Write a Py-ART Grid object to a GeoTIFF file.

.. autosummary::
    :toctree: generated/

    write_grid_geotiff
    _get_rgb_values
    _create_sld

"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import shutil
from ..exceptions import MissingOptionalDependency
try:
    from osgeo import gdal
    IMPORT_FLAG = True
except ImportError:
    IMPORT_FLAG = False


def write_grid_geotiff(grid, filename, field, rgb=False, level=None,
                       cmap='viridis', vmin=0, vmax=75, color_levels=None,
                       warp=False, sld=False, use_doublequotes=False):
    """
    Write a Py-ART Grid object to a GeoTIFF file.

    The GeoTIFF can be the standard Azimuthal Equidistant projection used
    in Py-ART, or a lat/lon projection on a WGS84 sphere. The latter is
    typically more usable in web mapping applications. The GeoTIFF can
    contain a single float-point raster band, or three RGB byte raster bands.
    The former will require an SLD file for colorful display using standard
    GIS or web mapping software, while the latter will show colors
    "out-of-the-box" but lack actual data values. The function also can
    output an SLD file based on the user-specified inputs. User can specify
    the 2D vertical level to be output. If this is not specified, a 2D
    composite is created. User also can specify the field to output.

    This function requires GDAL Python libraries to be installed. These are
    available via conda; e.g., 'conda install gdal'

    Parameters
    ----------
    grid : pyart.core.Grid object
        Grid object to write to file.
    filename : str
        Filename for the GeoTIFF.
    field : str
        Field name to output to file.

    Other Parameters
    ----------------
    rbg : bool, optional
        True - Output 3-band RGB GeoTIFF

        False - Output single-channel, float-valued GeoTIFF. For display,
                likely will need an SLD file to provide a color table.

    level : int or None, optional
        Index for z-axis plane to output. None gives composite values
        (i.e., max in each vertical column).
    cmap : str or matplotlib.colors.Colormap object, optional
        Colormap to use for RGB output or SLD file.
    vmin : int or float, optional
        Minimum value to color for RGB output or SLD file.
    vmax : int or float, optional
        Maximum value to color for RGB output or SLD file.
    color_levels : int or None, optional
        Number of color levels in cmap. Useful for categorical colormaps
        with steps << 255 (e.g., hydrometeor ID).
    warp : bool, optional
        True - Use gdalwarp (called from command line using os.system)
               to warp to a lat/lon WGS84 grid.

        False - No warping will be performed. Output will be Az. Equidistant.

    sld : bool, optional
        True - Create a Style Layer Descriptor file (SLD) mapped to vmin/vmax
               and cmap. File is named same as output TIFF, except for .sld
               extension.

        False - Don't do this.

    use_doublequotes : bool, optional
        True - Use double quotes in the gdalwarp call (requires warp=True),
               which may help if that command is producing and error like:
               'Translating source or target SRS failed'

        False - Use single quotes instead

    """
    if not IMPORT_FLAG:
        raise MissingOptionalDependency(
            'GDAL not detected, GeoTIFF output failure!')

    if field not in grid.fields.keys():
        raise KeyError('Failed -', field, 'field not found in Grid object.')

    # Determine whether filename template already contains a suffix
    # If not, append an appropriate one.
    if '.' not in filename:
        name = filename
        end = 'tif'
        ofile = name + "." + end
    else:
        ofile = filename
    nz, ny, nx = grid.fields[field]['data'].shape
    dist = max(grid.x['data'])
    rangestep = grid.x['data'][1] - grid.x['data'][2]
    lat = grid.origin_latitude['data'][0]
    lon = grid.origin_longitude['data'][0]
    # Check if masked array; if so, fill missing data
    filled = np.ma.filled(grid.fields[field]['data'], fill_value=-32768)
    if level is None:
        data = np.amax(filled, 0)
    else:
        data = filled[level]
    data = data.astype(float)
    data[data == -32768] = np.nan

    iproj = 'PROJCS["unnamed",GEOGCS["WGS 84",DATUM["unknown",' + \
        'SPHEROID["WGS84",6378137,298.257223563]],' + \
        'PRIMEM["Greenwich",0],' + \
        'UNIT["degree",0.0174532925199433]],' + \
        'PROJECTION["Azimuthal_Equidistant"],' + \
        'PARAMETER["latitude_of_center",' + str(lat) + '],' + \
        'PARAMETER["longitude_of_center",' + str(lon) + '],' + \
        'PARAMETER["false_easting",0],' + \
        'PARAMETER["false_northing",0],' + \
        'UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'
    out_driver = gdal.GetDriverByName("GTiff")

    # Output dataset depends on rgb flag
    if not rgb:
        # Single-channel, floating-point output
        dst_options = ['COMPRESS=LZW', 'ALPHA=YES']
        dst_ds = out_driver.Create(
            ofile, data.shape[1], data.shape[0], 1,
            gdal.GDT_Float32, dst_options)
    else:
        # Assign data RGB levels based on value relative to vmax/vmin
        rarr, garr, barr = _get_rgb_values(
            data, vmin, vmax, color_levels, cmap)
        dst_ds = out_driver.Create(ofile, data.shape[1],
                                   data.shape[0], 3, gdal.GDT_Byte)

    # Common Projection and GeoTransform
    dst_ds.SetGeoTransform([-dist, -rangestep, 0, dist, 0, rangestep])
    dst_ds.SetProjection(iproj)

    # Final output depends on rgb flag
    if not rgb:
        dst_ds.GetRasterBand(1).WriteArray(data[::-1, :])
    else:
        dst_ds.GetRasterBand(1).WriteArray(rarr[::-1, :])
        dst_ds.GetRasterBand(2).WriteArray(garr[::-1, :])
        dst_ds.GetRasterBand(3).WriteArray(barr[::-1, :])
    dst_ds.FlushCache()
    dst_ds = None

    if sld:
        _create_sld(cmap, vmin, vmax, ofile, color_levels)

    if warp:
        # Warps TIFF to lat/lon WGS84 projection that is more useful
        # for web mapping applications. Likely changes array shape.
        if use_doublequotes:
            os.system('gdalwarp -q -t_srs \"+proj=longlat +ellps=WGS84 ' +
                      '+datum=WGS84 +no_defs\" ' + ofile + ' ' +
                      ofile + '_tmp.tif')
        else:
            os.system('gdalwarp -q -t_srs \'+proj=longlat +ellps=WGS84 ' +
                      '+datum=WGS84 +no_defs\' ' + ofile + ' ' +
                      ofile + '_tmp.tif')
        shutil.move(ofile+'_tmp.tif', ofile)


def _get_rgb_values(data, vmin, vmax, color_levels, cmap):
    """
    Get RGB values for later output to GeoTIFF, given a 2D data field,
    display min/max and color table info. Missing data get numpy.nan.
    Only called if rgb is True in write_grid_geotiff.

    Parameters
    ----------
    data : numpy.ndarray object, dtype int or float
        Two-dimensional data array
    vmin : int or float
        Minimum value to color for RGB output or SLD file.
    vmax : int or float
        Maximum value to color for RGB output or SLD file.
    color_levels : int
        Number of color levels in cmap. Useful for categorical colormaps
        with steps << 255 (e.g., hydrometeor ID).
    cmap : str or matplotlib.colors.Colormap object, optional
        Colormap to use for RGB output or SLD file.

    Returns
    -------
    rarr : numpy.ndarray object, dtype int
        Red channel indices (range = 0-255)
    barr : numpy.ndarray object, dtype int
        Blue channel indices (range = 0-255)
    garr : numpy.ndarray object, dtype int
        Green channel indices (range = 0-255)

    """
    frac = (data - vmin) / np.float(vmax-vmin)
    if color_levels is None:
        color_levels = 255
    index = (frac * color_levels).ravel()
    # Out-of-bounds values will be lowest/highest colors
    index[index < 0] = 0
    index[index > 255] = 255
    rarr = []
    garr = []
    barr = []
    cmap = plt.cm.get_cmap(cmap)
    for val in index:
        if not np.isnan(val):
            ind = np.int(np.round(val))
            r, g, b, t = cmap(ind)
            rarr.append(np.int(np.round(r * 255)))
            garr.append(np.int(np.round(g * 255)))
            barr.append(np.int(np.round(b * 255)))
        else:
            rarr.append(np.nan)
            garr.append(np.nan)
            barr.append(np.nan)
    rarr = np.reshape(rarr, data.shape)
    garr = np.reshape(garr, data.shape)
    barr = np.reshape(barr, data.shape)
    return rarr, garr, barr


def _create_sld(cmap, vmin, vmax, filename, color_levels=None):
    """
    Develop a Style Layer Descriptor file given a color table and
    user-specified min/max files. Output color info to that file.
    Only called if sld is True in write_grid_geotiff.

    Parameters
    ----------
    cmap : str or matplotlib.colors.Colormap object, optional
        Colormap to use for RGB output or SLD file.
    vmin : int or float
        Minimum value to color for RGB output or SLD file.
    vmax : int or float
        Maximum value to color for RGB output or SLD file.
    filename : str
        Template for SLD filename. The suffix (presumably .tif or .tiff)
        is removed and replaced with .sld. Thus, if provided a filename
        radar_reflectivity.tif, the output SLD file will be called
        radar_reflectivity.sld.

    Other Parameters
    ----------------
    color_levels : int or None, optional
        Number of color levels in cmap. Useful for categorical colormaps
        with steps << 255 (e.g., hydrometeor ID).

    """
    cmap = plt.cm.get_cmap(cmap)
    if color_levels is None:
        color_levels = 255
    name, end = filename.split('.')
    ofile = name + '.sld'
    fileobj = open(ofile, 'w')

    header = """<?xml version="1.0" encoding="UTF-8"?>
<sld:StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:sld="http://www.opengis.net/sld" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml" version="1.0.0">
    <sld:UserLayer>
        <sld:LayerFeatureConstraints>
            <sld:FeatureTypeConstraint/>
        </sld:LayerFeatureConstraints>
        <sld:UserStyle>
            <sld:Name>""" + str(name) + """</sld:Name>
            <sld:FeatureTypeStyle>
                <sld:Name>name</sld:Name>
                <sld:FeatureTypeName>Feature</sld:FeatureTypeName>
                <sld:Rule>
                    <sld:RasterSymbolizer>
                        <sld:ColorMap>"""
    fileobj.write(header)

    for i in np.arange(color_levels + 1):
        val = i * (vmax - vmin) / (color_levels) + vmin  # color_levels + 1
        rgbt = cmap(i)
        if i == 0:
            op = 0.0
        else:
            op = 1.0
        hexval = colors.rgb2hex(rgbt[0:3])
        wstr = '\n                            <sld:ColorMapEntry color=\"' + \
            str(hexval).upper() + '\" opacity=\"' + str(op) + \
            '\" quantity=\"' + str(val) + '\"/>'
        fileobj.write(wstr)

    footer = """
                        </sld:ColorMap>
                    </sld:RasterSymbolizer>
                </sld:Rule>
            </sld:FeatureTypeStyle>
        </sld:UserStyle>
    </sld:UserLayer>
</sld:StyledLayerDescriptor>
"""
    fileobj.write(footer)
    fileobj.close()
