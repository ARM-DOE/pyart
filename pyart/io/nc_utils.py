"""
pyart.io.nc_utils
=================

Utilities for saving mapped and radar co-ordinate radar data.

.. autosummary::
    :toctree: generated/

    is_moment
    fix_variables_csapr
    append_global_metatdata_csapr
    append_global_metatdata_xsapr
    append_coords
    fix_time
    dt_to_dict
    save_netcdf_cube
    noncf_append_global_metatdata
    copy_meta
    copy_axis
    save_pyGrid
    save_mdv_ncf

"""

import sys
import socket
import getpass
import datetime as dt

import netCDF4
import numpy as np

from common import COMMON2STANDARD, dms_to_d


def is_moment(varname, moment_fixes):
    moments = moment_fixes.keys()
    return True in [foo in varname for foo in moments]


def fix_variables_csapr(cgfile):
    debug = True
    print "go"
    csapr_names = COMMON2STANDARD
    moment_fixes = {
        'DBZ_F': {
            'units': 'dBZ',
            'standard_name': 'equivalent_reflectivity_factor',
            'long_name': 'equivalent_reflectivity_factor',
            'valid_max': 80.0,
            'valid_min': -45.0},
        'VEL_F': {
            'units': 'm/s',
            'standard_name': (
                'radial_velocity_of_scatterers_away_from_instrument'),
            'long_name': (
                'radial_velocity_of_scatterers_away_from_instrument'),
            'valid_max': 95.0,
            'valid_min': -95.0},
        'KDP_F': {
            'units': 'degrees/km',
            'standard_name': 'specific_differential_phase_hv',
            'long_name': 'specific_differential_phase_hv',
            'valid_max': 20.0,
            'valid_min': -10.0},
        'ZDR_F': {
            'units': 'dB',
            'standard_name': 'log_differential_reflectivity_hv',
            'long_name': 'log_differential_reflectivity_hv',
            'valid_max': 8.0,
            'valid_min': -6.0},
        'RHOHV_F': {
            'units': 'ratio',
            'standard_name': 'cross_correlation_ratio_hv',
            'long_name': 'cross_correlation_ratio_hv',
            'valid_max': 1.0,
            'valid_min': 0.0},
        'NCP_F': {
            'units': 'ratio',
            'standard_name': 'signal_quality',
            'long_name': 'signal_quality',
            'valid_max': 1.0,
            'valid_min': 0.0,
            'comment': 'Also know as Normalized Coherent Power'},
        'WIDTH_F': {
            'units': 'm/s',
            'standard_name': 'spectrum_width',
            'long_name': 'spectrum_width',
            'valid_max': 45.0,
            'valid_min': 0.0},
        'PHIDP_F': {
            'units': 'degrees',
            'standard_name': 'differential_phase_hv',
            'long_name': 'differential_phase_hv',
            'valid_max': 180.0,
            'valid_min': -180.0},
        'VEL_COR': {
            'units': 'm/s',
            'standard_name': (
                'radial_velocity_of_scatterers_away_from_instrument'),
            'long_name': (
                'radial_velocity_of_scatterers_away_from_instrument'),
            'valid_max': 45.0,
            'valid_min': -45.0},
        'PHIDP_UNF': {
            'units': 'degrees',
            'standard_name': 'differential_phase_hv',
            'long_name': 'differential_phase_hv',
            'valid_max': 480.0,
            'valid_min': 0.0},
        'DBZ_AC': {
            'units': 'dBZ',
            'standard_name': 'equivalent_reflectivity_factor',
            'long_name': 'equivalent_reflectivity_factor',
            'valid_max': 80.0,
            'valid_min': -45.0},
        'KDP_SOB': {
            'units': 'degrees/km',
            'standard_name': 'specific_differential_phase_hv',
            'long_name': 'specific_differential_phase_hv',
            'valid_max': 20.0,
            'valid_min': -1.0}
    }
    for long_variable_name in list(set(cgfile.variables.keys()) & set(csapr_names.values())):
        variable_name = csapr_names.keys()[csapr_names.values().index(
            long_variable_name)]
        if debug:
            print "doing ", variable_name
        if is_moment(variable_name, moment_fixes):
            print "doing ", variable_name
            for attr in moment_fixes[variable_name].keys():
                setattr(cgfile.variables[long_variable_name], attr,
                        moment_fixes[variable_name][attr])


def append_global_metatdata_csapr(cgfile):
    cgfile.Conventions = 'CF 1.5'
    cgfile.title = 'Radar moments mapped to a Cartesian grid'
    cgfile.history = 'Dealiasing done with the U Washington 4DD code, attenuation correction using the ZPHI method in high (45+dBz) reflectivities, phidp only elsewhere'
    cgfile.source = 'Mapped moments from the C Band Scanning ARM Radar. '
    cgfile.institution = 'Atmospheric Radiation Measurement Climate Facility, United States Department of Energy'
    cgfile.comment = """Evaluation mapped moments from the C-SAPR radar. Please note that this data is under
    active development and there may be errors or issues in the data. Please report any issues, errors or suggestions
    to the precipitation radar data translator scollis@anl.gov"""
    cgfile.site = 'SGP'


def append_global_metatdata_xsapr(cgfile):
    cgfile.Conventions = 'CF 1.5'
    cgfile.title = 'Radar moments mapped to a Cartesian grid'
    cgfile.history = 'Dealiasing done with the U Washington 4DD code'
    cgfile.source = 'Mapped moments from the X Band Scanning ARM Radar. '
    cgfile.institution = 'Atmospheric Radiation Measurement Climate Facility, United States Department of Energy'
    cgfile.comment = """Evaluation mapped moments from the X-SAPR radar. Please note that this data is under
    active development and there may be errors or issues in the data. Please report any issues, errors or suggestions
    to the precipitation radar data translator scollis@anl.gov"""
    cgfile.site = 'SGP'


def append_coords(cgfile, xar, yar, zar, radar_loc):
    Re = 6371.0 * 1000.0
    rad_at_radar = Re * np.sin(np.pi / 2.0 -
                               np.abs(radar_loc[0] * np.pi / 180.0))
    lons = radar_loc[1] + 360.0 * xar / (rad_at_radar * 2.0 * np.pi)
    lats = radar_loc[0] + 360.0 * yar / (Re * 2.0 * np.pi)
    levs = zar
    levvar = cgfile.createVariable('z_disp', 'float32', ('nz',))
    yvar = cgfile.createVariable('y_disp', 'float32', ('ny',))
    xvar = cgfile.createVariable('x_disp', 'float32', ('nx',))
    latvar = cgfile.createVariable('lat', 'float32', ('nx',))
    lonvar = cgfile.createVariable('lon', 'float32', ('ny',))
    print lats.shape
    latvar[:] = lats
    print lons.shape
    lonvar[:] = lons
    print xar.shape
    xvar[:] = xar
    print yar.shape
    yvar[:] = yar
    print levs.shape
    levvar[:] = levs
    xvar.axis = "X"
    xvar.long_name = "x-coordinate in Cartesian system"
    xvar.standard_name = xvar.long_name
    xvar.units = "m"
    xvar.comment = "X Displacement from the central facility"
    yvar.axis = "Y"
    yvar.comment = "Y Displacement from the central facility"
    yvar.long_name = "y-coordinate in Cartesian system"
    yvar.standard_name = yvar.long_name
    yvar.units = "m"
    levvar.long_name = "height"
    levvar.standard_name = levvar.long_name
    levvar.units = "meter"
    lonvar.long_name = "longitude"
    lonvar.standard_name = lonvar.long_name
    lonvar.units = "degrees_east"
    latvar.long_name = "latitude"
    latvar.standard_name = latvar.long_name
    latvar.units = "degrees_north"
    levvar.positive = "up"


def fix_time(cgfile, time):
    times = cgfile.createVariable('time', 'double', ('time', ))
    times.units = 'seconds since 1970-01-01 00:00:00.0'
    times.calendar = 'gregorian'
    times.standard_name = 'time'
    print "time here is"
    print time
    print netCDF4.date2num(time, units=times.units, calendar=times.calendar)
    times[0] = netCDF4.date2num(time, units=times.units,
                                calendar=times.calendar)
    print times


def dt_to_dict(dt, **kwargs):
    pref = kwargs.get('pref', '')
    return dict([(pref+key, getattr(dt, key)) for key in
                ['year', 'month', 'day', 'hour', 'minute', 'second']])


def save_netcdf_cube(grids, output_file, myfile, xr, yr, zr, origin, **kwargs):
    meta = kwargs.get('meta', 'csapr_ac')
    print "opening ", output_file
    ofile = netCDF4.Dataset(output_file, 'w', format='NETCDF3_CLASSIC')
    nz, ny, nx = grids[grids.keys()[0]].shape
    ofile.createDimension('x', nx)
    ofile.createDimension('y', ny)
    ofile.createDimension('z', nz)
    ofile.createDimension('time', 1)
    csapr_names = COMMON2STANDARD
    available_pars = list(set(csapr_names.keys()) & set(grids.keys()))
    print available_pars
    vvars = [ofile.createVariable(csapr_names[par], np.float,
             ('time', 'z', 'y', 'x'), fill_value=-999) for par in
             available_pars]
    for i in range(len(available_pars)):
        vvars[i][0, :, :, :] = np.ma.masked_array(
            grids[available_pars[i]], np.isnan(grids[available_pars[i]]))
    fix_variables_csapr(ofile)
    if meta == 'csapr_ac':
        append_global_metatdata_csapr(ofile)
    elif meta == 'xsapr':
        append_global_metatdata_xsapr(ofile)
    xar = np.linspace(xr[0], xr[1], nx)
    yar = np.linspace(yr[0], yr[1], ny)
    zar = np.linspace(zr[0], zr[1], nz)
    if "radar_info" in dir(myfile):  # this is a mdv object
        radar_loc = [myfile.radar_info['latitude_deg'],
                     myfile.radar_info['longitude_deg']]
        fix_time(ofile, myfile.times['time_begin'])
    else:
        # this is an RSL object
        radar_loc = [dms_to_d((myfile.contents.h.latd,
                     myfile.contents.h.latm, myfile.contents.h.lats)),
                     dms_to_d((myfile.contents.h.lond, myfile.contents.h.lonm,
                               myfile.contents.h.lons))]
        fix_time(ofile, myfile.contents.datetime)
    append_coords(ofile, xar, yar, zar, origin)
    ofile.close()


def noncf_append_global_metatdata(cgfile):
    cgfile.Conventions = 'None, this is not a cfradial file'
    cgfile.title = 'Radar moments converted to a NetCDF cube (nsweeps, nrays, ngates)'
    cgfile.history = 'Saved by save_mdv_ncf'
    cgfile.source = 'Moments from the C Band Scanning ARM Radar. '
    cgfile.institution = 'Atmospheric Radiation Measurement Climate Facility, United States Department of Energy'
    cgfile.comment = """Moments from the C-SAPR radar. Please note that this data is under
    active development and there may be errors or issues in the data. Please report any issues, errors or suggestions
    to the precipitation radar data translator scollis@anl.gov"""
    cgfile.site = 'SGP'


def copy_meta(myfile, ofile):
    for key in myfile.radar_info.keys():
        setattr(ofile, key, myfile.radar_info[key])


def copy_axis(myfile, cgfile):
    xvar = cgfile.createVariable('x', 'float32',
                                 ('time', 'nsweeps', 'nrays', 'ngates'))
    yvar = cgfile.createVariable('y', 'float32',
                                 ('time', 'nsweeps', 'nrays', 'ngates'))
    zvar = cgfile.createVariable('z', 'float32',
                                 ('time', 'nsweeps', 'nrays', 'ngates'))
    xvar[0, :, :, :] = myfile.carts['x']
    yvar[0, :, :, :] = myfile.carts['y']
    zvar[0, :, :, :] = myfile.carts['z']
    azvar = cgfile.createVariable('azimtuth', 'float32', ('nrays',))
    rnvar = cgfile.createVariable('range', 'float32', ('ngates',))
    elvar = cgfile.createVariable('elevation', 'float32', ('nsweeps',))
    azvar[:] = myfile.az_deg
    rnvar[:] = myfile.range_km * 1000.0
    elvar[:] = myfile.el_deg
    xvar.axis = "X"
    xvar.long_name = "x-coordinate in Cartesian system"
    xvar.units = "m"
    yvar.axis = "Y"
    yvar.long_name = "y-coordinate in Cartesian system"
    yvar.units = "m"
    zvar.long_name = "height"
    zvar.units = "meter"
    zvar.positive = "up"
    azvar.long_name = "Azimuth of antenna"
    azvar.units = "degrees"
    elvar.long_name = "Elevation of antenna"
    elvar.units = "degrees"
    rnvar.long_name = "Range from antenna"
    rnvar.units = "meter"


def save_pyGrid(ncfobj, pygrid, **kwargs):
    """Saves a pyGrid object to a CF and ARM standard netcdf file

    usage: save_pyGrid(netCDF4_object, pyGrid_object)

    """

    #first create the time dimension

    ncfobj.createDimension('time', None)

    #Axes dimensions
    #grab the dimensions from the first moment field

    nz, ny, nx = pygrid.fields[pygrid.fields.keys()[0]]['data'].shape
    ncfobj.createDimension('nz', nz)
    ncfobj.createDimension('ny', ny)
    ncfobj.createDimension('nx', nx)

    #axes
    #Generate the variables

    vvars = [ncfobj.createVariable(
        key, np.float, ('time', 'nz', 'ny', 'nx'),
        fill_value=pygrid.fields[key]['_FillValue'],
        zlib=kwargs.get('zlib', False)) for key in
        pygrid.fields.keys()]

    # loop and populate attributes

    for field in pygrid.fields.keys():
        for meta in pygrid.fields[field].keys():
            if meta != 'data' and meta != '_FillValue':
                setattr(ncfobj.variables[field], meta,
                        pygrid.fields[field][meta])
    akeys = pygrid.axes.keys()
    akeys.sort()  # makes sure time comes first

    #For ARM compliance we want alt, lat and lon to be at the end

    for mkey in ['lat', 'lon', 'alt']:
        try:
            akeys.remove(mkey)
            akeys.append(mkey)
        except ValueError:
            print(mkey, " not existing")

    dims_lookup = {'time': 'time', 'x_disp': 'nx', 'y_disp': 'ny',
                   'z_disp': 'nz', 'time_end': 'time', 'time_start': 'time',
                   'lat': 'time', 'lon': 'time', 'alt': 'time'}
    avars = [ncfobj.createVariable(key, np.float, (dims_lookup[key], ))
             for key in akeys]

    # loop and populate attributes

    for axis in akeys:
        metakeys = pygrid.axes[axis].keys()

        #again, reorder to meet ARM standards..

        for mkey in ['units', 'long_name']:
            try:
                metakeys.remove(mkey)
                metakeys.insert(0, mkey)
            except ValueError:
                print(mkey, " not existing")
        for meta in metakeys:
            if meta != 'data':
                setattr(ncfobj.variables[axis],
                        meta, pygrid.axes[axis][meta])

    # global metadata
    if 'Conventions' in pygrid.metadata.keys():
        ncfobj.Conventions = pygrid.metadata['conventions']
    else:
        ncfobj.Conventions = 'CF-1.5'
    if 'process_version' in pygrid.metadata.keys():
        ncfobj.process_version = pygrid.metadata['process_version']
    for meta in pygrid.metadata.keys():
        if meta != 'history' or meta != 'process_version':
            setattr(ncfobj, meta, pygrid.metadata[meta])
    ncfobj.history = pygrid.metadata['history']

    #now populate data.. we leave this until last to speed up..

    for i in range(len(pygrid.fields.keys())):
        vvars[i][0, :, :, :] = pygrid.fields[
            pygrid.fields.keys()[i]]['data'][:, :, :]

    for i in range(len(akeys)):
        if 'shape' in dir(pygrid.axes[akeys[i]]['data']):
            avars[i][:] = pygrid.axes[akeys[i]]['data']
            print akeys[i], "is array"
        else:
            avars[i][:] = np.array([pygrid.axes[akeys[i]]['data']])
            print np.array([pygrid.axes[akeys[i]]['data']])
            print akeys[i], "is not array"


def save_mdv_ncf(myfile, output_file, parms, **kwargs):
    debug = kwargs.get('debug', False)
    if debug:
        print "opening ", output_file
    ofile = netCDF4.Dataset(output_file, 'w', format='NETCDF3_CLASSIC')
    if debug:
        print 'Appending global metadata'
    noncf_append_global_metatdata(ofile)
    if debug:
        print 'Appending radar metadata'
    copy_meta(myfile, ofile)
    nsweeps, nrays, ngates = (myfile.read_a_field(0)).shape
    if debug:
        print "Creating dimensions"
    ofile.createDimension('nsweeps', nsweeps)
    ofile.createDimension('nrays', nrays)
    ofile.createDimension('ngates', ngates)
    ofile.createDimension('time', 1)
    if debug:
        print "Creating axes"
    copy_axis(myfile, ofile)
    if debug:
        print "Creating variables"
    csapr_names = COMMON2STANDARD
    available_pars = list(set(csapr_names.keys()) & set(parms))
    vvars = [ofile.createVariable(csapr_names[par], np.float,
             ('time', 'nsweeps', 'nrays', 'ngates'), fill_value=-999) for
             par in available_pars]
    if debug:
        print 'fixing variables'
    fix_variables_csapr(ofile)
    if debug:
        print "Filling in data"
    for i in range(len(parms)):
        if debug:
            print "Doing ", parms[i]
        mdv_data = myfile.read_a_field(myfile.fields.index(parms[i]))
        vvars[i][0, :, :, :] = np.ma.masked_array(
            mdv_data, np.isnan(mdv_data))
    if debug:
        print "fixing time"
    fix_time(ofile, myfile.times['time_begin'])
    ofile.close()
