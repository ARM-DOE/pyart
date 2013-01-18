"""
Utilities for converting to RSL radar objects
"""

import ctypes

from netCDF4 import num2date

import py4dd


def mdv_to_rsl(myfile, trans):
    """
    Convert a mdv object into a RSL Radar object

    Parameters
    ----------

    """
    radar = py4dd.RSL_new_radar(40)

    for mdv_fld in trans.keys():
        mdv_data = myfile.read_a_field(myfile.fields.index(mdv_fld))
        rsl_index = getattr(py4dd.fieldTypes(), trans[mdv_fld])
        print "Transfering ", mdv_fld, " to ", trans[mdv_fld], \
            " which has an RSL index of ", rsl_index
        nsweeps, nrays, ngates = mdv_data.shape
        print nsweeps, nrays, ngates
        vol = py4dd.RSL_new_volume(nsweeps)
        radar.contents.v[rsl_index] = vol
        vol.contents.h.field_type = trans[mdv_fld]
        vol.contents.h.nsweeps = nsweeps
        vol.contents.h.calibr_const = 0
        field_f, field_invf = py4dd.conversion_functions[trans[mdv_fld]]
        vol.contents.h.f = field_f
        vol.contents.h.invf = field_invf
        for i_s in range(nsweeps):
            print "Sweep: ", i_s, "of ", nsweeps
            sweep = py4dd.RSL_new_sweep(nrays)
            vol.contents.sweep[i_s] = sweep
            sweep.contents.h.field_type = trans[mdv_fld]
            sweep.contents.h.sweep_num = i_s + 1  # one-indexed
            sweep.contents.h.elev = myfile.el_deg[i_s]
            sweep.contents.h.beam_width = 1.0  # change this
            sweep.contents.h.nrays = nrays
            sweep.contents.h.f = field_f
            sweep.contents.h.invf = field_invf
            sweep_time = myfile.times['time_begin']
            if myfile.scan_type == 'rhi':
                sweep_az = [myfile.az_deg[0]] * nrays
            else:
                sweep_az = myfile.az_deg
            #print sweep_az
            sweep_beam_width = 1.0
            sweep_gate_width = myfile.field_headers[1]['grid_dx']*1000
            sweep_nyquist = myfile.radar_info['unambig_vel_mps']
            for i_r in range(nrays):
                ray = py4dd.RSL_new_ray(ngates)
                sweep.contents.ray[i_r] = ray
                ray.contents.h.nbins = ngates
                ray.contents.h.year = sweep_time.year
                ray.contents.h.month = sweep_time.month
                ray.contents.h.day = sweep_time.day
                ray.contents.h.hour = sweep_time.hour
                ray.contents.h.minute = sweep_time.minute
                ray.contents.h.sec = sweep_time.second
                ray.contents.h.azimuth = sweep_az[i_r]
                ray.contents.h.ray_num = i_r + 1  # one-indexed
                ray.contents.h.elev = myfile.el_deg[i_s]
                ray.contents.h.elev_num = i_s + 1  # one-indexed
                ray.contents.h.range_bin1 = int(
                    myfile.radar_info['start_range_km'] * 1000.0)
                ray.contents.h.gate_size = int(
                    myfile.field_headers[1]['grid_dx'] * 1000.0)
                ray.contents.h.fix_angle = myfile.el_deg[i_s]
                ray.contents.h.lat = myfile.radar_info['latitude_deg']
                ray.contents.h.lon = myfile.radar_info['longitude_deg']
                ray.contents.h.alt = int(
                    myfile.radar_info['altitude_km']*1000.0)
                ray.contents.h.beam_width = 1.0
                ray.contents.h.nyq_vel = myfile.radar_info['unambig_vel_mps']
                ray.contents.h.f = field_f
                ray.contents.h.invf = field_invf
                for i_g in range(ngates):
                    float_value = ctypes.c_float(mdv_data[i_s, i_r, i_g])
                    range_value = ray.contents.h.invf(float_value)
                    ray.contents.range[i_g] = range_value
    return radar


def radar_to_rsl(myradar, trans):
    """
    Convert a Radar object into a RSL Radar object

    Parameters
    ----------

    """

    radar = py4dd.RSL_new_radar(40)
    for radar_fld in trans.keys():
        radar_data = myradar.fields[radar_fld]['data']
        rsl_index = getattr(py4dd.fieldTypes(), trans[radar_fld])
        print "Transfering ", radar_fld, " to ", trans[radar_fld],\
            " which has an RSL index of ", rsl_index
        nsweeps, nrays, ngates = myradar.nsweeps, myradar.naz, myradar.ngates
        print nsweeps, nrays, ngates
        vol = py4dd.RSL_new_volume(nsweeps)
        radar.contents.v[rsl_index] = vol
        vol.contents.h.field_type = trans[radar_fld]
        vol.contents.h.nsweeps = nsweeps
        vol.contents.h.calibr_const = 0
        field_f, field_invf = py4dd.conversion_functions[trans[radar_fld]]
        vol.contents.h.f = field_f
        vol.contents.h.invf = field_invf
        for i_s in range(nsweeps):
            print "Sweep: ", i_s, "of ", nsweeps
            sweep_start = myradar.sweep_info['sweep_start_ray_index'][
                'data'][i_s]
            sweep_end = myradar.sweep_info['sweep_end_ray_index']['data'][i_s]
            sweep = py4dd.RSL_new_sweep(nrays)
            vol.contents.sweep[i_s] = sweep
            sweep.contents.h.field_type = trans[radar_fld]
            sweep.contents.h.sweep_num = i_s + 1  # one-indexed
            sweep.contents.h.elev = myradar.elevation['data'][sweep_start]
            sweep.contents.h.beam_width = 1.0  # change this
            sweep.contents.h.nrays = nrays
            sweep.contents.h.f = field_f
            sweep.contents.h.invf = field_invf
            sweep_time = num2date(myradar.time['data'][sweep_start],
                                  myradar.time['units'],
                                  calendar=myradar.time['calendar'])
            sweep_az = myradar.azimuth['data'][sweep_start:sweep_end]
            #print sweep_az
            sweep_beam_width = 1.0
            sweep_gate_width = (myradar.range['data'][1] -
                                myradar.range['data'][0])
            #myfile.field_headers[1]['grid_dx']*1000
            sweep_nyquist = myradar.inst_params['nyquist_velocity'][
                'data'][sweep_start]
            for i_r in range(nrays):
                ray = py4dd.RSL_new_ray(ngates)
                sweep.contents.ray[i_r] = ray
                ray.contents.h.nbins = ngates
                ray.contents.h.year = sweep_time.year
                ray.contents.h.month = sweep_time.month
                ray.contents.h.day = sweep_time.day
                ray.contents.h.hour = sweep_time.hour
                ray.contents.h.minute = sweep_time.minute
                ray.contents.h.sec = sweep_time.second
                ray.contents.h.azimuth = myradar.azimuth['data'][
                    sweep_start+i_r]
                ray.contents.h.ray_num = i_r + 1  # one-indexed
                ray.contents.h.elev = myradar.elevation['data'][
                    sweep_start+i_r]
                ray.contents.h.elev_num = i_s + 1  # one-indexed
                ray.contents.h.range_bin1 = int(myradar.range['data'][0])
                ray.contents.h.gate_size = int(sweep_gate_width)
                ray.contents.h.fix_angle = myradar.sweep_info[
                    'fixed_angle']['data'][i_s]
                ray.contents.h.lat = myradar.location['latitude']['data']
                ray.contents.h.lon = myradar.location['longitude']['data']
                ray.contents.h.alt = int(myradar.location['altitude']['data'])
                ray.contents.h.beam_width = 1.0
                ray.contents.h.nyq_vel = myradar.inst_params[
                    'nyquist_velocity']['data'][sweep_start+i_r]
                ray.contents.h.f = field_f
                ray.contents.h.invf = field_invf
                for i_g in range(ngates):
                    float_value = ctypes.c_float(
                        radar_data[sweep_start + i_r, i_g])
                    range_value = ray.contents.h.invf(float_value)
                    ray.contents.range[i_g] = range_value
    return radar
