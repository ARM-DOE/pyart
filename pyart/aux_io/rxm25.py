"""
Routines for Ridgeline Instruments RXM-25 formatted NetCDF files.

"""

import datetime

try:
    import pytz
    _PYTZ_AVAILABLE = True
except ImportError:
    _PYTZ_AVAILABLE = False

import netCDF4
import numpy as np

import pyart
from ..config import get_metadata
from ..core.radar import Radar
from ..exceptions import MissingOptionalDependency
from ..testing import make_empty_ppi_radar

def read_rxm25(filename, cfradial_outfile=None, heading=None):
    """
    Read in Ridgeline Instruments RXM-25 formatted NetCDF data.

    Parameters
    ----------
    filename : str
        Name of Ridgeline Instruments (RLI) RXM-25 formatted NetCDF file 
        from which to read data.
    cfradial_outfile : str, optional
        If file is to be converted to CF-Radial format, specify the output 
        filename here.
    heading : float, optional
        If a heading offset exists, enter it here (in degrees).

    Returns
    -------
    radar : Radar
        Radar object.

    """
    if not _PYTZ_AVAILABLE:
        raise MissingOptionalDependency(
            "Pytz is required to use read_rxm25 but is " +
            "not installed")

    data  = netCDF4.Dataset(filename, 'r')

    ngates = data.dimensions['Gate'].size
    rays_per_sweep = data.dimensions['Radial'].size
    radar = make_empty_ppi_radar(ngates, rays_per_sweep, 1)

    # Time needs to be converted from nss1970 to nss1989 and added to
    # Radar object.
    nineteen89 = datetime.datetime(1989, 1, 1, 0, 0, 1, tzinfo = pytz.utc)
    baseTime = np.array(
        [datetime.datetime.fromtimestamp(
            t, tz=pytz.UTC) for t in data.variables['Time'][:]])
    radar.time['data'] = np.array(
        [t.total_seconds() for t in baseTime - nineteen89])

    if heading is not None:
        radar.heading = heading
        radar.azimuth['data'] = np.mod(
            data['Azimuth'][:] - radar.heading, 360.)
    else:
        radar.azimuth['data'] = data['Azimuth'][:]

    radar.longitude['data'] = np.array([data.Longitude], dtype='float64')
    radar.latitude['data'] = np.array([data.Latitude], dtype='float64')
    radar.elevation['data'] = data['Elevation'][:]
    radar.altitude['data'] = np.array([data.Height], dtype='float64')

    fixed_agl_data = np.empty((1, ), dtype='float32')
    fixed_agl_data[:] = np.mean(
        radar.elevation['data'][:rays_per_sweep])

    radar.fixed_angle['data'] = fixed_agl_data

    radar.range['data'] = np.linspace(
        data['StartRange'][0]/1000,
        (ngates - 1)*data['GateWidth'][0]/1000 + data['StartRange'][0]/1000,
        ngates)

    ref = data['Reflectivity'][:]
    norm_pow = data['NormalizedCoherentPower'][:]
    spec_w = data['SpectralWidth'][:]
    vel = data['Velocity'][:]
    corr_ref = data['CorrectedReflectivity'][:]
    diff_ref = data['DifferentialReflectivity'][:]
    diff_phase = data['DifferentialPhase'][:]
    spec_phase = data['SpecificPhase'][:]
    corr_diff_ref = data['CorrectedDifferentialReflectivity'][:]
    sig_noise = data['SignalToNoiseRatio'][:]
    rain_rate = data['RainfallRate'][:]
    cross_ra = data['CrossPolCorrelation'][:]

    fields = {
        'reflectivity': get_metadata('reflectivity'),
        'normalized_coherent_power': get_metadata('normalized_coherent_power'),
        'spectral_width': get_metadata('spectral_width'),
        'velocity': get_metadata('velocity'),
        'corrected_reflectivity': get_metadata('correct_reflectivity'),
        'differential_reflectivity': get_metadata('differential_reflectivity'),
        'differential_phase': get_metadata('differential_phase'),
        'specific_differential_phase': get_metadata('specific_differential_phase'),
        'corrected_differential_reflectivity': get_metadata('corrected_differential_reflectivity'),
        'signal_to_noise_ratio': get_metadata('signal_to_noise_ratio'),
        'rain_rate': get_metadata('rain_rate'),
        'cross_correlation_ratio': get_metadata('cross_correlation_ratio')}
    
    radar.fields = fields
    radar.fields['reflectivity']['data'] = ref
    radar.fields['normalized_coherent_power']['data'] = norm_pow
    radar.fields['spectral_width']['data'] = spec_w
    radar.fields['velocity']['data'] = vel
    radar.fields['corrected_reflectivity']['data'] = corr_ref
    radar.fields['differential_reflectivity']['data'] = diff_ref
    radar.fields['differential_phase']['data'] = diff_phase
    radar.fields['specific_differential_phase']['data'] = spec_phase
    radar.fields['corrected_differential_reflectivity']['data'] = corr_diff_ref
    radar.fields['signal_to_noise_ratio']['data'] = sig_noise
    radar.fields['rain_rate']['data'] = rain_rate
    radar.fields['cross_correlation_ratio']['data'] = cross_ra

    radar.metadata['instrument_name'] = 'RXM-25'

    if cfradial_outfile is not None:
        pyart.io.write_cfradial(
            cfradial_outfile, radar, arm_time_variables=True)

    return radar
