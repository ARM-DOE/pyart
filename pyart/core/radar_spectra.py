"""
A general central radial scanning (or dwelling) instrument class.

"""

import copy
import sys

import numpy as np

try:
    import xarray
    _XARRAY_AVAILABLE = True
except ImportError:
    _XARRAY_AVAILABLE = False

from ..config import get_metadata
from ..lazydict import LazyLoadDict
from .transforms import antenna_vectors_to_cartesian, cartesian_to_geographic


class RadarSpectra(object):
    """
    A class for storing antenna coordinate radar spectra data.
    The structure of the Radar class is based on the CF/Radial Data file
    format. Global attributes and variables (section 4.1 and 4.3) are
    represented as a dictionary in the metadata attribute. Other required and
    optional variables are represented as dictionaries in a attribute with the
    same name as the variable in the CF/Radial standard. When a optional
    attribute not present the attribute has a value of None. The data for a
    given variable is stored in the dictionary under the 'data' key. Moment
    field data is stored as a dictionary of dictionaries in the fields
    attribute. Sub-convention variables are stored as a dictionary of
    dictionaries under the meta_group attribute.

    Refer to the attribute section for information on the parameters.

    Attributes
    ----------
    time : dict
        Time at the center of each ray.
    range : dict
        Range to the center of each gate (bin).
    fields : dict of dicts
        Spectra fields.
    metadata : dict
        Metadata describing the instrument and data.
    scan_type : str
        Type of scan, one of 'ppi', 'rhi', 'sector' or 'other'. If the scan
        volume contains multiple sweep modes this should be 'other'.
    latitude : dict
        Latitude of the instrument.
    longitude : dict
        Longitude of the instrument.
    altitude : dict
        Altitude of the instrument, above sea level.
    altitude_agl : dict or None
        Altitude of the instrument above ground level. If not provided this
        attribute is set to None, indicating this parameter not available.
    sweep_number : dict
        The number of the sweep in the volume scan, 0-based.
    sweep_mode : dict
        Sweep mode for each mode in the volume scan.
    fixed_angle : dict
        Target angle for thr sweep. Azimuth angle in RHI modes, elevation
        angle in all other modes.
    sweep_start_ray_index : dict
        Index of the first ray in each sweep relative to the start of the
        volume, 0-based.
    sweep_end_ray_index : dict
        Index of the last ray in each sweep relative to the start of the
        volume, 0-based.
    rays_per_sweep : LazyLoadDict
        Number of rays in each sweep. The data key of this attribute is
        create upon first access from the data in the sweep_start_ray_index and
        sweep_end_ray_index attributes. If the sweep locations needs to be
        modified, do this prior to accessing this attribute or use
        :py:func:`init_rays_per_sweep` to reset the attribute.
    target_scan_rate : dict or None
        Intended scan rate for each sweep. If not provided this attribute is
        set to None, indicating this parameter is not available.
    rays_are_indexed : dict or None
        Indication of whether ray angles are indexed to a regular grid in
        each sweep. If not provided this attribute is set to None, indicating
        ray angle spacing is not determined.
    ray_angle_res : dict or None
        If rays_are_indexed is not None, this provides the angular resolution
        of the grid. If not provided or available this attribute is set to
        None.
    azimuth : dict
        Azimuth of antenna, relative to true North. Azimuth angles are
        recommended to be expressed in the range of [0, 360], but other
        representations are not forbidden.
    elevation : dict
        Elevation of antenna, relative to the horizontal plane. Elevation
        angles are recommended to be expressed in the range of [-180, 180],
        but other representations are not forbidden.
    gate_x, gate_y, gate_z : LazyLoadDict
        Location of each gate in a Cartesian coordinate system assuming a
        standard atmosphere with a 4/3 Earth's radius model. The data keys of
        these attributes are create upon first access from the data in the
        range, azimuth and elevation attributes. If these attributes are
        changed use :py:func:`init_gate_x_y_z` to reset.
    gate_longitude, gate_latitude : LazyLoadDict
        Geographic location of each gate. The projection parameter(s) defined
        in the `projection` attribute are used to perform an inverse map
        projection from the Cartesian gate locations relative to the radar
        location to longitudes and latitudes. If these attributes are changed
        use :py:func:`init_gate_longitude_latitude` to reset the attributes.
    projection : dic or str
        Projection parameters defining the map projection used to transform
        from Cartesian to geographic coordinates. The default dictionary sets
        the 'proj' key to 'pyart_aeqd' indicating that the native Py-ART
        azimuthal equidistant projection is used. This can be modified to
        specify a valid pyproj.Proj projparams dictionary or string.
        The special key '_include_lon_0_lat_0' is removed when interpreting
        this dictionary. If this key is present and set to True, which is
        required when proj='pyart_aeqd', then the radar longitude and
        latitude will be added to the dictionary as 'lon_0' and 'lat_0'.
    gate_altitude : LazyLoadDict
        The altitude of each radar gate as calculated from the altitude of the
        radar and the Cartesian z location of each gate. If this attribute
        is changed use :py:func:`init_gate_altitude` to reset the attribute.
    scan_rate : dict or None
        Actual antenna scan rate. If not provided this attribute is set to
        None, indicating this parameter is not available.
    antenna_transition : dict or None
        Flag indicating if the antenna is in transition, 1 = yes, 0 = no.
        If not provided this attribute is set to None, indicating this
        parameter is not available.
    rotation : dict or None
        The rotation angle of the antenna. The angle about the aircraft
        longitudinal axis for a vertically scanning radar.
    tilt : dict or None
        The tilt angle with respect to the plane orthogonal (Z-axis) to
        aircraft longitudinal axis.
    roll : dict or None
        The roll angle of platform, for aircraft right wing down is positive.
    drift : dict or None
        Drift angle of antenna, the angle between heading and track.
    heading : dict or None
        Heading (compass) angle, clockwise from north.
    pitch : dict or None
        Pitch angle of antenna, for aircraft nose up is positive.
    georefs_applied : dict or None
        Indicates whether the variables have had georeference calculation
        applied.  Leading to Earth-centric azimuth and elevation angles.
    instrument_parameters : dict of dicts or None
        Instrument parameters, if not provided this attribute is set to None,
        indicating these parameters are not avaiable. This dictionary also
        includes variables in the radar_parameters CF/Radial subconvention.
    radar_calibration : dict of dicts or None
        Instrument calibration parameters. If not provided this attribute is
        set to None, indicating these parameters are not available
    ngates : int
        Number of gates (bins) in a ray.
    nrays : int
        Number of rays in the volume.
    nsweeps : int
        Number of sweep in the volume.

    """
    def __init__(self, time, _range, fields, metadata, scan_type,
                 latitude, longitude, altitude,
                 sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
                 sweep_end_ray_index,
                 azimuth, elevation, altitude_agl=None,
                 target_scan_rate=None, rays_are_indexed=None,
                 ray_angle_res=None,
                 scan_rate=None, antenna_transition=None,
                 instrument_parameters=None,
                 radar_calibration=None, georefs_applied=None
                ):
        if not _XARRAY_AVAILABLE:
            raise MissingOptionalDependency(
                "Xarray is required to use RadarSpectra but is "
                "not installed!")
        self.field_names = ['spectra']

        self.ds = xr.Dataset(
            data_vars={
                'spectra': (('time', 'range', 'npulses_max'), fields),
                'scan_type': scan_type,
                'latitude': latitude,
                'longitude': longitude,
                'altitude': altitude,
                'sweep_number': sweep_number,
                'sweep_mode': sweep_mode,
                'fixed_angle': fixed_angle,
                'sweep_start_ray_index': sweep_start_ray_index,
                'sweep_end_ray_index': sweep_end_ray_index,
                'azimuth': azimuth,
                'elevation': elevation},
            coords={'time': time,
                    'range': _range,
                    'npulses_max': npulses_max},
            attrs=metadata)

        self.ds.attrs['projection'] = {'proj': 'pyart_aeqd', '_include_lon_0_lat_0': True}

        # initalize attributes with lazy load dictionaries
        self.init_rays_per_sweep()
        self.init_gate_x_y_z()
        self.init_gate_longitude_latitude()
        self.init_gate_altitude()

    @property
    def fields(self):
        field_dict = {}
        for key in self.field_names:
            if key in self.ds.variables.keys():
                field_dict[key] = self.ds[key]  
        return xr.Dataset(field_dict)

    @property
    def time(self):
        return self.ds.time
    
    @property
    def range(self):
        return self.ds.range

    @property
    def npulses_max(self):
        return self.ds.npulses_max

    @property
    def latitude(self):
        return self.ds.latitude

    @property
    def longitude(self):
        return self.ds.longitude

    @property
    def altitude(self):
        return self.ds.altitude

    @property
    def fixed_angle(self):
        return self.ds.fixed_angle

    @property
    def sweep_mode(self):
        return self.ds.sweep_mode

    @property
    def sweep_number(self):
        return self.ds.sweep_number

    @property
    def scan_type(self):
        return self.ds.scan_type

    @property
    def elevation(self):
        return self.ds.elevation

    @property
    def azimuth(self):
        return self.ds.azimuth

    @property
    def sweep_start_ray_index(self):
        return self.ds.sweep_start_ray_index

    @property
    def sweep_end_ray_index(self):
        return self.ds.sweep_end_ray_index

    @property
    def rays_per_sweep(self):
        return self.ds.rays_per_sweep

    @property
    def gate_x(self):
        return self.ds.gate_x

    @property
    def gate_y(self):
        return self.ds.gate_y

    @property
    def gate_z(self):
        return self.ds.gate_z

    @property
    def gate_latitude(self):
        return self.ds.gate_latitude

    @property
    def gate_longitude(self):
        return self.ds.gate_longitude

    @property
    def gate_altitude(self):
        return self.ds.gate_altitude

    def init_rays_per_sweep(self):
        """ Initialize or reset the rays_per_sweep attribute. """
        lazydic = LazyLoadDict(get_metadata('rays_per_sweep'))
        lazydic.set_lazy('data', _rays_per_sweep_data_factory(self.ds))
        self.ds['rays_per_sweep'] = xr.DataArray(lazydic['data'], attrs=lazydic)

    def init_gate_x_y_z(self):
        """ Initialize or reset the gate_{x, y, z} attributes. """
        gate_x = LazyLoadDict(get_metadata('gate_x'))
        gate_x.set_lazy('data', _gate_data_factory(self.ds, 0))
        self.ds['gate_x'] = xr.DataArray(gate_x['data'], dims=('time', 'range'),
                                         attrs=gate_x)

        gate_y = LazyLoadDict(get_metadata('gate_y'))
        gate_y.set_lazy('data', _gate_data_factory(self.ds, 1))
        self.ds['gate_y'] = xr.DataArray(gate_y['data'], dims=('time', 'range'),
                                         attrs=gate_y)

        gate_z = LazyLoadDict(get_metadata('gate_z'))
        gate_z.set_lazy('data', _gate_data_factory(self.ds, 2))
        self.ds['gate_z'] = xr.DataArray(gate_z['data'], dims=('time', 'range'),
                                         attrs=gate_z)

    def init_gate_longitude_latitude(self):
        """
        Initialize or reset the gate_longitude and gate_latitude attributes.
        """
        gate_longitude = LazyLoadDict(get_metadata('gate_longitude'))
        gate_longitude.set_lazy('data', _gate_lon_lat_data_factory(self.ds, 0))
        self.ds['gate_longitude'] = xr.DataArray(gate_longitude['data'],
                                                 dims=('time', 'range'),
                                                 attrs=gate_longitude)

        gate_latitude = LazyLoadDict(get_metadata('gate_latitude'))
        gate_latitude.set_lazy('data', _gate_lon_lat_data_factory(self.ds, 1))
        self.ds['gate_latitude'] = xr.DataArray(gate_latitude['data'],
                                                dims=('time', 'range'),
                                                attrs=gate_latitude)

    def init_gate_altitude(self):
        """ Initialize the gate_altitude attribute. """
        gate_altitude = LazyLoadDict(get_metadata('gate_altitude'))
        gate_altitude.set_lazy('data', _gate_altitude_data_factory(self.ds))
        self.ds['gate_altitude'] = xr.DataArray(gate_altitude['data'],
                                                dims=('time', 'range'),
                                                attrs=gate_altitude)

def _rays_per_sweep_data_factory(radar):
    """ Return a function which returns the number of rays per sweep. """
    def _rays_per_sweep_data():
        """ The function which returns the number of rays per sweep. """
        return (radar.sweep_end_ray_index.values -
                radar.sweep_start_ray_index.values + 1)
    return _rays_per_sweep_data

def _gate_data_factory(radar, coordinate):
    """ Return a function which returns the Cartesian locations of gates. """
    def _gate_data():
        """ The function which returns the Cartesian locations of gates. """
        ranges = radar.range.values
        azimuths = radar.azimuth.values
        elevations = radar.elevation.values
        cartesian_coords = antenna_vectors_to_cartesian(
            ranges, azimuths, elevations, edges=False)

        # load x, y, and z data except for the coordinate in question
        if coordinate != 0:
            radar['gate_x'] = xr.DataArray(cartesian_coords[0],
                                           dims=('time', 'range'))
        if coordinate != 1:
            radar['gate_y'] = xr.DataArray(cartesian_coords[1],
                                           dims=('time', 'range'))
        if coordinate != 2:
            radar['gate_z'] = xr.DataArray(cartesian_coords[2],
                                           dims=('time', 'range'))
        return cartesian_coords[coordinate]
    return _gate_data

def _gate_lon_lat_data_factory(radar, coordinate):
    """ Return a function which returns the geographic locations of gates. """
    def _gate_lon_lat_data():
        """ The function which returns the geographic locations gates. """
        x = radar.gate_x.values
        y = radar.gate_y.values
        projparams = radar.projection.copy()
        if projparams.pop('_include_lon_0_lat_0', False):
            projparams['lon_0'] = radar.longitude.values
            projparams['lat_0'] = radar.latitude.values
        geographic_coords = cartesian_to_geographic(x, y, projparams)
        # set the other geographic coordinate
        if coordinate == 0:
            radar['gate_latitude'] = xr.DataArray(geographic_coords[1],
                                                  dims=('time', 'range'))
        else:
            radar['gate_longitude'] = xr.DataArray(geographic_coords[0],
                                                   dims=('time', 'range'))
        return geographic_coords[coordinate]
    return _gate_lon_lat_data

def _gate_altitude_data_factory(radar):
    """ Return a function which returns the gate altitudes. """
    def _gate_altitude_data():
        """ The function which returns the gate altitudes. """
        try:
            return radar.altitude.values + radar.gate_z.values
        except ValueError:
            return np.mean(radar.altitude.values) + radar.gate_z.values
    return _gate_altitude_data
