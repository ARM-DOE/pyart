"""
A general central radial scanning (or dwelling) instrument class.

"""

import copy
import sys
import warnings

import numpy as np

try:
    import xarray as xr
    _XARRAY_AVAILABLE = True
except ImportError:
    _XARRAY_AVAILABLE = False

from ..config import get_metadata
from ..exceptions import MissingOptionalDependency
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
                 azimuth, elevation, npulses_max, velocity_bins,
                 altitude_agl=None,
                 target_scan_rate=None, rays_are_indexed=None,
                 ray_angle_res=None,
                 scan_rate=None, antenna_transition=None,
                 instrument_parameters=None,
                 radar_calibration=None, georefs_applied=None
                ):
        warnings.warn("Radar Spectra object is in early development, "
                      "errors may arise, use at your own risk! ")
        if not _XARRAY_AVAILABLE:
            raise MissingOptionalDependency(
                "Xarray is required to use RadarSpectra but is "
                "not installed!")
        self.field_names = ['spectra']

        self.ds = xr.Dataset(
            data_vars={
                'spectra': (('time', 'range', 'npulses_max'), fields),
                'velocity_bins': velocity_bins,
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

        self.ds['ngates'] = len(_range.values)
        self.ds['nrays'] = len(time.values)
        self.ds['nsweeps'] = len(sweep_number.values)
        self.ds.attrs['projection'] = {'proj': 'pyart_aeqd',
                                       '_include_lon_0_lat_0': True}

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
    def velocity_bins(self):
        return self.ds.velocity_bins

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

    @property
    def ngates(self):
        return self.ds.ngates

    @property
    def nrays(self):
        return self.ds.nrays

    @property
    def nsweeps(self):
        return self.ds.nsweeps

    @property
    def projection(self):
        return self.ds.attrs['projection']

    def init_rays_per_sweep(self):
        """ Initialize or reset the rays_per_sweep attribute. """
        _rays_per_sweep_data_factory(self.ds)

    def init_gate_x_y_z(self):
        """ Initialize or reset the gate_{x, y, z} attributes. """
        _gate_data_factory(self.ds)

    def init_gate_longitude_latitude(self):
        """
        Initialize or reset the gate_longitude and gate_latitude attributes.
        """
        _gate_lon_lat_data_factory(self.ds)

    def init_gate_altitude(self):
        """ Initialize the gate_altitude attribute. """
        _gate_altitude_data_factory(self.ds)

    def _check_sweep_in_range(self, sweep):
        """ Check that a sweep number is in range. """
        if sweep < 0 or sweep >= self.nsweeps:
            raise IndexError('Sweep out of range: ', sweep)

        # get methods

    def get_start(self, sweep):
        """ Return the starting ray index for a given sweep. """
        self._check_sweep_in_range(sweep)
        return self.sweep_start_ray_index.values[sweep]

    def get_end(self, sweep):
        """ Return the ending ray for a given sweep. """
        self._check_sweep_in_range(sweep)
        return self.sweep_end_ray_index.values[sweep]

    def get_start_end(self, sweep):
        """ Return the starting and ending ray for a given sweep. """
        return self.get_start(sweep), self.get_end(sweep)

    def get_slice(self, sweep):
        """ Return a slice for selecting rays for a given sweep. """
        start, end = self.get_start_end(sweep)
        return slice(start, end+1)

    def check_field_exists(self, field_name):
        """
        Check that a field exists in the fields dictionary.
        If the field does not exist raise a KeyError.

        Parameters
        ----------
        field_name : str
            Name of field to check.

        """
        if field_name not in self.fields.keys():
            raise KeyError('Field not available: ' + field_name)

    # Iterators
    def iter_start(self):
        """ Return an iterator over the sweep start indices. """
        return (s for s in self.sweep_start_ray_index.values)

    def iter_end(self):
        """ Return an iterator over the sweep end indices. """
        return (s for s in self.sweep_end_ray_index.values)

    def iter_start_end(self):
        """ Return an iterator over the sweep start and end indices. """
        return ((s, e) for s, e in zip(self.iter_start(), self.iter_end()))

    def iter_slice(self):
        """ Return an iterator which returns sweep slice objects. """
        return (slice(s, e+1) for s, e in self.iter_start_end())

    def iter_field(self, field_name):
        """ Return an iterator which returns sweep field data. """
        self.check_field_exists(field_name)
        return (self.fields[field_name].values[s] for s in self.iter_slice())

    def iter_azimuth(self):
        """ Return an iterator which returns sweep azimuth data. """
        return (self.azimuth.values[s] for s in self.iter_slice())

    def iter_elevation(self):
        """ Return an iterator which returns sweep elevation data. """
        return (self.elevation.values[s] for s in self.iter_slice())

    def to_vpt(self):
        """ Returns a simple Radar object in VPT scan type with spectra moments
        such as reflectivity and mean velocity. """
        from ..testing import make_empty_ppi_radar
        from ..retrieve import spectra_moments
        from ..util import to_vpt

        rng_len = len(self.range.values)
        time_len = len(self.time.values)
        vpt_radar = make_empty_ppi_radar(
            ngates=rng_len, rays_per_sweep=time_len, nsweeps=1)

        fields = spectra_moments(self)

        rng_dict = get_metadata('range')
        rng_dict['data'] = self.range.values

        time_dict = get_metadata('time')
        time_dict['data'] = self.time.values

        vpt_radar.range = rng_dict
        vpt_radar.time = time_dict
        vpt_radar.fields = fields
        vpt_radar.metadata['instrument_name'] = 'KAZR'
        to_vpt(vpt_radar)
        return vpt_radar


def _rays_per_sweep_data_factory(radar):
    """ Return a function which returns the number of rays per sweep. """
    rays_per_sweep_dict = get_metadata('rays_per_sweep')
    rays_per_sweep = (
        radar.sweep_end_ray_index.values -
        radar.sweep_start_ray_index.values + 1)
    radar['rays_per_sweep'] = xr.DataArray(np.array(rays_per_sweep),
                                           attrs=rays_per_sweep_dict)


def _gate_data_factory(radar):
    """ Return a function which returns the Cartesian locations of gates. """
    ranges = radar.range.values
    azimuths = radar.azimuth.values
    elevations = radar.elevation.values
    cartesian_coords = antenna_vectors_to_cartesian(
        ranges, azimuths, elevations, edges=False)

    # load x, y, and z data except for the coordinate in question
    gate_x_dict = get_metadata('gate_x')
    radar['gate_x'] = xr.DataArray(cartesian_coords[0],
                                   dims=('time', 'range'),
                                   attrs=gate_x_dict)
    gate_y_dict = get_metadata('gate_y')
    radar['gate_y'] = xr.DataArray(cartesian_coords[1],
                                   dims=('time', 'range'),
                                   attrs=gate_y_dict)
    gate_z_dict = get_metadata('gate_z')
    radar['gate_z'] = xr.DataArray(cartesian_coords[2],
                                   dims=('time', 'range'),
                                   attrs=gate_z_dict)

def _gate_lon_lat_data_factory(radar):
    """ Return a function which returns the geographic locations of gates. """
    x = radar.gate_x.values
    y = radar.gate_y.values
    projparams = radar.projection.copy()
    if projparams.pop('_include_lon_0_lat_0', False):
        projparams['lon_0'] = radar.longitude.values
        projparams['lat_0'] = radar.latitude.values
    geographic_coords = cartesian_to_geographic(x, y, projparams)
    # set the other geographic coordinate
    gate_latitude_dict = get_metadata('gate_latitude')
    radar['gate_latitude'] = xr.DataArray(geographic_coords[1],
                                          dims=('time', 'range'),
                                          attrs=gate_latitude_dict)
    gate_longitude_dict = get_metadata('gate_longitude')
    radar['gate_longitude'] = xr.DataArray(geographic_coords[0],
                                           dims=('time', 'range'),
                                           attrs=gate_longitude_dict)

def _gate_altitude_data_factory(radar):
    """ Return a function which returns the gate altitudes. """
    try:
        alt = radar.altitude.values + radar.gate_z.values
    except ValueError:
        alt = np.mean(radar.altitude.values) + radar.gate_z.values
    gate_altitude_dict = get_metadata('gate_altitude')
    radar['gate_altitude'] = xr.DataArray(alt, dims=('time', 'range'),
                                          attrs=gate_altitude_dict)
