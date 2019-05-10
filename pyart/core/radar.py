"""
pyart.core.radar
================

A general central radial scanning (or dwelling) instrument class.

.. autosummary::
    :toctree: generated/

    _rays_per_sweep_data_factory
    _gate_data_factory
    _gate_lon_lat_data_factory
    _gate_altitude_data_factory

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    Radar


"""
from __future__ import print_function

import copy
import sys

import numpy as np

from ..config import get_metadata
from ..lazydict import LazyLoadDict
from .transforms import antenna_vectors_to_cartesian, cartesian_to_geographic


class Radar(object):
    """
    A class for storing antenna coordinate radar data.

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
        Moment fields.
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

                 azimuth, elevation,

                 altitude_agl=None,
                 target_scan_rate=None, rays_are_indexed=None,
                 ray_angle_res=None,

                 scan_rate=None, antenna_transition=None,

                 instrument_parameters=None,
                 radar_calibration=None,

                 rotation=None, tilt=None, roll=None, drift=None, heading=None,
                 pitch=None, georefs_applied=None,

                 ):

        if 'calendar' not in time:
            time['calendar'] = 'gregorian'
        self.time = time
        self.range = _range

        self.fields = fields
        self.metadata = metadata
        self.scan_type = scan_type

        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude
        self.altitude_agl = altitude_agl  # optional

        self.sweep_number = sweep_number
        self.sweep_mode = sweep_mode
        self.fixed_angle = fixed_angle
        self.sweep_start_ray_index = sweep_start_ray_index
        self.sweep_end_ray_index = sweep_end_ray_index

        self.target_scan_rate = target_scan_rate  # optional
        self.rays_are_indexed = rays_are_indexed    # optional
        self.ray_angle_res = ray_angle_res  # optional

        self.azimuth = azimuth
        self.elevation = elevation
        self.scan_rate = scan_rate  # optional
        self.antenna_transition = antenna_transition  # optional
        self.rotation = rotation  # optional
        self.tilt = tilt  # optional
        self.roll = roll  # optional
        self.drift = drift  # optional
        self.heading = heading  # optional
        self.pitch = pitch  # optional
        self.georefs_applied = georefs_applied  # optional

        self.instrument_parameters = instrument_parameters  # optional
        self.radar_calibration = radar_calibration  # optional

        self.ngates = len(_range['data'])
        self.nrays = len(time['data'])
        self.nsweeps = len(sweep_number['data'])
        self.projection = {'proj': 'pyart_aeqd', '_include_lon_0_lat_0': True}

        # initalize attributes with lazy load dictionaries
        self.init_rays_per_sweep()
        self.init_gate_x_y_z()
        self.init_gate_longitude_latitude()
        self.init_gate_altitude()

    def __getstate__(self):
        """ Return object's state which can be pickled. """
        state = self.__dict__.copy()  # copy the objects state
        # Remove unpicklable entries (those which are lazily loaded
        del state['rays_per_sweep']
        del state['gate_x']
        del state['gate_y']
        del state['gate_z']
        del state['gate_longitude']
        del state['gate_latitude']
        del state['gate_altitude']
        return state

    def __setstate__(self, state):
        """ Restore unpicklable entries from pickled object. """
        self.__dict__.update(state)
        self.init_rays_per_sweep()
        self.init_gate_x_y_z()
        self.init_gate_longitude_latitude()
        self.init_gate_altitude()

    # Attribute init/reset method
    def init_rays_per_sweep(self):
        """ Initialize or reset the rays_per_sweep attribute. """
        lazydic = LazyLoadDict(get_metadata('rays_per_sweep'))
        lazydic.set_lazy('data', _rays_per_sweep_data_factory(self))
        self.rays_per_sweep = lazydic

    def init_gate_x_y_z(self):
        """ Initialize or reset the gate_{x, y, z} attributes. """
        gate_x = LazyLoadDict(get_metadata('gate_x'))
        gate_x.set_lazy('data', _gate_data_factory(self, 0))
        self.gate_x = gate_x

        gate_y = LazyLoadDict(get_metadata('gate_y'))
        gate_y.set_lazy('data', _gate_data_factory(self, 1))
        self.gate_y = gate_y

        gate_z = LazyLoadDict(get_metadata('gate_z'))
        gate_z.set_lazy('data', _gate_data_factory(self, 2))
        self.gate_z = gate_z

    def init_gate_longitude_latitude(self):
        """
        Initialize or reset the gate_longitude and gate_latitude attributes.
        """
        gate_longitude = LazyLoadDict(get_metadata('gate_longitude'))
        gate_longitude.set_lazy('data', _gate_lon_lat_data_factory(self, 0))
        self.gate_longitude = gate_longitude

        gate_latitude = LazyLoadDict(get_metadata('gate_latitude'))
        gate_latitude.set_lazy('data', _gate_lon_lat_data_factory(self, 1))
        self.gate_latitude = gate_latitude

    def init_gate_altitude(self):
        """ Initialize the gate_altitude attribute. """
        gate_altitude = LazyLoadDict(get_metadata('gate_altitude'))
        gate_altitude.set_lazy('data', _gate_altitude_data_factory(self))
        self.gate_altitude = gate_altitude

    # private functions for checking limits, etc.
    def _check_sweep_in_range(self, sweep):
        """ Check that a sweep number is in range. """
        if sweep < 0 or sweep >= self.nsweeps:
            raise IndexError('Sweep out of range: ', sweep)
        return

    # public check functions
    def check_field_exists(self, field_name):
        """
        Check that a field exists in the fields dictionary.

        If the field does not exist raise a KeyError.

        Parameters
        ----------
        field_name : str
            Name of field to check.

        """
        if field_name not in self.fields:
            raise KeyError('Field not available: ' + field_name)
        return

    # Iterators

    def iter_start(self):
        """ Return an iterator over the sweep start indices. """
        return (s for s in self.sweep_start_ray_index['data'])

    def iter_end(self):
        """ Return an iterator over the sweep end indices. """
        return (s for s in self.sweep_end_ray_index['data'])

    def iter_start_end(self):
        """ Return an iterator over the sweep start and end indices. """
        return ((s, e) for s, e in zip(self.iter_start(), self.iter_end()))

    def iter_slice(self):
        """ Return an iterator which returns sweep slice objects. """
        return (slice(s, e+1) for s, e in self.iter_start_end())

    def iter_field(self, field_name):
        """ Return an iterator which returns sweep field data. """
        self.check_field_exists(field_name)
        return (self.fields[field_name]['data'][s] for s in self.iter_slice())

    def iter_azimuth(self):
        """ Return an iterator which returns sweep azimuth data. """
        return (self.azimuth['data'][s] for s in self.iter_slice())

    def iter_elevation(self):
        """ Return an iterator which returns sweep elevation data. """
        return (self.elevation['data'][s] for s in self.iter_slice())

    # get methods

    def get_start(self, sweep):
        """ Return the starting ray index for a given sweep. """
        self._check_sweep_in_range(sweep)
        return self.sweep_start_ray_index['data'][sweep]

    def get_end(self, sweep):
        """ Return the ending ray for a given sweep. """
        self._check_sweep_in_range(sweep)
        return self.sweep_end_ray_index['data'][sweep]

    def get_start_end(self, sweep):
        """ Return the starting and ending ray for a given sweep. """
        return self.get_start(sweep), self.get_end(sweep)

    def get_slice(self, sweep):
        """ Return a slice for selecting rays for a given sweep. """
        start, end = self.get_start_end(sweep)
        return slice(start, end+1)

    def get_field(self, sweep, field_name, copy=False):
        """
        Return the field data for a given sweep.

        When used with :py:func:`get_gate_x_y_z` this method can be used to
        obtain the data needed for plotting a radar field with the correct
        spatial context.

        Parameters
        ----------
        sweep : int
            Sweep number to retrieve data for, 0 based.
        field_name : str
            Name of the field from which data should be retrieved.
        copy : bool, optional
            True to return a copy of the data. False, the default, returns
            a view of the data (when possible), changing this data will
            change the data in the underlying Radar object.

        Returns
        -------
        data : array
            Array containing data for the requested sweep and field.

        """
        self.check_field_exists(field_name)
        s = self.get_slice(sweep)
        data = self.fields[field_name]['data'][s]
        if copy:
            return data.copy()
        else:
            return data

    def get_azimuth(self, sweep, copy=False):
        """
        Return an array of azimuth angles for a given sweep.

        Parameters
        ----------
        sweep : int
            Sweep number to retrieve data for, 0 based.
        copy : bool, optional
            True to return a copy of the azimuths. False, the default, returns
            a view of the azimuths (when possible), changing this data will
            change the data in the underlying Radar object.

        Returns
        -------
        azimuths : array
            Array containing the azimuth angles for a given sweep.

        """
        s = self.get_slice(sweep)
        azimuths = self.azimuth['data'][s]
        if copy:
            return azimuths.copy()
        else:
            return azimuths

    def get_elevation(self, sweep, copy=False):
        """
        Return an array of elevation angles for a given sweep.

        Parameters
        ----------
        sweep : int
            Sweep number to retrieve data for, 0 based.
        copy : bool, optional
            True to return a copy of the elevations. False, the default,
            returns a view of the elevations (when possible), changing this
            data will change the data in the underlying Radar object.

        Returns
        -------
        azimuths : array
            Array containing the elevation angles for a given sweep.

        """
        s = self.get_slice(sweep)
        elevation = self.elevation['data'][s]
        if copy:
            return elevation.copy()
        else:
            return elevation

    def get_gate_x_y_z(self, sweep, edges=False, filter_transitions=False):
        """
        Return the x, y and z gate locations in meters for a given sweep.

        With the default parameter this method returns the same data as
        contained in the gate_x, gate_y and gate_z attributes but this method
        performs the gate location calculations only for the specified sweep
        and therefore is more efficient than accessing this data through these
        attribute.

        When used with :py:func:`get_field` this method can be used to obtain
        the data needed for plotting a radar field with the correct spatial
        context.

        Parameters
        ----------
        sweep : int
            Sweep number to retrieve gate locations from, 0 based.
        edges : bool, optional
            True to return the locations of the gate edges calculated by
            interpolating between the range, azimuths and elevations.
            False (the default) will return the locations of the gate centers
            with no interpolation.
        filter_transitions : bool, optional
            True to remove rays where the antenna was in transition between
            sweeps. False will include these rays. No rays will be removed
            if the antenna_transition attribute is not available (set to None).

        Returns
        -------
        x, y, z : 2D array
            Array containing the x, y and z, distances from the radar in
            meters for the center (or edges) for all gates in the sweep.

        """
        azimuths = self.get_azimuth(sweep)
        elevations = self.get_elevation(sweep)

        if filter_transitions and self.antenna_transition is not None:
            sweep_slice = self.get_slice(sweep)
            valid = self.antenna_transition['data'][sweep_slice] == 0
            azimuths = azimuths[valid]
            elevations = elevations[valid]

        return antenna_vectors_to_cartesian(
            self.range['data'], azimuths, elevations, edges=edges)

    def get_gate_lat_lon_alt(self, sweep, reset_gate_coords=False,
                             filter_transitions=False):
        """
        Return the longitude, latitude and altitude gate locations.
        Longitude and latitude are in degrees and altitude in meters.

        With the default parameter this method returns the same data as
        contained in the gate_latitude, gate_longitude and gate_altitude
        attributes but this method performs the gate location calculations
        only for the specified sweep and therefore is more efficient than
        accessing this data through these attribute. If coordinates have
        at all, please use the reset_gate_coords parameter.

        Parameters
        ----------
        sweep : int
            Sweep number to retrieve gate locations from, 0 based.
        reset_gate_coords : bool, optional
            Optional to reset the gate latitude, gate longitude and gate
            altitude attributes before using them in this function. This
            is useful when the geographic coordinates have changed and gate
            latitude, gate longitude and gate altitude need to be reset.
        filter_transitions : bool, optional
            True to remove rays where the antenna was in transition between
            sweeps. False will include these rays. No rays will be removed
            if the antenna_transition attribute is not available (set to None).

        Returns
        -------
        lat, lon, alt : 2D array
            Array containing the latitude, longitude and altitude,
            for all gates in the sweep.

        """
        s = self.get_slice(sweep)

        if reset_gate_coords:
            gate_latitude = LazyLoadDict(get_metadata('gate_latitude'))
            gate_latitude.set_lazy('data', _gate_lon_lat_data_factory(self, 1))
            self.gate_latitude = gate_latitude

            gate_longitude = LazyLoadDict(get_metadata('gate_longitude'))
            gate_longitude.set_lazy('data', _gate_lon_lat_data_factory(self, 0))
            self.gate_longitude = gate_longitude

            gate_altitude = LazyLoadDict(get_metadata('gate_altitude'))
            gate_altitude.set_lazy('data', _gate_altitude_data_factory(self))
            self.gate_altitude = gate_altitude

        lat = self.gate_latitude['data'][s]
        lon = self.gate_longitude['data'][s]
        alt = self.gate_altitude['data'][s]

        if filter_transitions and self.antenna_transition is not None:
            valid = self.antenna_transition['data'][s] == 0
            lat = lat[valid]
            lon = lon[valid]
            alt = alt[valid]

        return lat, lon, alt

    def get_nyquist_vel(self, sweep, check_uniform=True):
        """
        Return the Nyquist velocity in meters per second for a given sweep.

        Raises a LookupError if the Nyquist velocity is not available, an
        Exception is raised if the velocities are not uniform in the sweep
        unless check_uniform is set to False.

        Parameters
        ----------
        sweep : int
            Sweep number to retrieve data for, 0 based.
        check_uniform : bool
            True to check to perform a check on the Nyquist velocities that
            they are uniform in the sweep, False will skip this check and
            return the velocity of the first ray in the sweep.

        Returns
        -------
        nyquist_velocity : float
            Array containing the Nyquist velocity in m/s for a given sweep.

        """
        s = self.get_slice(sweep)
        try:
            nyq_vel = self.instrument_parameters['nyquist_velocity']['data'][s]
        except:
            raise LookupError('Nyquist velocity unavailable')
        if check_uniform:
            if np.any(nyq_vel != nyq_vel[0]):
                raise Exception('Nyquist velocities are not uniform in sweep')
        return float(nyq_vel[0])

    # Methods

    def info(self, level='standard', out=sys.stdout):
        """
        Print information on radar.

        Parameters
        ----------
        level : {'compact', 'standard', 'full', 'c', 's', 'f'}, optional
            Level of information on radar object to print, compact is
            minimal information, standard more and full everything.
        out : file-like, optional
            Stream to direct output to, default is to print information
            to standard out (the screen).

        """
        if level == 'c':
            level = 'compact'
        elif level == 's':
            level = 'standard'
        elif level == 'f':
            level = 'full'

        if level not in ['standard', 'compact', 'full']:
            raise ValueError('invalid level parameter')

        self._dic_info('altitude', level, out)
        self._dic_info('altitude_agl', level, out)
        self._dic_info('antenna_transition', level, out)
        self._dic_info('azimuth', level, out)
        self._dic_info('elevation', level, out)

        print('fields:', file=out)
        for field_name, field_dic in self.fields.items():
            self._dic_info(field_name, level, out, field_dic, 1)

        self._dic_info('fixed_angle', level, out)

        if self.instrument_parameters is None:
            print('instrument_parameters: None', file=out)
        else:
            print('instrument_parameters:', file=out)
            for name, dic in self.instrument_parameters.items():
                self._dic_info(name, level, out, dic, 1)

        self._dic_info('latitude', level, out)
        self._dic_info('longitude', level, out)

        print('nsweeps:', self.nsweeps, file=out)
        print('ngates:', self.ngates, file=out)
        print('nrays:', self.nrays, file=out)

        if self.radar_calibration is None:
            print('radar_calibration: None', file=out)
        else:
            print('radar_calibration:', file=out)
            for name, dic in self.radar_calibration.items():
                self._dic_info(name, level, out, dic, 1)

        self._dic_info('range', level, out)
        self._dic_info('scan_rate', level, out)
        print('scan_type:', self.scan_type, file=out)
        self._dic_info('sweep_end_ray_index', level, out)
        self._dic_info('sweep_mode', level, out)
        self._dic_info('sweep_number', level, out)
        self._dic_info('sweep_start_ray_index', level, out)
        self._dic_info('target_scan_rate', level, out)
        self._dic_info('time', level, out)

        # Airborne radar parameters
        if self.rotation is not None:
            self._dic_info('rotation', level, out)
        if self.tilt is not None:
            self._dic_info('tilt', level, out)
        if self.roll is not None:
            self._dic_info('roll', level, out)
        if self.drift is not None:
            self._dic_info('drift', level, out)
        if self.heading is not None:
            self._dic_info('heading', level, out)
        if self.pitch is not None:
            self._dic_info('pitch', level, out)
        if self.georefs_applied is not None:
            self._dic_info('georefs_applied', level, out)

        # always print out all metadata last
        self._dic_info('metadata', 'full', out)

    def _dic_info(self, attr, level, out, dic=None, ident_level=0):
        """ Print information on a dictionary attribute. """
        if dic is None:
            dic = getattr(self, attr)

        ilvl0 = '\t' * ident_level
        ilvl1 = '\t' * (ident_level + 1)

        if dic is None:
            print(str(attr) + ': None', file=out)
            return

        # make a string summary of the data key if it exists.
        if 'data' not in dic:
            d_str = 'Missing'
        elif not isinstance(dic['data'], np.ndarray):
            d_str = '<not a ndarray>'
        else:
            data = dic['data']
            t = (data.dtype, data.shape)
            d_str = '<ndarray of type: %s and shape: %s>' % t

        # compact, only data summary
        if level == 'compact':
            print(ilvl0 + str(attr) + ':', d_str, file=out)

        # standard, all keys, only summary for data
        elif level == 'standard':
            print(ilvl0 + str(attr) + ':', file=out)
            print(ilvl1 + 'data:', d_str, file=out)
            for key, val in dic.items():
                if key == 'data':
                    continue
                print(ilvl1 + key + ':', val, file=out)

        # full, all keys, full data
        elif level == 'full':
            print(str(attr) + ':', file=out)
            if 'data' in dic:
                print(ilvl1 + 'data:', dic['data'], file=out)
            for key, val in dic.items():
                if key == 'data':
                    continue
                print(ilvl1 + key + ':', val, file=out)

        return

    def add_field(self, field_name, dic, replace_existing=False):
        """
        Add a field to the object.

        Parameters
        ----------
        field_name : str
            Name of the field to add to the dictionary of fields.
        dic : dict
            Dictionary contain field data and metadata.
        replace_existing : bool, optional
            True to replace the existing field with key field_name if it
            exists, loosing any existing data. False will raise a ValueError
            when the field already exists.

        """
        # check that the field dictionary to add is valid
        if field_name in self.fields and replace_existing is False:
            err = 'A field with name: %s already exists' % (field_name)
            raise ValueError(err)
        if 'data' not in dic:
            raise KeyError("dic must contain a 'data' key")
        if dic['data'].shape != (self.nrays, self.ngates):
            t = (self.nrays, self.ngates)
            err = "'data' has invalid shape, should be (%i, %i)" % t
            raise ValueError(err)
        # add the field
        self.fields[field_name] = dic
        return

    def add_field_like(self, existing_field_name, field_name, data,
                       replace_existing=False):
        """
        Add a field to the object with metadata from a existing field.

        Note that the data parameter is not copied by this method.
        If data refers to a 'data' array from an existing field dictionary, a
        copy should be made within or prior to using this method. If this is
        not done the 'data' key in both field dictionaries will point to the
        same NumPy array and modification of one will change the second. To
        copy NumPy arrays use the copy() method. See the Examples section
        for how to create a copy of the 'reflectivity' field as a field named
        'reflectivity_copy'.

        Parameters
        ----------
        existing_field_name : str
            Name of an existing field to take metadata from when adding
            the new field to the object.
        field_name : str
            Name of the field to add to the dictionary of fields.
        data : array
            Field data. A copy of this data is not made, see the note above.
        replace_existing : bool, optional
            True to replace the existing field with key field_name if it
            exists, loosing any existing data. False will raise a ValueError
            when the field already exists.

        Examples
        --------
        >>> radar.add_field_like('reflectivity', 'reflectivity_copy',
        ...                      radar.fields['reflectivity']['data'].copy())

        """
        if existing_field_name not in self.fields:
            err = 'field %s does not exist in object' % (existing_field_name)
            raise ValueError(err)
        dic = {}
        for k, v in self.fields[existing_field_name].items():
            if k != 'data':
                dic[k] = v
        dic['data'] = data
        return self.add_field(field_name, dic,
                              replace_existing=replace_existing)

    def extract_sweeps(self, sweeps):
        """
        Create a new radar contains only the data from select sweeps.

        Parameters
        ----------
        sweeps : array_like
            Sweeps (0-based) to include in new Radar object.

        Returns
        -------
        radar : Radar
            Radar object which contains a copy of data from the selected
            sweeps.

        """

        # parse and verify parameters
        sweeps = np.array(sweeps, dtype='int32')
        if np.any(sweeps > (self.nsweeps - 1)):
            raise ValueError('invalid sweeps indices in sweeps parameter')
        if np.any(sweeps < 0):
            raise ValueError('only positive sweeps can be extracted')

        def mkdic(dic, select):
            """ Make a dictionary, selecting out select from data key """
            if dic is None:
                return None
            d = dic.copy()
            if 'data' in d and select is not None:
                d['data'] = d['data'][select].copy()
            return d

        # create array of rays which select the sweeps selected and
        # the number of rays per sweep.
        ray_count = (self.sweep_end_ray_index['data'] -
                     self.sweep_start_ray_index['data'] + 1)[sweeps]
        ssri = self.sweep_start_ray_index['data'][sweeps]
        rays = np.concatenate(
            [range(s, s+e) for s, e in zip(ssri, ray_count)]).astype('int32')

        # radar location attribute dictionary selector
        if len(self.altitude['data']) == 1:
            loc_select = None
        else:
            loc_select = sweeps

        # create new dictionaries
        time = mkdic(self.time, rays)
        _range = mkdic(self.range, None)

        fields = {}
        for field_name, dic in self.fields.items():
            fields[field_name] = mkdic(dic, rays)
        metadata = mkdic(self.metadata, None)
        scan_type = str(self.scan_type)

        latitude = mkdic(self.latitude, loc_select)
        longitude = mkdic(self.longitude, loc_select)
        altitude = mkdic(self.altitude, loc_select)
        altitude_agl = mkdic(self.altitude_agl, loc_select)

        sweep_number = mkdic(self.sweep_number, sweeps)
        sweep_mode = mkdic(self.sweep_mode, sweeps)
        fixed_angle = mkdic(self.fixed_angle, sweeps)
        sweep_start_ray_index = mkdic(self.sweep_start_ray_index, None)
        sweep_start_ray_index['data'] = np.cumsum(
            np.append([0], ray_count[:-1]), dtype='int32')
        sweep_end_ray_index = mkdic(self.sweep_end_ray_index, None)
        sweep_end_ray_index['data'] = np.cumsum(ray_count, dtype='int32') - 1
        target_scan_rate = mkdic(self.target_scan_rate, sweeps)

        azimuth = mkdic(self.azimuth, rays)
        elevation = mkdic(self.elevation, rays)
        scan_rate = mkdic(self.scan_rate, rays)
        antenna_transition = mkdic(self.antenna_transition, rays)

        # instrument_parameters
        # Filter the instrument_parameter dictionary based size of leading
        # dimension, this might not always be correct.
        if self.instrument_parameters is None:
            instrument_parameters = None
        else:
            instrument_parameters = {}
            for key, dic in self.instrument_parameters.items():
                if dic['data'].ndim != 0:
                    dim0_size = dic['data'].shape[0]
                else:
                    dim0_size = -1
                if dim0_size == self.nsweeps:
                    fdic = mkdic(dic, sweeps)
                elif dim0_size == self.nrays:
                    fdic = mkdic(dic, rays)
                else:   # keep everything
                    fdic = mkdic(dic, None)
                instrument_parameters[key] = fdic

        # radar_calibration
        # copy all field in radar_calibration as is except for
        # r_calib_index which we filter based upon time.  This might
        # leave some indices in the "r_calib" dimension not referenced in
        # the r_calib_index array.
        if self.radar_calibration is None:
            radar_calibration = None
        else:
            radar_calibration = {}
            for key, dic in self.radar_calibration.items():
                if key == 'r_calib_index':
                    radar_calibration[key] = mkdic(dic, rays)
                else:
                    radar_calibration[key] = mkdic(dic, None)

        return Radar(time, _range, fields, metadata, scan_type,
                     latitude, longitude, altitude,
                     sweep_number, sweep_mode, fixed_angle,
                     sweep_start_ray_index, sweep_end_ray_index,
                     azimuth, elevation,
                     altitude_agl=altitude_agl,
                     target_scan_rate=target_scan_rate,
                     scan_rate=scan_rate,
                     antenna_transition=antenna_transition,
                     instrument_parameters=instrument_parameters,
                     radar_calibration=radar_calibration)


def _rays_per_sweep_data_factory(radar):
    """ Return a function which returns the number of rays per sweep. """
    def _rays_per_sweep_data():
        """ The function which returns the number of rays per sweep. """
        return (radar.sweep_end_ray_index['data'] -
                radar.sweep_start_ray_index['data'] + 1)
    return _rays_per_sweep_data


def _gate_data_factory(radar, coordinate):
    """ Return a function which returns the Cartesian locations of gates. """
    def _gate_data():
        """ The function which returns the Cartesian locations of gates. """
        ranges = radar.range['data']
        azimuths = radar.azimuth['data']
        elevations = radar.elevation['data']
        cartesian_coords = antenna_vectors_to_cartesian(
            ranges, azimuths, elevations, edges=False)
        # load x, y, and z data except for the coordinate in question
        if coordinate != 0:
            radar.gate_x['data'] = cartesian_coords[0]
        if coordinate != 1:
            radar.gate_y['data'] = cartesian_coords[1]
        if coordinate != 2:
            radar.gate_z['data'] = cartesian_coords[2]
        return cartesian_coords[coordinate]
    return _gate_data


def _gate_lon_lat_data_factory(radar, coordinate):
    """ Return a function which returns the geographic locations of gates. """
    def _gate_lon_lat_data():
        """ The function which returns the geographic locations gates. """
        x = radar.gate_x['data']
        y = radar.gate_y['data']
        projparams = radar.projection.copy()
        if projparams.pop('_include_lon_0_lat_0', False):
            projparams['lon_0'] = radar.longitude['data'][0]
            projparams['lat_0'] = radar.latitude['data'][0]
        geographic_coords = cartesian_to_geographic(x, y, projparams)
        # set the other geographic coordinate
        if coordinate == 0:
            radar.gate_latitude['data'] = geographic_coords[1]
        else:
            radar.gate_longitude['data'] = geographic_coords[0]
        return geographic_coords[coordinate]
    return _gate_lon_lat_data


def _gate_altitude_data_factory(radar):
    """ Return a function which returns the gate altitudes. """
    def _gate_altitude_data():
        """ The function which returns the gate altitudes. """
        try:
            return radar.altitude['data'] + radar.gate_z['data']
        except ValueError:
            return np.mean(radar.altitude['data']) + radar.gate_z['data']
    return _gate_altitude_data
