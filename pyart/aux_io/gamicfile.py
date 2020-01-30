"""
GAMICFile class and utility functions.

"""

import numpy as np
import h5py


class GAMICFile(object):
    """
    A class to read GAMIC files.

    Parameters
    ----------
    filename : str
        Filename of GAMIC HDF5 file.

    Attributes
    ----------
    nsweeps : int
        Number of sweeps (or scans) in the file.
    rays_per_sweep : array of int32
        Number of rays in each sweep.
    total_rays : int
        Total number of rays in all sweeps.
    start_ray, end_ray : array of int32
        Index of the first (start) and last (end) ray in each sweep, 0-based.
    _hfile : HDF5 file
        Open HDF5 file object from which data is read.
    _scans : list
        Name of the HDF5 group for each scan.

    """

    def __init__(self, filename):
        """ initialize object. """
        self._hfile = h5py.File(filename, 'r')
        self.nsweeps = self._hfile['what'].attrs['sets']
        self._scans = ['scan%i' % (i) for i in range(self.nsweeps)]
        self.rays_per_sweep = self.how_attrs('ray_count', 'int32')
        self.total_rays = sum(self.rays_per_sweep)
        # starting and ending ray for each sweep
        self.start_ray = np.cumsum(np.append([0], self.rays_per_sweep[:-1]))
        self.end_ray = np.cumsum(self.rays_per_sweep) - 1
        return

    def close(self):
        """ Close the file. """
        self._hfile.close()

    # file checks
    def is_file_complete(self):
        """ True if all scans in file, False otherwise. """
        # check than scan0, scan1, ... scan[nsweeps-1] are present
        for scan in self._scans:
            if scan not in self._hfile:
                return False
        return True

    def is_file_single_scan_type(self):
        """ True is all scans are the same scan type, False otherwise. """
        scan_type = self._hfile['scan0/what'].attrs['scan_type']
        for scan in self._scans:
            if self._hfile[scan]['what'].attrs['scan_type'] != scan_type:
                return False
        return True

    # attribute look up
    def where_attr(self, attr, dtype):
        """ Return an array containing a attribute from the where group. """
        return np.array([self._hfile['where'].attrs[attr]], dtype=dtype)

    def how_attr(self, attr, dtype):
        """ Return an array containing a attribute from the how group. """
        return np.array([self._hfile['how'].attrs[attr]], dtype=dtype)

    def is_attr_in_group(self, group, attr):
        """ True is attribute is present in the group, False otherwise. """
        return attr in self._hfile[group].attrs

    def raw_group_attr(self, group, attr):
        """ Return an attribute from a group with no reformatting. """
        return self._hfile[group].attrs[attr]

    def raw_scan0_group_attr(self, group, attr):
        """ Return an attribute from the scan0 group with no reformatting. """
        return self._hfile['/scan0'][group].attrs[attr]

    # scan/sweep based attribute lookup
    def how_attrs(self, attr, dtype):
        """ Return an array of an attribute for each scan's how group. """
        return np.array([self._hfile[s]['how'].attrs[attr]
                         for s in self._scans], dtype=dtype)

    def how_ext_attrs(self, attr):
        """
        Return a list of an attribute in each scan's how/extended group.
        """
        return [float(self._hfile[s]['how']['extended'].attrs[attr])
                for s in self._scans]

    def what_attrs(self, attr, dtype):
        """ Return a list of an attribute for each scan's what group. """
        return np.array([self._hfile[s]['what'].attrs[attr]
                         for s in self._scans], dtype=dtype)

    # misc looping
    def moment_groups(self):
        """ Return a list of groups under scan0 where moments are stored. """
        return [k for k in self._hfile['/scan0'] if k.startswith('moment_')]

    def moment_names(self, scan0_groups):
        """ Return a list of moment names for a list of scan0 groups. """
        return [self._hfile['/scan0'][k].attrs['moment'].decode('utf-8') for
                k in scan0_groups]

    def is_field_in_ray_header(self, field):
        """ True if field is present in ray_header, False otherwise. """
        return field in self._hfile[self._scans[0]]['ray_header'].dtype.names

    def ray_header(self, field, dtype):
        """ Return an array containing a ray_header field for each sweep. """
        data = np.empty((self.total_rays, ), dtype=dtype)
        for scan, start, end in zip(self._scans, self.start_ray, self.end_ray):
            data[start:end+1] = self._hfile[scan]['ray_header'][field]
        return data

    def moment_data(self, group, dtype):
        """ Read in moment data from all sweeps. """
        ngates = int(self._hfile['/scan0/how'].attrs['bin_count'])
        data = np.ma.zeros((self.total_rays, ngates), dtype=dtype)
        data[:] = np.ma.masked      # volume data initially all masked
        for scan, start, end in zip(self._scans, self.start_ray, self.end_ray):
            # read in sweep data if field exists in scan.
            if group in self._hfile[scan]:
                sweep_data = _get_gamic_sweep_data(self._hfile[scan][group])
                data[start:end+1, :sweep_data.shape[1]] = sweep_data[:]
        return data

    def sweep_expand(self, arr, dtype='float32'):
        """ Expand an sweep indexed array to be ray indexed """
        return np.repeat(arr, self.rays_per_sweep).astype(dtype)


def _get_gamic_sweep_data(group):
    """ Get GAMIC HDF5 sweep data from an HDF5 group. """
    dyn_range_min = group.attrs['dyn_range_min']
    dyn_range_max = group.attrs['dyn_range_max']
    raw_data = group[:]
    fmt = group.attrs['format']
    if fmt == b'UV16':
        # unsigned 16-bit integer data, 0 indicates a masked value
        assert raw_data.dtype == np.uint16
        scale = (dyn_range_max - dyn_range_min) / 65535.
        offset = dyn_range_min
        sweep_data = np.ma.masked_array(
            raw_data * scale + offset, mask=(raw_data == 0), dtype='float32')
    elif fmt == b'UV8':
        # unsigned 8-bit integer data, 0 indicates a masked value
        assert raw_data.dtype == np.uint8
        scale = (dyn_range_max - dyn_range_min) / 255.
        offset = dyn_range_min
        sweep_data = np.ma.masked_array(
            raw_data * scale + offset, mask=(raw_data == 0), dtype='float32')
    else:
        raise NotImplementedError('GAMIC data format: %s', fmt)
    return sweep_data
