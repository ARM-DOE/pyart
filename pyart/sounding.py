"""
Class for handling sounding data

"""

import datetime

from scipy.interpolate import interp1d
from numpy import array, where
from pylab import date2num, num2date, datestr2num


__all__ = ['sounding']


class sounding:
    def __init__(self, ob):
        # populate attributes with sounding data, initially this will
        # only work with a netcdf variable object (from Sceintific.IO)
        # but more objects can be added by simply adding elif..
        # PLEASE always populate height in the values['alt'] position and
        # append values['date_list'] and datetime
        # datetime and date_list[index] are datetime objects
        # check if it is a netcdf variable list
        if 'getValue' in dir(ob[ob.keys()[0]]):
            #this is a netcdf variables object
            self.datetime = num2date(datestr2num('19700101') +
                                     ob['base_time'].getValue() /
                                     (24.0*60.0*60.0))
            values = {}
            units = {}
            longname = {}
            for var in ob.keys():
                values.update({var: ob[var][:]})
                try:
                    units.update({var: ob[var].units})
                except AttributeError:
                    units.update({var: "no units"})
                try:
                    longname.update({var: ob[var].long_name})
                except AttributeError:
                    longname.update({var: "no longname"})
                values.update({'date_list': num2date(date2num(self.datetime) +
                               values['time_offset'] / (24.0 * 60.0 * 60.0))})
                units.update({'date_list': 'unitless (object)'})
                self.values = values
                self.units = units
                self.long_name = longname

    def __getattr__(self, attr):
        # This is to allow values in self.values['key'] to be
        # returned as self.key
        if attr in self.values.keys():
            return self.values[attr]
        elif attr in ['values', 'units', 'ydaysdate', 'interpolate_values',
                      'interp_time_base', 'long_name']:
            #need to add to this as new attributes are made
            return self.__getattribute__(attr)

    def ydaysdate(self):
        # return the date as an integer in the form YYDDDhhmm this is
        # for compatibility with some VAD programs and dealias code
        try:  # assume it is a single datetime
            dayofyear = (datetime.datetime(
                self.datetime.year, self.datetime.month, self.datetime.day) -
                datetime.datetime(self.datetime.year, 01, 01)).days + 1
            juldate = ((self.datetime.year - int(("%(d)d" % {'d':
                        self.datetime.year})[0:2]) * 100) * 1000 + dayofyear)
            fulljuldate = (juldate * 10000 + self.datetime.hour * 100 +
                           self.datetime.minute)
        except AttributeError:  # it is a list
            dayofyear = (datetime.datetime(
                self.datetime[0].year, self.datetime[0].month,
                self.datetime[0].day) -
                datetime.datetime(self.datetime[0].year, 01, 01)).days + 1
            juldate = ((self.datetime[0].year - int(("%(d)d" % {'d':
                        self.datetime[0].year})[0:2]) * 100) * 1000 +
                       dayofyear)
            fulljuldate = (juldate * 10000 + self.datetime[0].hour * 100 +
                           self.datetime[0].minute)
        return fulljuldate

    def interpolate_values(self, vals):
        # Destructively (ie inplace) interpolate the sounding onto the
        # values vals
        for key in self.values.keys():
            if key != 'alt' and key != 'date_list':
                if array(self.values[key]).shape == array(self.values['alt']).shape:
                    f = self.interp1d(self.values['alt'], self.values[key],
                                      kind='linear')
                    self.values[key] = f(vals)
        f = self.interp1d(self.values['alt'],
                          date2num(self.values['date_list']), kind='linear')
        self.values['date_list'] = num2date(f(vals))
        self.values['alt'] = vals

    def as_osra_dict(self):
        # for backwards compatibility with OSRA (open source radar functions)
        ret_dict = {}
        lookup = {'alt(m)': 'alt', 'press(hPa)': 'pres', 'wspd(m/s)': 'wspd',
                  'tdry(degs)': 'tdry',  'wdir(degs)': 'deg',
                  'date_list': 'date_list'}
        for key in lookup.keys():
            ret_dict.update({key: self.values[lookup[key]]})
        return ret_dict

    def interp_time_base(self, other_sounding, time_desired):
        # interpolates the current sounding and the other sounding according
        # to the time at launch
        # fist check that time_desired fits in both
        time_between_sondes = (date2num(other_sounding.datetime) -
                               date2num(self.datetime))
        second_wt = ((date2num(time_desired) - date2num(self.datetime)) /
                     time_between_sondes)
        first_wt = (date2num(other_sounding.datetime) -
                    date2num(time_desired)) / time_between_sondes
        if date2num(self.datetime) > date2num(time_desired) > date2num(other_sounding.datetime):
            order = 'before'  # the desired time is before self
        elif date2num(self.datetime) < date2num(time_desired) < date2num(other_sounding.datetime):
            order = 'after'  # the desired time is after self
        else:
            print "time desired is outside of range"
            return

        min_ht = array([self.alt.min(), other_sounding.alt.min()]).max()
        max_ht = array([self.alt.max(), other_sounding.alt.max()]).min()
        idx_min = where(abs(other_sounding.alt-min_ht) ==
                        abs(other_sounding.alt-min_ht).min())[0][0]
        idx_max = where(abs(other_sounding.alt-max_ht) ==
                        abs(other_sounding.alt-max_ht).min())[0][0]
        myhts = other_sounding.alt[idx_min:idx_max]
        self.interpolate_values(myhts)
        altshape = self.alt.shape

        for key in self.values.keys():
            if array(self.values[key]).shape == altshape and key != 'date_list':
                self.values[key] = (
                    self.values[key] * first_wt +
                    other_sounding.values[key][idx_min:idx_max] * second_wt)
        self.values['date_list'] = num2date(date2num(self.values[
            'date_list']) * first_wt + date2num(other_sounding.values[
                'date_list'][idx_min:idx_max]) * second_wt)
        self.datetime = (num2date(date2num(self.datetime) * first_wt +
                         date2num(other_sounding.datetime) * second_wt))
