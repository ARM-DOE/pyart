#example code of how to use the grid plotter plus how to add reanalysis...

import pyart
import matplotlib
from netCDF4 import num2date, date2num
import numpy as np
from netCDF4 import Dataset

def reanal_oplot(mapobj, grid, **kwargs):
    clevs = kwargs.get('clevs', np.arange(900,1100.,1.))
    grid_date=num2date(grid.axes['time']['data'], grid.axes['time']['units'])
    datestr=grid_date[0].strftime('%Y%m%d')
    
    address='http://nomads.ncdc.noaa.gov/dods/NCEP_NARR_DAILY/201105/20110520/narr-a_221_'+datestr+'_0000_000'
    data=Dataset(address)
    lats = data.variables['lat'][:]
    lons = data.variables['lon'][:] #(360.0-data.variables['lon'][:])*-1.0

    want_ts=abs(data.variables['time'][:] - \
            date2num(grid_date[0], data.variables['time'].units)).argmin()
    
    prmsl = 0.01*data.variables['prmsl'][want_ts]
    # find x,y of map projection grid.
    lons, lats = np.meshgrid(lons, lats)
    x, y = mapobj(lons, lats)
    cs = mapobj.contour(x,y,prmsl,clevs,colors='k',linewidths=1.)

fname='/data/20110520100000_nex_3d.nc'
grid=pyart.io.grid.read_grid(fname)
font = {'size'   : 16}
matplotlib.rc('font', **font)
f=plt.figure(figsize=[15,8])
pyart.graph.mapplotgrid_3p(grid, f, vmin=-8, vmax=64, level=5, debug=True, lat_cut=36.5, lon_cut=-98.5, overplot=reanal_oplot)
plt.savefig('/home/user/test.png', dpi=200)