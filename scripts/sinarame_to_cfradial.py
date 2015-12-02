import pyart
from netcdftime import utime
from datetime import  datetime
import glob, os, sys
import numpy as np

path=sys.argv[1]
path=path+'/'
files=glob.glob(path+'*.H5')

#fieldnames=['CM','KDP','PHIDP','RHOHV','TH','TV','VRAD','WRAD']

for i in np.arange(len(files)):
    basename=os.path.basename(files[i])
    bs=basename.split('_')
    base1='{b1}_{b2}_{b3}_{fn}_{b4}'.format(b1=bs[0],b2=bs[1],b3=bs[2],fn=bs[3],b4=bs[4])
    file='{path}{base1}'.format(path=path,base1=base1)
    radar=pyart.aux_io.read_sinarame_h5(file,file_field_names=True)
    radar.fields[bs[3]]=radar.fields.pop(radar.fields.keys()[0])

    cal_temps = u"gregorian" 
    cdftime = utime(radar.time['units'])

    time1=cdftime.num2date(radar.time['data'][0]).strftime('%Y%m%d_%H%M%S')
    time2=cdftime.num2date(radar.time['data'][-1]).strftime('%Y%m%d_%H%M%S')

    cffile='cfrad.{time1}.0000_to_{time2}.0000_{b1}_{fld}_PPI.nc'.format(time1=time1,time2=time2,b1=bs[0],\
                                                                         fld=bs[3])
    radar._DeflateLevel=5
    print('writing to {path}{cffile}'.format(path=path,cffile=cffile))
    pyart.io.write_cfradial(path+cffile,radar,format='NETCDF4_CLASSIC')

