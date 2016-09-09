__author__ = 'wolfensb'

import numpy as np
import matplotlib.pyplot as plt
from cosmo_pol.radar import pyart_wrapper
from pyart.config import get_field_name
from pyart.core import transforms
from scipy import ndimage,interpolate,spatial

KE = 4/3.

def get_earth_radius(latitude):
    # Two extreme earth radius
    a=6378.1370*1000
    b=6356.7523*1000

    earth_radius=np.sqrt(((a**2*np.cos(latitude))**2+
        (b**2*np.sin(latitude))**2)/((a*np.cos(latitude))**2+
        (b*np.sin(latitude))**2))
        
    return earth_radius

def polar_to_carthesian(radar, sweep, variable_name, cart_res = 75,
                        max_range = 20000, transform = None,
                        inverse_transform = None):  

    is_ppi = radar.sweep_mode['data'][0] == 'ppi'

   # Get angles and distances of radar data
    if is_ppi:
        theta = part.azimuth['data']
    else:
        theta = part.elevation['data']
        
    r = radar.range['data']
    
    # Get data to be interpolated
    pol_data = radar.get_field(sweep,variable_name)
    r = r[r<max_range]
    
    # Cut data at max_range
    pol_data = pol_data[:,r<max_range]
    
    # If needed perform transform (linearize)
    if transform:
        pol_data = transform(pol_data)
        
    # Unmask array
    pol_data = np.ma.masked_values(pol_data, np.nan)

    # One specificity of using the kd-tree is that we need to pad the array
    # with nans at large ranges and angles smaller and larger
    pol_data = np.pad(pol_data,pad_width=((1,1),(0,1)),
                      mode='constant',constant_values=np.nan)
    # Get dim of data
    [N,M] = pol_data.shape

    # We need to pad theta and r as well    
    theta = np.hstack([theta[0]-0.1,theta,theta[-1]+0.1])
    r = np.hstack([r,r[-1]+0.1])
    
    r_grid_p, theta_grid_p = np.meshgrid(r,theta)
    
    # Generate regular cartesian grid
    if is_ppi:
        x_vec = np.arange(-max_range-cart_res,max_range+cart_res,cart_res)
        y_vec = np.arange(-max_range-cart_res,max_range+cart_res,cart_res)
    else:
        x_vec = np.arange((max_range-cart_res)*np.cos(np.radians(np.max(theta)))
                          ,max_range+cart_res,cart_res)
        y_vec = np.arange(0,max_range+cart_res,cart_res)
        
    x_grid_c, y_grid_c = np.meshgrid(x_vec, y_vec)
    
    if is_ppi:
        theta_grid_c= np.degrees(np.arctan2(-x_grid_c, -y_grid_c)+np.pi)
        r_grid_c = (np.sqrt(x_grid_c**2 + y_grid_c**2))
    else:
        theta_grid_c = np.degrees(-(np.arctan2(x_grid_c, y_grid_c)-np.pi/2))
        E = get_earth_radius(radar.latitude['data'])
        r_grid_c = -KE*E*np.sin(np.radians(theta_grid_c)) - \
            np.sqrt(2*((y_grid_c + KE*E)**2 + (KE*E)**2))
        print(r_grid_c)
    
    # Kd-tree triangulation

    kdtree = spatial.cKDTree(np.vstack((r_grid_p.ravel(),theta_grid_p.ravel())).T)
    _,idx_pts = kdtree.query(np.vstack((r_grid_c.ravel(),
                                        theta_grid_c.ravel())).T, k=1)
                                        
    cart_data = pol_data.ravel()[idx_pts]
    cart_data = np.reshape(cart_data,x_grid_c.shape)
    
    # If needed perform inverse transform
    if inverse_transform:
        cart_data = inverse_transform(cart_data)
        
    return (x_grid_c,y_grid_c), cart_data
#    
if __name__ == '__main__':
    file_rad ='/home/daniel/MXPol-polar-20140322-120940-RHI-048_0.nc'

    part = pyart_wrapper.PyradMXPOL(file_rad)
    loc,c = polar_to_carthesian(part,0,'Zh',max_range=20000,transform = lambda x: 10**(0.1*x),
                                inverse_transform = lambda x: 10*np.log10(x))
    plt.contourf(loc[0],loc[1],c)
#    coords, vals = polar_to_carthesian(part,0,variable_name = 'Zh',max_range=20000)
#    plt.figure()

#        if(variableName in ['Ph','Pv','SNRh','SNRv','Zh','Zv','MZh','MZv']):
#            varLin=10**(0.1*var)
#            if(typeInterp=='nearest'):
#                if pol2cart == -1 and cart2pol == -1:                
#                    data, pol2cart, cart2pol =radarBinInterpolation(varPol['range'], varPol['elevations'], x,z,varLin,xVec, zVec,'RHI')
#                else:
#                    data=radarBinInterpolationFast(varLin,pol2cart,cart2pol)


#    
#def projectRHI(varPol, listVariables, coords_radar, typethetaInterp='nearest', interp_res=75, maxRange=20000,maxHeight=-1, pol2cart=-1, cart2pol=-1):
#    # Projects a RHI in polar coordinates (range-elevation) on a Cartesian grid
#    # inputs are the name of the NetCDF input file, a list of variables to project ['Zh','Rhohv'], for example, the resolution
#    # (grid size) of the Cartesian grid, the maximum elevation angle to be projected,
#    #  the maximum range from the radar to be projected and the maximum height above radar to be projected
#
#    latitude=coords_radar[0]
#    
#    
#    maxElevRad=np.max(varPol['elevations'])*DEG2RAD  
#    
#    # Construct meshgrid of range-elevation
#
#    elevGrid, rGrid=np.meshgrid(varPol['elevations']*DEG2RAD,varPol['range'])
#    rVec=rGrid.ravel()
#    elevVec=elevGrid.ravel()
#
#    # Project x and y coordinates
#    x, y, z = sph2cart(0,elevVec,rVec)
#
#    x=y # due to how we defined sph2cart, x will be only 0 and y will be the horizontal coordinates, we swap this
#    # Project z coordinate
#    ke=4/3
#    # Two extreme earth radius
#    a=6378.1370*1000
#    b=6356.7523*1000
#
#    EarthRadius=np.sqrt(((a**2*np.cos(latitude))**2+(b**2*np.sin(latitude))**2)/((a*np.cos(latitude))**2+(b*np.sin(latitude))**2))
#    z=np.sqrt(rVec**2 + (ke*EarthRadius)**2+2*rVec*ke*EarthRadius*np.sin(elevVec))-ke*EarthRadius
#
#    # Generate cartesian grid
#    xMin=-maxRange*np.sin(max(maxElevRad,np.pi/2)-np.pi/2)
#    xMax=maxRange
#    zMin=0
#    if(maxHeight==-1):
#      try:
#          zMax=max(z[~np.isnan(varPol[listVariables[0]].ravel())])
#      except:
#          zMax=10000 # Basically there's nothing on the scan. do that to avoid crash
#    else:
#      zMax=maxHeight
#
#    xVec=np.arange(xMin, xMax, interp_res)
#    zVec=np.arange(zMin, zMax, interp_res)
#
#    xGrid, zGrid = np.meshgrid(xVec, zVec)
#
#    # Create output dic
#
#    output={}
#    output['variables']={}
#    output['coordinates']={}
#    output['attributes']={}
#    
#    # Create Delaunay triangulation
#    for variableName in listVariables:
#        var=varPol[variableName]
#        if(variableName in ['Ph','Pv','SNRh','SNRv','Zh','Zv','MZh','MZv']):
#            varLin=10**(0.1*var)
#            if(typeInterp=='nearest'):
#                if pol2cart == -1 and cart2pol == -1:                
#                    data, pol2cart, cart2pol =radarBinInterpolation(varPol['range'], varPol['elevations'], x,z,varLin,xVec, zVec,'RHI')
#                else:
#                    data=radarBinInterpolationFast(varLin,pol2cart,cart2pol)
#            else:
#                data = griddata((x, z), varLin.ravel(),(xGrid, zGrid), method=typeInterp)
#            dataInterp=10*np.log10(data)
#        elif(variableName=='Sw'):
#            varLin=var**2
#            if(typeInterp=='nearest'):
#                if pol2cart == -1 and cart2pol == -1:                
#                    data, pol2cart, cart2pol =radarBinInterpolation(varPol['range'], varPol['elevations'], x,z,varLin,xVec, zVec,'RHI')
#                else:
#                    data=radarBinInterpolationFast(varLin,pol2cart,cart2pol)
#            else:
#                data,  = griddata((x, z), varLin.ravel(),(xGrid, zGrid), method=typeInterp)
#            dataInterp=np.sqrt(data)
#        else:
#            if(typeInterp=='nearest'): 
#                if pol2cart == -1 and cart2pol == -1:                
#                    data, pol2cart, cart2pol =radarBinInterpolation(varPol['range'], varPol['elevations'], x,z,var,xVec, zVec,'RHI')
#                else:
#                    data=radarBinInterpolationFast(var,pol2cart,cart2pol)
#            else:
#                data = griddata((x, z), var.ravel(),(xGrid, zGrid), method=typeInterp)
#            dataInterp=data
#        #if(typeInterp=='nearest'):
#        #    correctArtifacts(dataInterp, xGrid, yGrid, maxRangeProj, (np.nanmin(azimuth_rad), np.nanmax(azimuth_rad)))
#            
#        output['variables'][variableName]=dataInterp
#
#    output['attributes']['resolution']=interp_res
#    output['coordinates']['Xvec']=xVec
#    output['coordinates']['Zvec']=zVec
#    output['attributes']['maxRange']=maxRange
#    output['attributes']['maxElev']=np.max(varPol['elevations'])
#
#    return output, pol2cart, cart2pol
#
#
#def getCoordsPPI(filename,maxRange):
#    varPol=readRadData(filename, maxRange)
#    if(len(varPol['azimuth'])==1):
#        print 'Dimension of input file not consistent...'
#        print 'Please give a PPI scan as input!'
#        print 'Aborting...'
#        return
#    else:
#        azimuth=varPol['azimuth']
#
#    azimuth_rad=azimuth*DEG2RAD
#
#    elevation=varPol['elevation']
#    elevation_pol=elevation*DEG2RAD
#    range=varPol['range']
#    maxRangeProj=maxRange*np.sin(np.pi/2-(elevation_pol))
#
#    aziGrid, rGrid=np.meshgrid(azimuth_rad,range)
#
#    rVec=rGrid.ravel()
#    aziVec=aziGrid.ravel()
#
#    # Project and get x, y and z coordinates
#    x, y, z = sph2cart(aziVec,elevation_pol,rVec)
#    return x,y,z
#
#
#
#def projectPPI(filename, listVariables, typeInterp, interp_res, maxRange):
#    # Projects a PPI in polar coordinates (range-azimuth) on a Cartesian grid
#    # inputs are the name of the NetCDF input file, a list of variables to project ['Zh','Rhohv'], for example, the resolution
#    # (grid size) of the Cartesian grid, and the maximum range from the radar to be projected
#
#    varPol=readRadData(filename, listVariables)
#    if(len(varPol['azimuth'])==1):
#        print 'Dimension of input file not consistent...'
#        print 'Please give a PPI scan as input!'
#        print 'Aborting...'
#        return
#    else:
#        azimuth=varPol['azimuth']
#    if(azimuth[0]>azimuth[1]):
#      azimuth = azimuth[::-1]
#
#    azimuth_rad=azimuth*DEG2RAD
#
#    elevation=varPol['elevation']
#    elevation_pol=elevation*DEG2RAD
#    range=varPol['range']
#    maxRangeProj=maxRange*np.sin(np.pi/2-(elevation_pol))
#
#    aziGrid, rGrid=np.meshgrid(azimuth_rad,range)
#
#    rVec=rGrid.ravel()
#    aziVec=aziGrid.ravel()
#
#    # Project and get x, y and z coordinates
#    x, y, z = sph2cart(aziVec,elevation_pol,rVec)
#
#    # Now we generate a regular grid and assign to each pixel the correct value
#    # of the projected points
#    xVec= np.arange(start=-maxRangeProj,stop=maxRangeProj,step=interp_res)
#    yVec= xVec
#    [xGrid, yGrid] = np.meshgrid(xVec, yVec)
#
#    R=np.sqrt(xGrid**2+yGrid**2)
#    zGrid=R*np.sin(elevation_pol)
#
#      # Create output dic
#    output={}
#    output['variables']={}
#    output['coordinates']={}
#    output['attributes']={}
#
#    pol2cart=-1
#    cart2pol=-1
#    
#    # Create Delaunay triangulation
#    for variableName in listVariables:
#        var=varPol[variableName]
#        if(variableName in ['Z','Pv','SNRh','SNRv','MZh','MZv', 'Zv', 'Zh']):
#            varLin=10**(0.1*var)
#            if(typeInterp=='nearest'):
#                if pol2cart == -1 and cart2pol == -1:
#                    data, pol2cart, cart2pol =radarBinInterpolation(range, azimuth, x,y,varLin,xVec, yVec,'PPI')
#                else:
#                    data=radarBinInterpolationFast(varLin,pol2cart,cart2pol)
#            else:
#                data = griddata((x, y), varLin,(xGrid, yGrid), method=typeInterp)
#            dataInterp=10*np.log10(data)
#        elif(variableName=='Sw'):
#            varLin=var**2
#            if(typeInterp=='nearest'):
#                if pol2cart == -1 and cart2pol == -1:
#                    data, pol2cart, cart2pol =radarBinInterpolation(range, azimuth, x,y,varLin,xVec, yVec,'PPI')
#                else:
#                    data=radarBinInterpolationFast(varLin,pol2cart,cart2pol)
#            else:
#                 data  = griddata((x, y), varLin,(xGrid, yGrid), method=typeInterp)
#            dataInterp=np.sqrt(data)
#        else:
#
#            if(typeInterp=='nearest'):
#                if pol2cart == -1 and cart2pol == -1:
#                    data, pol2cart, cart2pol =radarBinInterpolation(range, azimuth, x,y,var,xVec, yVec,'PPI')
#                else:
#                    data=radarBinInterpolationFast(var,pol2cart,cart2pol)
#            else:
#                data  = griddata((x, y), var,(xGrid, yGrid), method=typeInterp)
#            dataInterp=data
#
#        output['variables'][variableName]=dataInterp
#
#    output['attributes']['resolution']=interp_res
#    output['coordinates']['Xvec']=xVec
#    output['coordinates']['Yvec']=yVec
#    output['coordinates']['Z']=zGrid
#    output['attributes']['maxRange']=maxRange
#    
#
#    return output
#
#def radarBinInterpolationFast(radvar, pol2cart, cart2pol):
#
#     # Interpolates radar variable radvar sampled at polar coordinates p1, p2 that correspond to cartesian coordinates
#    # x1, x2 to a regular cartesian grid given by vectors X1_vec and X2_vec
#    # scanType is the type of scan, PPI or RHI
#
#    radvar_vec=radvar.ravel()
#    
##==============================================================================
## 
##     # From polar to cartesian
##     ####################################################################################################################
#    radvar_cart=np.zeros(cart2pol['x'].shape)
#    numSamplesPerGridPoint=np.zeros(cart2pol['x'].shape)
## 
#    for vec in pol2cart:
#        numSamplesPerGridPoint[vec[0],vec[1]]=numSamplesPerGridPoint[vec[0],vec[1]]+1
#        radvar_cart[vec[0],vec[1]]=radvar_cart[vec[0],vec[1]]+radvar_vec[vec[2]]
## 
#    numSamplesPerGridPoint[numSamplesPerGridPoint==0]=float('nan')
#    radvar_cart=radvar_cart/numSamplesPerGridPoint
##==============================================================================
#
#    # From cartesian to polar
#    ####################################################################################################################
#    # Now assign the value of the closest radar bin to the missing pixels
#
#    N,M = cart2pol['x'].shape
#    
#    for i in range(0,N):
#       for j in range(0,M):
#           idx=[int(cart2pol['x'][i,j]),int(cart2pol['y'][i,j])]
#           if idx != [-1, -1]:
#               radvar_cart[i,j]=radvar[idx[0],idx[1]]
#    
#    return radvar_cart
#    
#    
#def radarBinInterpolation( p1, p2, x1, x2, radvar, X1_vec, X2_vec, scanType):
#    # Interpolates radar variable radvar sampled at polar coordinates p1, p2 that correspond to cartesian coordinates
#    # x1, x2 to a regular cartesian grid given by vectors X1_vec and X2_vec
#    # scanType is the type of scan, PPI or RHI
#
#    radvar_vec=radvar.ravel()
#    # Get interpolation, radial and angular resolutions
#    resInterp=X1_vec[1]-X1_vec[0]
#    resRange=p1[1]-p1[0]
#
#    # Get limits of coordinates
#    x1_0=X1_vec[0]
#    x2_0=X2_vec[0]
#    x1_end=X1_vec[-1]
#    x2_end=X2_vec[-1]
#
#    p1_0=p1[0]
#    p2_0=p2[0]
#    p1_end=p1[-1]
#    p2_end=p2[-1]
#
#    coordsBins=np.vstack((x1,x2)).T
#
#    # From polar to cartesian
#    ####################################################################################################################
#    radvar_cart=np.zeros((len(X2_vec), len(X1_vec)))
#
#    i=0
#    numSamplesPerGridPoint=np.zeros((len(X2_vec), len(X1_vec)))
#
#    pol2cart=[]
#    
#    for coord in coordsBins:
#        if(coord[0]>=x1_0 and coord[0]<=x1_end and coord[1]>=x2_0 and coord[1]<=x2_end):
#            idx=(np.floor((coord[1]-x2_0)/resInterp), np.floor((coord[0]-x1_0)/resInterp))
#            numSamplesPerGridPoint[idx]=numSamplesPerGridPoint[idx]+1
#
#            radvar_cart[idx]=radvar_cart[idx]+radvar_vec[i]
#            pol2cart.append([idx[0], idx[1],i])
#        i=i+1
#
#    numSamplesPerGridPoint[numSamplesPerGridPoint==0]=float('nan')
#    radvar_cart=radvar_cart/numSamplesPerGridPoint
#    
#    # From cartesian to polar
#    ####################################################################################################################
#    # Now assign the value of the closest radar bin to the missing pixels
#    [X1_Grid, X2_Grid] = np.meshgrid(X1_vec, X2_vec)
#
#    Rgrid=np.sqrt(X1_Grid**2+X2_Grid**2)
#    if(scanType=='RHI'):
#      AngGrid=-(np.arctan2(X1_Grid, X2_Grid)-np.pi/2)*1/DEG2RAD
#    elif(scanType=='PPI'):
#      AngGrid=(np.arctan2(-X1_Grid, -X2_Grid)+np.pi)*1/DEG2RAD
#    else:
#      print 'INVALID SCAN TYPE!!'
#      return []
#    # Find out which angles are not used
#    if(p2_0<p2_end):
#      anglesUsed=np.logical_and(AngGrid>=p2_0,AngGrid<=p2_end)
#    if(p2_0>p2_end):
#      anglesUsed=np.logical_or(AngGrid>=p2_0,AngGrid<=p2_end)
#    N=len(X2_vec)
#    M=len(X1_vec)
#    
#    cart2pol={}
#    cart2pol['x']=np.ones((N,M))*-1
#    cart2pol['y']=np.ones((N,M))*-1
#
#    
#    for i in range(0,N):
#       for j in range(0,M):
#            if(anglesUsed[i,j]):
#              if(np.isnan(radvar_cart[i,j])):
#                  if(Rgrid[i,j]>=p1_0 and Rgrid[i,j]<=p1_end and checkAng(AngGrid[i,j],p2_0,p2_end)):
#                      idx=(np.floor((Rgrid[i,j]-p1_0)/resRange),np.argmin(abs(AngGrid[i,j]-p2)))
#                      radvar_cart[i,j]=radvar[idx]
#                      cart2pol['x'][i,j]=idx[0]
#                      cart2pol['y'][i,j]=idx[1]
#                    
#    return radvar_cart, pol2cart, cart2pol
#
#def checkAng(ang, ang_0, ang_end):
#  out=False
#  if(ang_0<ang_end):
#    if(ang>ang_0 and ang<=ang_end):
#      return True
#  else:
#    if(ang>ang_0):
#      return True
#    elif(ang<ang_end):
#      return True
#
#

#
#def getRHI_polar(filenames, varnames, azimuth, maxRange=-1):
#    RHI_pol={}
#    maxLen=0
#    
#    lists={}
#    for var in varnames:
#        lists[var]=[]
#        
#    for idx, fname in enumerate(filenames):
#        f=h5py.File(fname,'r')
#        for var in varnames:
#            data=f['moments'][var][azimuth-1,:]
#            data=np.asarray(data)
#            data=data.astype(float)
#            clut=f['moments']['CLUT'][azimuth-1,:]
#            data[clut>100]=float('nan')
#            if maxRange>0:
#                data=data[0:np.floor(maxRange/RAD_RES)]
#            lists[var].append(int2float_radar(data, var, idx))
#            leng=len(lists[var][idx])
#            maxLen=np.max([maxLen,leng])
#    
#    for var in varnames:
#        RHI_pol[var]=np.zeros((maxLen, len(lists[var])))*float('nan')
#        for i in range(0, len(lists[var])):
#            RHI_pol[var][0:len(lists[var][i]),i]=lists[var][i]
#    
#    RHI_pol['elevations']=np.asarray(ELEVATION_ANGLES)
#    RHI_pol['azimuth']=azimuth
#    RHI_pol['range']=np.asarray(np.arange(0,maxLen)*RAD_RES)
#    
#    return RHI_pol
#
#
#def readRadData(filename, variableList=[]):
#    # This function reads a netcdf containing processed radar data in polar coordinates
#    # Input  : filename (complete path of the file
#    #          variableList, list of variables to be projected
#    # Output : varPol, dictionary containing the projected variables, the azimuth and the range
#
#    varPol={}
#
#    ncid=netCDF4.Dataset(filename,'r')
#
#    # Check if PPI or RHI and get azimuth or elevation
#    if('PPI' in filename):
#       elevation=np.mean(ncid.variables['Elevation'][:])
#       azimuth = ncid.variables['Azimuth'][:]
#    elif('RHI' in filename):
#       elevation=ncid.variables['Elevation'][:]
#       azimuth = np.mean(ncid.variables['Azimuth'][:])
#
#    # Get range
#    range = ncid.variables['Range'][:]
#    # Get variables in polar coordinates
#    for varname in variableList:
#        data=[]
#
#        if('PPI' in filename):
#            data=ncid.variables[varname][:]
#        elif('RHI' in filename):
#            data=ncid.variables[varname][:]
#
#        data[data==-9.99000000e+04]=float('nan')
#        varPol[varname]=data
#
#    varPol['coordRad']=(ncid.__getattr__('Longitude-value'), ncid.__getattr__('Latitude-value'))
#    varPol['range']=range
#    varPol['azimuth']=azimuth
#    varPol['elevation']=elevation
#
#    # Close netcdf
#    ncid.close()
#
#    return varPol
#
#def sph2cart(az,elev,r):
#    z = r * np.sin(elev)
#    rcoselev = r * np.cos(elev)
#    x = rcoselev * np.sin(az)
#    y = rcoselev * np.cos(az)
#
#    return x,y,z




