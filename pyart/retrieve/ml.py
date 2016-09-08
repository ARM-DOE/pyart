# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 11:18:25 2016

@author: wolfensb
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, signal
from scipy.ndimage.filters import convolve

# Constants
MAXTHICKNESS_ML = 1000
LOWMLBOUND = 0.7
UPMLBOUND = 1.3
SIZEFILT_M = 75
THRESHLENGTH = 250 # Interpolate holes of 250 m max


def normalizeImage(im,minval,maxval):
    newMax=1
    newMin=0

    out=(im-minval)*(newMax-newMin)/(maxval-minval)+newMin
    out[im>maxval]=newMax
    out[im<minval]=newMin

    return out

def gradient2D(im):
    # dim = 1 = gradient along the rows (Y)
    # dim = 2 = gradient along the column (X)

    Gx=np.array([[1,0,-1],[2, 0, -2],[1, 0, -1]])
    Gy=np.array([[1,2,1],[0, 0, 0],[-1, -2, -1]])

    out={}

    out['Gx']=convolve(im,Gx,  mode='reflect')
    out['Gy']=convolve(im,Gy,  mode='reflect')

    return out
    
def convolve_with_nan(input_array,kernel,boundary='mirror'):
    if isinstance(input_array,np.ma.masked_array):
        input_array = np.ma.masked(input_array,np.nan)

    # Flat function with NaNs for comparison.
    on = np.ones(input_array.shape)
    
    # Find all the NaNs in the input.
    n = np.isnan(a)
    
    # Replace NaNs with zero, both in 'a' and 'on'.
    a[n] = 0.
    on[n] = 0.
    
    # Check that the filter does not have NaNs.
    if(np.isnan(kernel).any()):
        print('Kernel contains NaN values.')
        print('Aborting...')
        return None
    
    # Calculate what a 'flat' function looks like after convolution.
    flat = convolve(on,kernel,mode=boundary)
    # Do the actual convolution
    c = convolve(a,kernel,mode=boundary)/flat;
    return c    

def mean_filter(input_array, shape=(3,3), boundary='extend'):
    kernel=np.ones(shape)
    kernel=kernel/np.sum(kernel.ravel())
    out=convolve_with_nan(input_array, kernel, boundary)

    return out
    
def detectML(radInst, threshold=0.02, interp=True,checkMinLength=True):
    resolution=radInst.attributes['interp_res']
    if(resolution!=25 and resolution!=75):
      print('[ERROR] The ML detection algorithm is designed to work for a Cartesian grid of 25x25 or 75x75 cell size, please choose the' \
            ' resolution accordingly when you create a class instance')
      output={}
      output['validMLFlag']= False
      return output

    imZh=normalizeImage(radInst.radar_variables['Zh'],(10,60))
    imRhohv=normalizeImage(radInst.radar_variables['Rhohv'],(0.75,1))
    # Combine
    ImageInput=(1-imRhohv)*imZh

    ImageInput[np.isnan(ImageInput)]=0

    sizeFilt=np.floor(SIZEFILT_M/resolution)

    gradient=gradient2D(mean_filter(ImageInput, (sizeFilt,sizeFilt)))

    gradientZ=gradient['Gy']

    bottomML,topML=processMapML(gradientZ,  radInst.radar_variables, threshold)

    medianBottomHeight=np.median(bottomML[~np.isnan(bottomML)])
    medianTopHeight=np.median(topML[~np.isnan(topML)])

    if(not np.isnan(medianBottomHeight)):
        gradientZ[0:np.floor(LOWMLBOUND*medianBottomHeight),:]=float('nan')
    if(not np.isnan(medianTopHeight)):
        gradientZ[np.floor(UPMLBOUND*medianTopHeight):gradientZ.shape[0]-1,:]=float('nan')

    bottomML,topML=processMapML(gradientZ, radInst.radar_variables, threshold, (0.65,1))


    mapML=np.zeros(gradientZ.shape)

    # New part
    ###################################################################
    # Compute gradient of Zh only
    imZh[np.isnan(radInst.radar_variables['Zh'])]=float('nan')

    gradient_Zh=gradient2D(mean_filter(imZh,(sizeFilt,sizeFilt)))

    gradientZ_Zh=gradient_Zh['Gy']

    gradientZ_Zh[np.isnan(imZh)]=float('nan')

    # We need to have an idea of the top of the ML everywhere so we interpolate
    # linearly
    if(sum(~np.isnan(topML))>10): # Check if boundaries of a ML were detected
      topML_fill=interpolate.InterpolatedUnivariateSpline(np.array(np.where(~np.isnan(topML))[0]),topML[~np.isnan(topML)])(np.array(np.where(np.isnan(topML))[0]))
      topML_interp=topML
      topML_interp[np.isnan(topML_interp)]=topML_fill

      for i in range(0,len(topML)-1):
        gradientZ_Zh[0:topML_interp[i],i]=float('nan')
        # We cut the gradient above as soon as we go from lower to higher
        # values (gradient line inflexion point)

        _max = signal.find_peaks_cwt(gradientZ_Zh[:,i],np.arange(1,6))

        if(_max):
          loc=_max[0]
          gradientZ_Zh[loc:-1,i]=float('nan')

      # Threshold gradientZ
      gradientZ_Zh[np.abs(gradientZ_Zh)<=threshold]=float('nan')

      topML=processMapML_onlyZh(gradientZ_Zh)

      medianBottomHeight=np.median(bottomML[~np.isnan(bottomML)])
      medianTopHeight=np.median(topML[~np.isnan(topML)])

      if(~np.isnan(medianTopHeight)):
        gradientZ_Zh[np.floor(UPMLBOUND*medianTopHeight):-1,:]=float('nan')

      topML=processMapML_onlyZh(gradientZ_Zh)

      # Put NaN at pixels which have a ML thickness larger than max_thickness  (1000)
      # Also put NaN where either the top or the bottom are undefined

      bad_pixels=np.where((topML-bottomML)>MAXTHICKNESS_ML/resolution)
      topML[bad_pixels]=float('nan')
      bottomML[bad_pixels]=float('nan')
      topML[np.isnan(bottomML)]=float('nan')
      bottomML[np.isnan(topML)]=float('nan')

      medianBottomHeight=np.median(bottomML[~np.isnan(bottomML)])
      medianTopHeight=np.median(topML[~np.isnan(topML)])

    # End of new part
    ###################################################################

    if(interp and np.sum(~np.isnan(bottomML))>=2):
        # Calculate subsequences
        sub=calcSubInd(bottomML)
        sub2interp={}
        sub2interp['lengths']=sub['lengths'][np.logical_and(sub['values']==0,sub['lengths']<=THRESHLENGTH/resolution)]
        sub2interp['idx']=sub['idx'][np.logical_and(sub['values']==0 , sub['lengths']<=THRESHLENGTH/resolution)]
        index2interp=[]
        for k in range(0,len(sub2interp['idx'])):
            index2interp.extend(range(sub2interp['idx'][k], sub2interp['idx'][k]+sub2interp['lengths'][k]))
            # Interpolate all subsequences of less than threshLength [m] using
            # piecewise cubic hermite interpolation
        index2interp=np.array(index2interp)
        if(index2interp.size):
            bottomML[index2interp]= interpolate.pchip(np.array(np.where(~np.isnan(bottomML))[0]),np.array(bottomML[~np.isnan(bottomML)]))(index2interp)
            topML[index2interp]= interpolate.pchip(np.array(np.where(~np.isnan(topML))[0]),topML[~np.isnan(topML)])(index2interp)

    # Finds the missing subsequences and calculate their characteristics
    sub=calcSubInd(bottomML)

    subSeqNaN=sub['lengths'][sub['values']==False]
    # In the calculation of the number of NaNs remove the ones due to the removal of the high elevation angles
    midML=(medianTopHeight+medianBottomHeight)/2
    numNaN=sum(subSeqNaN)-midML*np.tan(radInst.attributes['maxElev']/180*np.pi)

    # Check if ML is valid
    # 1) check if medianBottomHeight and medianTopHeight are defined
    if(np.isnan(medianBottomHeight+medianTopHeight)):
        invalidML=True
    else:
        invalidML=False
        # 2) Check how many values in the data are defined at the height of the ML
        lineVal=radInst.radar_variables['Rhohv'][np.round(midML),:]
        # Check if ML is long enough
        if(checkMinLength):
            if(numNaN>0.5*sum(~np.isnan(lineVal))):
                invalidML=True


    # If ML is invalid, just fill topML and bottomML with NaNs
    if(invalidML):
        topML=float('nan')*np.zeros((1,gradientZ.shape[1]))
        bottomML=topML
    else:
        for j in range(0,len(topML)-1):
            if(not np.isnan(topML[j]) and not np.isnan(bottomML[j])):
                mapML[np.round(bottomML[j]):np.round(topML[j]),j]=1

    output={}
    output['validMLFlag']= not invalidML
    output['MapML']=np.array(mapML)
    output['BottomML']=np.array((bottomML)*resolution)
    output['TopML']=np.array((topML)*resolution)

    return output

def processMapML(gradientZ, radar_variables, threshold ,thresholdRhohvMin = 0, threshRhohvMax = np.Inf):
    Ncols=gradientZ.shape[1]
    bottomML=np.zeros((Ncols))*float('nan')
    topML=np.zeros((Ncols))*float('nan')

    for j in range(0,Ncols-1):
        gradLine=gradientZ[:,j]
        gradNoNan=gradLine
        gradNoNan=gradNoNan[~np.isnan(gradNoNan)]
        indBot=float('nan')
        indTop=float('nan')
        if(gradNoNan.size):

            sortedGrad=np.sort(gradNoNan)
            maxVal=sortedGrad[-1]
            minVal=sortedGrad[0]
            indTop=np.where(gradLine==minVal)
            indTop=indTop[0][0]+2
            indBot=np.where(gradLine==maxVal)
            indBot=indBot[0][0]

            if(not indBot.size or not indTop.size or indTop<=indBot):
                rhohvThreshCond=False
                rhohvNaNCond=False
            else:
                rhohvThreshCond=np.nanmax(radar_variables['Rhohv'][indBot:indTop,j])<=threshRhohvMax\
                and np.nanmin(radar_variables['Rhohv'][indBot:indTop,j])>=thresholdRhohvMin

                rhohvNaNCond=sum(np.isnan(radar_variables['Rhohv'][indBot:indTop,j]))==0

            if(minVal<=-threshold and maxVal >= threshold and rhohvNaNCond and rhohvThreshCond):
                bottomML[j]=indBot-1
                topML[j]=indTop+1

    return bottomML, topML

def processMapML_onlyZh(gradientZ):
    Ncols=gradientZ.shape[1]
    topML=np.zeros((Ncols))*float('nan')

    for j in range(0,Ncols-1):
        gradLine=gradientZ[:,j]
        gradNoNan=gradLine
        gradNoNan=gradNoNan[~np.isnan(gradNoNan)]
        if(gradNoNan.size):

            sortedGrad=np.sort(gradNoNan)
            minVal=sortedGrad[0]
            indTop=np.where(gradLine==minVal)
            indTop=indTop[0][0]
            topML[j]=indTop+1

    return topML


def calcSubInd(inputVec):
    sub={}
    sub['values']=[]
    sub['lengths']=[]
    sub['idx']=[]

    for l in range(0, len(inputVec)-1):
        if(l==0):
            sub['idx'].append(l)
            sub['values'].append(~np.isnan(inputVec[l]))
        if(~np.isnan(inputVec[l]) != sub['values'][-1]):
            sub['values'].append(~np.isnan(inputVec[l]))
            sub['lengths'].append(l-sub['idx'][-1])
            sub['idx'].append(l)
    sub['lengths'].append(l+1-sub['idx'][-1])

    sub['lengths']=np.array(sub['lengths'])
    sub['idx']=np.array(sub['idx'])
    sub['values']=np.array(sub['values'])
    return sub
    

if __name__ == '__main__':
    import numpy as np
    a = np.random.rand(100,100)
    a[10:20,10:20] = np.nan
    b = convolve_with_nan(a,np.ones((10,10))/100.)
#    b2 = convolve(a,np.ones((5,5))/25.)
    
#    plt.imshow(b)
#    plt.figure()
#    plt.imshow(b2)
    
    plt.figure()
    plt.plot(b[:,1])
    
    p = signal.find_peaks_cwt(-b[:,1],np.arange(1,6))