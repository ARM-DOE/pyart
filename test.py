
# test file KBGM 10:30 UTC 7 feb 2020
import pyart
import nexradaws
import numpy as np
from scipy import ndimage
import itertools
import matplotlib.pyplot as plt
import time
import scipy

conn = nexradaws.NexradAwsInterface()
availscans = conn.get_avail_scans('2020', '02', '07', 'KBGM')
availscans = conn.get_avail_scans('2021', '02', '07', 'KOKX')

scan = availscans[120]
scan = availscans[122]
path = scan.create_filepath(basepath='s3://noaa-nexrad-level2/', keep_aws_structure=True)[1]
path = path.replace("\\", '/')
radar = pyart.io.read_nexrad_archive(path)

# extract the lowest
radar = radar.extract_sweeps([0])

# interpolate to grid
grid = pyart.map.grid_from_radars(
    (radar,), grid_shape=(1, 251, 251),
    grid_limits=((0, 10000), (-200000.0, 200000.0), (-200000.0, 200000.0)),
    fields=['reflectivity'])

dx = grid.x['data'][1] - grid.x['data'][0]
dy = grid.y['data'][1] - grid.y['data'][0]

refl_og = grid.fields['reflectivity']['data'][0,:,:]
refl_og = np.ma.masked_invalid(refl_og)
refl_og = np.ma.masked_less(refl_og, 5)

refl_linear = 10 ** (refl_og / 10)  # mm6/m3

snow_rate = (refl_linear/57.3)**(1/1.67)
refl = snow_rate

refl = refl_og


# trying scipy filter
refl_bkg = ndimage.uniform_filter(refl, size=11) # doesn't work because of NaNs
refl_bkg = ndimage.generic_filter(refl, np.nanmean, mode='constant', cval=np.nan, size=11)
refl_bkg = np.ma.masked_where(refl.mask, refl_bkg)

bkg_mask_bool = bkg_mask_array.astype(bool)
time1 = time.time()
refl_bkg_scipy = ndimage.generic_filter(refl.filled(np.nan), function=np.nanmean, mode='constant',
                                        footprint=bkg_mask_bool, cval=np.nan)
refl_bkg_scipy = np.ma.masked_where(refl.mask, refl_bkg_scipy)
time2 = time.time()
print(time2-time1)

# Create idealized arrays
refl_test = np.empty_like(refl)
refl_test[:] = 30
refl_test[115:125,115:125] = 50

# create idealized array (gradient grid with gradient background)
# create background gradient
#bkg_gradient = np.linspace(0, 50, np.shape(ref)[0])
bkg_gradient = np.linspace(0, 50, np.shape(refl)[0])

# block background
bkg_gradient = 5 * np.round(bkg_gradient/5)

# constant background
bkg_gradient = np.empty_like(bkg_gradient)
bkg_gradient[:] = 25
# repeat array over other dimension
refl_test = np.array(np.shape(refl)[1]*[bkg_gradient])
# given a 240 x 240 array, let's set the cores (5x5)
core_gradient = np.linspace(25, 50, 5)
dn = int(np.floor(np.shape(refl_test)[0]/5))
indices = np.arange(int(dn/2), int(dn/2+5*dn), dn)
for ind in itertools.product(indices, indices):
    core_ind = np.where(indices==ind[0])[0][0]
    refl_test[ind[0]-3:ind[0]+3, ind[1]-3:ind[1]+3] = core_gradient[core_ind]
    #ref_test[ind[0], ind[1]] = core_gradient[core_ind]
    #print(ind[0])

refl = refl_test.T
refl = np.ma.masked_invalid(refl)

# plot
fig = plt.figure()
ax = plt.axes()
cs = ax.pcolormesh(refl, vmin=0, vmax=50)
cbar = plt.colorbar(cs)
ax.set_aspect('equal')
plt.show()

fig = plt.figure()
ax = plt.axes()
cs = ax.pcolormesh(ze_bkg, vmin=0, vmax=50)
cbar = plt.colorbar(cs)
ax.set_aspect('equal')
#fig.savefig("Q:\\My Drive\\phd\\winter_storms\\test.png", dpi=800)
plt.show()

fig = plt.figure()
ax = plt.axes()
cs = ax.pcolormesh(conv_core_array, vmin=0, vmax=3)
cbar = plt.colorbar(cs)
ax.set_aspect('equal')
plt.show()

fig = plt.figure()
ax = plt.axes()
cs = ax.pcolormesh(conv_strat_array, vmin=0, vmax=2)
cbar = plt.colorbar(cs)
ax.set_aspect('equal')
plt.show()

fig = plt.figure()
ax = plt.axes()
cs = ax.pcolormesh(convRadiuskm, vmin=1, vmax=5)
cbar = plt.colorbar(cs)
ax.set_aspect('equal')
plt.show()


fig = plt.figure()
ax = plt.axes()
cs = ax.pcolormesh(zeDiff, vmin=0, vmax=50)
cbar = plt.colorbar(cs)
ax.set_aspect('equal')
#fig.savefig("Q:\\My Drive\\phd\\winter_storms\\test.png", dpi=800)
plt.show()

fig = plt.figure()
ax = plt.axes()
cs = ax.pcolormesh(refl, vmin=0, vmax=12)
cbar = plt.colorbar(cs)
ax.set_aspect('equal')
plt.show()

fig = plt.figure()
ax = plt.axes()
cs = ax.pcolormesh(ze_bkg_test, vmin=0, vmax=10)#, cmap='magma_r')
cbar = plt.colorbar(cs)
ax.set_aspect('equal')
fig.savefig("Q:\\My Drive\\phd\\winter_storms\\nexrad_stitching\\convsf_imgs\\testing\\ze_bkg_10km_kbgm.png", dpi=600, bbox_inches='tight')
#fig.savefig("Q:\\My Drive\\phd\\winter_storms\\test.png", dpi=800)
plt.show()

fig = plt.figure()
ax = plt.axes()
cs = ax.pcolormesh(conv_core_array, vmin=0, vmax=3)
cbar = plt.colorbar(cs)
ax.set_aspect('equal')
plt.show()

# plot
fig = plt.figure()
ax = plt.axes()
cs = ax.pcolormesh(snow_rate, vmin=0, vmax=10)
cbar = plt.colorbar(cs)
ax.set_aspect('equal')
plt.show()

fig = plt.figure()
ax = plt.axes()
cs = ax.pcolormesh(refl_bkg_scipy, vmin=0, vmax=10)
cbar = plt.colorbar(cs)
ax.set_aspect('equal')
plt.show()

fig = plt.figure()
ax = plt.axes()
cs = ax.pcolormesh(bkg_diff, vmin=-5e-7, vmax=5e-7, cmap='RdBu_r')
cbar = plt.colorbar(cs)
ax.set_aspect('equal')
plt.show()


