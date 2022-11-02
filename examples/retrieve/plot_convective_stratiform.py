"""
=======================================
Convective-Stratiform classification examples
=======================================
This example shows how to use the updated convective stratiform classifcation algorithm. We show 3 examples,
a summer convective example, an example from Hurricane Ian, and an example from a winter storm. 

"""

print(__doc__)

# Author: Laura Tomkins (lmtomkin@ncsu.edu)
# License: BSD 3 clause


import pyart
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs

######################################
# **Example with summer convection**
#


# read in file
filename = pyart.testing.get_test_data('swx_20120520_0641.nc')
radar = pyart.io.read(filename)

# extract the lowest sweep
radar = radar.extract_sweeps([0])

# interpolate to grid
grid = pyart.map.grid_from_radars(
    (radar,), grid_shape=(1, 201, 201),
    grid_limits=((0, 10000), (-50000.0, 50000.0), (-50000.0, 50000.0)),
    fields=['reflectivity_horizontal'])

# get dx dy
dx = grid.x['data'][1] - grid.x['data'][0]
dy = grid.y['data'][1] - grid.y['data'][0]

# convective stratiform classification
convsf_dict = pyart.retrieve.conv_strat(grid, dx, dy, refl_field='reflectivity_horizontal', always_core_thres=40,
                                        bkg_rad_km=20, use_cosine=True, max_diff=5, zero_diff_cos_val=55,
                                        weak_echo_thres=10, max_conv_rad_km=2)

# add to grid object
# mask zero values (no surface echo)
convsf_masked = np.ma.masked_equal(convsf_dict['convsf']['data'], 0)
# mask 3 values (weak echo)
convsf_masked = np.ma.masked_equal(convsf_masked, 3)
# add dimension to array to add to grid object
convsf_dict['convsf']['data'] = convsf_masked[None,:,:]
# add field
grid.add_field('convsf', convsf_dict['convsf'], replace_existing=True)

# create plot using GridMapDisplay
# plot variables
display = pyart.graph.GridMapDisplay(grid)
magma_r_cmap = plt.get_cmap('magma_r')
ref_cmap = mcolors.LinearSegmentedColormap.from_list('ref_cmap', magma_r_cmap(np.linspace(0, 0.9, magma_r_cmap.N)))
projection = ccrs.LambertConformal(central_latitude=radar.latitude['data'][0],
                                   central_longitude=radar.longitude['data'][0])
######################################
# You'll notice that the convective stratiform field has less data around the edges compared to the reflectivity
# field. This is because the function is designed to only compute the background radius where the footprint has 75%
# of data, so along the edges there is not enough data to fulfill this requirement. The footprint percentage can be
# changed using the variable, calc_thres.

# plot
plt.figure(figsize=(10,4))
ax1=plt.subplot(1,2,1, projection=projection)
display.plot_grid('reflectivity_horizontal', vmin=5, vmax=45, cmap=ref_cmap, projection=projection,
                  transform=ccrs.PlateCarree(), ax=ax1)
ax2=plt.subplot(1,2,2, projection=projection)
display.plot_grid('convsf', vmin=0, vmax=2, cmap=plt.get_cmap('viridis', 3), projection=projection, ax=ax2,
                  transform=ccrs.PlateCarree(), ticks=[1/3, 1, 5/3], ticklabs=['', 'Stratiform', 'Convective'])
plt.show()


######################################
# In addition to the default convective-stratiform classification, the function also returns an underestimate (
# convsf_under) and an overestimate (convsf_over) to take into consideration the uncertainty when choosing
# classification parameters. The under and overestimate use the same parameters, but vary the input field by a
# certain value (default is 5 dBZ). The estimation can be turned off (estimate_flag=False), but we recommend keeping
# it turned on.

# mask weak echo and no surface echo
convsf_masked = np.ma.masked_equal(convsf_dict['convsf']['data'], 0)
convsf_masked = np.ma.masked_equal(convsf_masked, 3)
convsf_dict['convsf']['data'] = convsf_masked
# underest.
convsf_masked = np.ma.masked_equal(convsf_dict['convsf_under']['data'], 0)
convsf_masked = np.ma.masked_equal(convsf_masked, 3)
convsf_dict['convsf_under']['data'] = convsf_masked
# overest.
convsf_masked = np.ma.masked_equal(convsf_dict['convsf_over']['data'], 0)
convsf_masked = np.ma.masked_equal(convsf_masked, 3)
convsf_dict['convsf_over']['data'] = convsf_masked

# Plot each estimation
plt.figure(figsize=(10,4))
ax1=plt.subplot(131)
ax1.pcolormesh(convsf_dict['convsf']['data'], vmin=0, vmax=2, cmap=plt.get_cmap('viridis', 3))
ax1.set_title('Best estimate')
ax1.set_aspect('equal')
ax2=plt.subplot(132)
ax2.pcolormesh(convsf_dict['convsf_under']['data'], vmin=0, vmax=2, cmap=plt.get_cmap('viridis', 3))
ax2.set_title('Underestimate')
ax2.set_aspect('equal')
ax3=plt.subplot(133)
ax3.pcolormesh(convsf_dict['convsf_over']['data'], vmin=0, vmax=2, cmap=plt.get_cmap('viridis', 3))
ax3.set_title('Overestimate')
ax3.set_aspect('equal')
plt.show()

######################################
# **Tropical example**
# Let's get a NEXRAD file from Hurricane Ian

# Read in file
nexrad_file = 's3://noaa-nexrad-level2/2022/09/28/KTBW/KTBW20220928_190142_V06'
radar = pyart.io.read_nexrad_archive(nexrad_file)

# extract the lowest sweep
radar = radar.extract_sweeps([0])

# interpolate to grid
grid = pyart.map.grid_from_radars(
    (radar,), grid_shape=(1, 201, 201),
    grid_limits=((0, 10000), (-200000.0, 200000.0), (-200000.0, 200000.0)),
    fields=['reflectivity'])

# get dx dy
dx = grid.x['data'][1] - grid.x['data'][0]
dy = grid.y['data'][1] - grid.y['data'][0]

# convective stratiform classification
convsf_dict = pyart.retrieve.conv_strat(grid, dx, dy, refl_field='reflectivity', always_core_thres=40,
                                        bkg_rad_km=20, use_cosine=True, max_diff=3, zero_diff_cos_val=55,
                                        weak_echo_thres=5, max_conv_rad_km=2)

# add to grid object
# mask zero values (no surface echo)
convsf_masked = np.ma.masked_equal(convsf_dict['convsf']['data'], 0)
# mask 3 values (weak echo)
convsf_masked = np.ma.masked_equal(convsf_masked, 3)
# add dimension to array to add to grid object
convsf_dict['convsf']['data'] = convsf_masked[None,:,:]
# add field
grid.add_field('convsf', convsf_dict['convsf'], replace_existing=True)

# create plot using GridMapDisplay
# plot variables
display = pyart.graph.GridMapDisplay(grid)
magma_r_cmap = plt.get_cmap('magma_r')
ref_cmap = mcolors.LinearSegmentedColormap.from_list('ref_cmap', magma_r_cmap(np.linspace(0, 0.9, magma_r_cmap.N)))
projection = ccrs.LambertConformal(central_latitude=radar.latitude['data'][0],
                                   central_longitude=radar.longitude['data'][0])
# plot
plt.figure(figsize=(10,4))
ax1=plt.subplot(1,2,1, projection=projection)
display.plot_grid('reflectivity', vmin=5, vmax=45, cmap=ref_cmap, projection=projection,
                  transform=ccrs.PlateCarree(), ax=ax1, axislabels_flag=False)
ax2=plt.subplot(1,2,2, projection=projection)
display.plot_grid('convsf', vmin=0, vmax=2, cmap=plt.get_cmap('viridis', 3), projection=projection,
                  axislabels_flag=False, transform=ccrs.PlateCarree(), ticks=[1/3, 1, 5/3],
                  ticklabs=['', 'Stratiform', 'Convective'], ax=ax2)
plt.show()

######################################
# **Winter storm example with image muting**
# Here is a final example of the convective stratiform classification using an example from a winter storm. Before
# doing the classification, we image mute the reflectivity to remove regions with melting or mixed precipitation. We
# then rescale the reflectivity to snow rate. We recommend using a rescaled reflectivity to do the classification,
# but if you do make sure to changed dB_averaging to False because this parameter is used to convert reflectivity to
# a linear value before averaging (set dB_averaging to True for reflectivity fields in dBZ units).

# Read in file
nexrad_file = 's3://noaa-nexrad-level2/2021/02/07/KOKX/KOKX20210207_161413_V06'
radar = pyart.io.read_nexrad_archive(nexrad_file)

# extract the lowest sweep
radar = radar.extract_sweeps([0])

# interpolate to grid
grid = pyart.map.grid_from_radars(
    (radar,), grid_shape=(1, 201, 201),
    grid_limits=((0, 10000), (-200000.0, 200000.0), (-200000.0, 200000.0)),
    fields=['reflectivity', 'cross_correlation_ratio'])

# image mute grid object
grid = pyart.util.image_mute_radar(grid, 'reflectivity', 'cross_correlation_ratio', 0.97, 20)

# convect non-muted reflectivity to snow rate
nonmuted_ref = grid.fields['nonmuted_reflectivity']['data'][0,:,:]
nonmuted_ref = np.ma.masked_invalid(nonmuted_ref)

nonmuted_ref_linear = 10 ** (nonmuted_ref / 10)  # mm6/m3
snow_rate = (nonmuted_ref_linear/57.3)**(1/1.67) #

# add to grid
snow_rate_dict = {
        'data': snow_rate[None,:,:],
        'standard_name': 'snow_rate',
        'long_name': 'Snow rate converted from linear reflectivity',
        'units': 'mm/hr',
        'valid_min': 0,
        'valid_max': 40500}
grid.add_field('snow_rate', snow_rate_dict, replace_existing=True)

# get dx dy
dx = grid.x['data'][1] - grid.x['data'][0]
dy = grid.y['data'][1] - grid.y['data'][0]

# convective stratiform classification
convsf_dict = pyart.retrieve.conv_strat(grid, dx, dy, refl_field='snow_rate', dB_averaging=False,always_core_thres=4,
                                        bkg_rad_km=40, use_cosine=True, max_diff=1.5, zero_diff_cos_val=5,
                                        weak_echo_thres=0, min_dBZ_used=0, max_conv_rad_km=1)

# add to grid object
# mask zero values (no surface echo)
convsf_masked = np.ma.masked_equal(convsf_dict['convsf']['data'], 0)
# mask 3 values (weak echo)
convsf_masked = np.ma.masked_equal(convsf_masked, 3)
# add dimension to array to add to grid object
convsf_dict['convsf']['data'] = convsf_masked[None,:,:]
# add field
grid.add_field('convsf', convsf_dict['convsf'], replace_existing=True)

# create plot using GridMapDisplay
# plot variables
display = pyart.graph.GridMapDisplay(grid)
magma_r_cmap = plt.get_cmap('magma_r')
ref_cmap = mcolors.LinearSegmentedColormap.from_list('ref_cmap', magma_r_cmap(np.linspace(0, 0.9, magma_r_cmap.N)))
projection = ccrs.LambertConformal(central_latitude=radar.latitude['data'][0],
                                   central_longitude=radar.longitude['data'][0])
# plot
plt.figure(figsize=(10,4))
ax1=plt.subplot(1,2,1, projection=projection)
display.plot_grid('snow_rate', vmin=0, vmax=10, cmap=plt.get_cmap('viridis'), projection=projection,
                  transform=ccrs.PlateCarree(), ax=ax1, axislabels_flag=False)
ax2=plt.subplot(1,2,2, projection=projection)
display.plot_grid('convsf', vmin=0, vmax=2, cmap=plt.get_cmap('viridis', 3), projection=projection,
                  axislabels_flag=False, transform=ccrs.PlateCarree(), ticks=[1/3, 1, 5/3],
                  ticklabs=['', 'Stratiform', 'Convective'], ax=ax2)
plt.show()

