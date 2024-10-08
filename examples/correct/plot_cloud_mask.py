"""
=====================================
Calculating and Plotting a Cloud Mask
=====================================

This example shows how to correct and plot reflectivity from an ARM
KAZR using a noise floor cloud mask.

"""
print(__doc__)

# Author: Adam Theisen and Zach Sherman
# License: BSD 3 clause

import matplotlib.pyplot as plt
import numpy as np
from open_radar_data import DATASETS

import pyart

############################
# **Read and plot raw data**
#
# First let's read and plot our dataset without any mask.

# Fetch and read in the ARM KAZR file.
filename = DATASETS.fetch("sgpkazrgeC1.a1.20190529.000002.cdf")
radar = pyart.aux_io.read_kazr(filename)

# Let's now take a look at reflectivity data prior to any corrections.
display = pyart.graph.RadarDisplay(radar)
display.plot("reflectivity_copol")
display.set_limits(xlim=(0, 55))
plt.show()

#################################################
# **Calculate cloud mask and plot corrected data**
#
# Now lets apply a mask by using the calc_cloud_mask function
# that will use a noise floor calculation from range and more
# to calculate the mask.

# First lets correct the data by calculating the mask.
cloud_mask_radar = pyart.correct.calc_cloud_mask(radar, "reflectivity_copol", "range")

# In this new radar object we should now have a new cloud mask field.
print(cloud_mask_radar.fields["cloud_mask_2"])

# Next we'll create a copy of the reflectivity field so we are not
# overwriting the original data.
cloud_mask_radar.add_field_like(
    "reflectivity_copol",
    "reflectivity_cloud_mask",
    cloud_mask_radar.fields["reflectivity_copol"]["data"].copy(),
)

# Now let's apply the mask to the copied reflectivity data.
cloud_mask_radar.fields["reflectivity_cloud_mask"]["data"][
    cloud_mask_radar.fields["cloud_mask_2"]["data"] == 0
] = np.nan

# And now we can plot the masked reflectivity field.
display = pyart.graph.RadarDisplay(cloud_mask_radar)
display.plot("reflectivity_cloud_mask")
display.set_limits(xlim=(0, 55))
plt.show()
