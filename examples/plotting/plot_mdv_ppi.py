"""
==========================
Plot a PPI from a MDV file
==========================

This is an example showing how to use Py-ART to plot a field from a MDV file.

"""

# Author: Jonathan J. Helmus (jhelmus@anl.gov)
# License: BSD 3 clause

import matplotlib.pyplot as plt
import pyart

filename = '110635.mdv'

# read in the data, create the display
mdvfile = pyart.io.mdv.MdvFile(filename)
display = pyart.graph.MdvDisplay(mdvfile)

# create the figure
fig = plt.figure(figsize=[5, 5])
ax = fig.add_subplot(111, frameon=False)
display.plot_ppi('DBZ_F', 0, mask_tuple=['NCP_F', 0.5],
                 colorbar_flag=False, title_flag=False,
                 axislabels_flag=False)
display.set_limits(ylim=[-120, 120], xlim=[-120, 120])
fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
fig.show()
