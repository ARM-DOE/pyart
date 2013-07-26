#! /usr/bin/env python
"""
Extract a single ray from a C-SAPR PPI volume and store in a .npz file.

"""

import numpy as np
import pyart

# read in the full radar
radar = pyart.io.read_mdv('095636.mdv')

# extract ray 191 from the requested fields
fields_to_extract = ['reflectivity_horizontal', 'norm_coherent_power',
                     'copol_coeff', 'dp_phase_shift', 'diff_phase']
ray_dic = dict()
for field_name in fields_to_extract:
    ray = radar.fields[field_name]['data'].data[191]
    ray = ray.reshape(1, -1)
    ray_dic[field_name] = ray

# save into a npz file
np.savez('example_rays.npz', **ray_dic)
