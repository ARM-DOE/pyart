#! /usr/bin/env python
"""
Extract a single ray from a C-SAPR PPI volume and store in a .npz file.

"""

import numpy as np
import pyart

# read in the full radar
radar = pyart.io.read('095636.mdv')

# extract ray 191 from the requested fields
fields_to_extract = ['reflectivity', 'normalized_coherent_power',
                     'cross_correlation_ratio', 'specific_differential_phase',
                     'differential_phase', 'differential_reflectivity']
ray_dic = dict()
for field_name in fields_to_extract:
    ray = radar.fields[field_name]['data'].data[191]
    ray = ray.reshape(1, -1)
    ray_dic[field_name] = ray

# save into a npz file
np.savez('example_rays.npz', **ray_dic)
