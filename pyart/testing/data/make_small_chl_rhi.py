#! /usr/bin/env python
"""
Make a small CSU-CHILL, CHL file containing the first ray in each sweep
from a two scan RHI volume.

All field data is left untouched.
"""

import struct
import numpy as np

# parameters
INFILE = 'CHL20120705_230123'
OUTFILE = 'example_chl_rhi.chl'

# open the input and output files
f = open(INFILE, 'rb')
out = open(OUTFILE, 'wb')

# FILE_HDR block (56 bytes)
out.write(f.read(56))

# 30 FIELD_SCALE block (232 bytes echo, 6960 bytes total)
out.write(f.read(6960))

# RADAR_INFO (128 bytes)
out.write(f.read(128))

# PROCESSOR_INFO (88 bytes)
out.write(f.read(88))

# UNKNOWN block (84 bytes)
out.write(f.read(84))

# SCAN_SEG block (140 bytes)
out.write(f.read(140))

# UNKNOWN block (128 bytes)
out.write(f.read(128))

# RAY_HDR and data (64056 bytes)
out.write(f.read(64056))

# skip to second scan
f.seek(2841480)

# RADAR_INFO block
out.write(f.read(128))

# SCAN_SEG block
out.write(f.read(140))

# UNKNOWN block
out.write(f.read(16))

# UNKNOWN block
out.write(f.read(128))

# UNKNOWN block
out.write(f.read(2072))

# RAY_HDR and data (64056 bytes)
out.write(f.read(64056))

# skip to SWEEP_BLOCK (44 bytes at end of file)
f.seek(5677548)
out.write(f.read(44))

# close the files
f.close()
out.close()
