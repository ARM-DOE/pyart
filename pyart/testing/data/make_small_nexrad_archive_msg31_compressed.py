#! /usr/bin/env python
"""
Make a small NEXRAD (WSR-88D) Level II file with compressed records.

The file created by this script only contains the first 254 records
(120 radials) of the full sized file.  It should not be used as a full sized
example, only for testing the decompression routines.
"""

f = open('Level2_KATX_20130717_1950.ar2v', 'rb')
o = open('example_nexrad_archive_msg31_compressed.ar2v', 'wb')

o.write(f.read(24))             # volume header
o.write(f.read(12527 + 4))      # first series of compressed records
o.write(f.read(105727 + 4))     # second series of compressed records
o.close()
