#! /usr/bin/env python
"""
Make a small MDV file containing a grid with X and Y axes in degrees.

Partial grid is taken from full sized file 000000.mdv which is contained in
the 200202.MASTER15.mdv.tar.gz file available online at:
http://www2.mmm.ucar.edu/imagearchive/WSI/mdv/
"""
# The MDV RLE decoding routine ends at the end of the data stream even if not
# all data point have been read.  Therefore a file truncated at a break in the
# RLE can still be read but data past the last data point will be filled with
# random data.  Here we end after the first non-keyed RLE byte.
infile = open('000000.mdv', 'rb')
outfile = open('example_mdv_grid.mdv', 'wb')
outfile.write(infile.read(8134))
infile.close()
outfile.close()
