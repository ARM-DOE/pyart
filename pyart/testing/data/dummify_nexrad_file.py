# /usr/bin/env python
# dummify (replace all data with a single value) a NEXRAD Level II file.

import numpy as np
import pyart.io.nexrad_level2 as nexrad
import struct

NEXRAD_FILE = 'KATX20130717_195021_V06'
OUTPUT_FILE = 'KATX20130717_195021_V06_DUMMY'

# sizes of structures
RECORD_SIZE = 2432
MSG_HEADER_SIZE = 16
DATA_BLOCK_SIZE = 28

fin = open(NEXRAD_FILE, 'rb')
out = open(OUTPUT_FILE, 'wb')

# read in the first 134 records as write as is
out.write(fin.read(24 + 12 + RECORD_SIZE * 134))

# read the rest of the file, and create a array of int8
buf = fin.read()
buf_length = len(buf)
buf_data = np.frombuffer(buf, dtype='i1')
fin.close()

# replace all radial data with a single value so it compresses well
pos = 0
while pos < buf_length:

    # retrieve the record (msg31)
    new_pos, record = nexrad._get_record_from_buf(buf, pos)
    if record['header']['type'] != 31:
        pos += RECORD_SIZE
        continue

    # dummy out all moments
    for bp_number in [4, 5, 6, 7, 8, 9]:
        bp_str = 'block_pointer_' + str(bp_number)
        block_pointer = record['msg_header'][bp_str]
        if block_pointer != 0:
            bpos = pos + block_pointer + MSG_HEADER_SIZE
            d = nexrad._unpack_from_buf(buf, bpos, nexrad.GENERIC_DATA_BLOCK)
            points = d['ngates']
            if d['data_name'] == 'PHI':
                points *= 2
            start = pos + MSG_HEADER_SIZE + DATA_BLOCK_SIZE + block_pointer
            buf_data[start:start+points] = 2

    pos = new_pos

# write out the modified data
out.write(buf_data.tostring())
out.close()
