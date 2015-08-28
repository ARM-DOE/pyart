cimport cython


@cython.boundscheck(False)
def _load_nn_field_data(object[:, :] data, int nfields, int npoints,
                        int[:] r_nums, int[:] e_nums, double[:, :] sdata):
    """
    _load_nn_field_data(data, nfields, npoints, r_nums, e_nums, sdata)

    Load the nearest neighbor field data into sdata
    """
    cdef unsigned int i, j, r_num, e_num

    for i in range(npoints):
        r_num = r_nums[i]
        e_num = e_nums[i]
        for j in range(nfields):
            # if we knew the dtype of data[j, r_num] we could speed this
            # up with a memory view, but we can't guarantee the dtype without
            # making a copy
            sdata[i, j] = data[j, r_num][e_num]
    return
