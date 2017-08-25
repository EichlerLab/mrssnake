import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.uint16

ctypedef np.uint16_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)

def create_depth_array(unsigned short [:, :] count, unsigned short [:, :] depth, unsigned int rlen=36):
    cdef unsigned int i, j, k
    for i in range(count.shape[0] - rlen):
        for j in range(count.shape[1]):
            if count[i, j] != 0:
                for k in range(rlen):
                    depth[i+k, j] += count[i, j]
