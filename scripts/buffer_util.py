import numpy as np, multiprocessing as mp


def get_np_view_of_buffer(buf, dtype, shape):
    '''
        given a piece of memory, view it as a numpy matrix
        of given dtype and shape  
    '''

    return np.frombuffer(buf, dtype=dtype).reshape(shape)


def allocate_shared_buffer(_dtype, shape):
    '''
        returns pointer to allocated, uninitiated piece of
        memory to be shared among processes
    '''
    dtype = np.dtype(_dtype)
    cdtype = np.ctypeslib.as_ctypes_type(dtype)
    shared_buffer = mp.RawArray(cdtype, shape[0]*shape[1])

    return shared_buffer
