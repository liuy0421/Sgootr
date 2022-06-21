'''
    error_correction.py <filtered_rc_cna.npz>    
'''
import sys, signal, numpy as np, multiprocessing as mp
from datetime import datetime


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


def _init(buf_methylated, buf_unmethylated, dtype, shape):
    '''
        make pointer to memory buffer shareable among workers
        and redirect SIGINT signals to master so that script
        is interruptable.
    '''
    global sb_methylated, sb_unmethylated

    sb_methylated = (buf_methylated, dtype, shape)
    sb_unmethylated = (buf_unmethylated, dtype, shape)
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def enumerate_jobs(N, M, C, e):
   
    entries = np.where(~np.isnan(N))

    return [(i, j, N[i,j], M[i,j], C[i,j], e) for \
            (i,j) in zip(entries[0], entries[1])]


def safe_log(x): 

    return np.log(x*(10**100)+10**(-300))    


def sequencing_error_correction(i,j,n,m,c,e):
    '''
        i: cell index
        j: site index
        n: number of methylated reads
        m: number of unmethylated reads
        c: copy number call at site j in cell i
        e: sequencing error parameter
    '''
    arr_n = get_np_view_of_buffer(*sb_methylated)
    arr_m = get_np_view_of_buffer(*sb_unmethylated)
    

    if c == 0:
        arr_n[i,j], arr_m[i,j] = np.nan, np.nan 
        return 1, 1

    L = [0 for gamma in range(c+1)]
    L[0] = safe_log((e**n)*((1-e)**m)) 
    L[c] = safe_log(((1-e)**n)*(e**m))

    for gamma in range(1,c):
        L[gamma] = safe_log(((gamma/c)**n)*((1-gamma/c)**m))

    status = L.index(max(L))

    if (status == 0 or status == c) and (n != 0 and m != 0):
        if n > m:
            arr_n[i,j], arr_m[i,j] = n+m, 0
            return 1, 0
        elif n < m:
            arr_n[i,j], arr_m[i,j] = 0, n+m
            return 1, 0
        else:
            arr_n[i,j], arr_m[i,j] = np.nan, np.nan
            return 1, 1

    arr_n[i,j], arr_m[i,j] = n,m
    
    return 0, 0



if __name__ == "__main__":

    f = open(snakemake.log[0],'w')
    sys.stderr = sys.stdout = f
    
    f.write('[{}] gmelin-larch is performing read error ' \
            'correction\n'.format(datetime.now())) 

    f.write('[{}] reading input\n'.format(datetime.now()))
    obj = np.load(snakemake.input[0], allow_pickle=True)
    _N, _M, _C = obj['n'], obj['m'], obj['cna']
    e = snakemake.params.e
    n_entries = np.sum(~np.isnan(_N))

    f.write('[{}] done reading input, enumerating ' \
            'jobs\n'.format(datetime.now()))
    jobs = enumerate_jobs(_N, _M, _C, e)
    f.write('[{}] joblist takes {} bytes\n'.format(datetime.now(), \
                                                sys.getsizeof(jobs)))

    f.write('[{}] allocating and initializing memory to be ' \
            'shared\n'.format(datetime.now()))  
    sb_methylated, sb_unmethylated = allocate_shared_buffer(_N.dtype, _N.shape), \
                                     allocate_shared_buffer(_N.dtype, _N.shape)
    arr_n, arr_m = get_np_view_of_buffer(sb_methylated, _N.dtype, _N.shape), \
                   get_np_view_of_buffer(sb_unmethylated, _N.dtype, _N.shape)
    arr_n.fill(np.nan)
    arr_m.fill(np.nan)

    f.write('[{}] done allocating shared memory, creating pools of ' \
            'workers and processing job list\n'.format(datetime.now()))
    mp.set_start_method('spawn')
    ps = mp.Pool(snakemake.threads, initializer=_init, \
                 initargs=(sb_methylated, sb_unmethylated, _N.dtype, _N.shape))

    try:
        results = ps.starmap(sequencing_error_correction, jobs)
    except (KeyboardInterrupt, SystemExit):
        ps.terminate()
        ps.join()
        sys.exit(1)

    ps.close()
    ps.join()    

    n_corrected, n_deleted = map(sum, zip(*results))

    f.write('[{}] gmelin-larch corrected {}/{} entries, among which {} were' \
            ' deleted\n'.format(datetime.now(), n_corrected, n_entries, n_deleted))
    np.savez(snakemake.output[0], n=arr_n, m=arr_m, cna=_C, rows=obj['rows'], cols=obj['cols'])                                                 

    f.write('[{}] DONE\n'.format(datetime.now()))
    f.close() 
