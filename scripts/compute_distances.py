'''
    compute_distances.py <site_mask.npz> <status_likelihoods.npz>
'''

import sys, signal, numpy as np, multiprocessing as mp
from datetime import datetime
from buffer_util import get_np_view_of_buffer, allocate_shared_buffer
from itertools import combinations


def _init(buf_p00, buf_p10, buf_p11, buf_pwd, dtype_p, shape_p, dtype_pwd, shape_pwd):

    global sb_p00, sb_p10, sb_p11, sb_pwd

    sb_p00 = (buf_p00, dtype_p, shape_p)
    sb_p10 = (buf_p10, dtype_p, shape_p)
    sb_p11 = (buf_p11, dtype_p, shape_p)
    sb_pwd = (buf_pwd, dtype_pwd, shape_pwd)
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def enumerate_jobs(n_cells, d0010, d0011, d1011):

    entries = combinations(range(n_cells), 2)

    return [(i, j, d0010, d0011, d1011) for (i,j) in entries]


def get_total_distance(i,j, d0010, d0011, d1011):

    arr_pwd = get_np_view_of_buffer(*sb_pwd)
    arr_p00 = get_np_view_of_buffer(*sb_p00)
    arr_p10 = get_np_view_of_buffer(*sb_p10) 
    arr_p11 = get_np_view_of_buffer(*sb_p11)

    p00_i, p00_j = arr_p00[i,:], arr_p00[j,:]
    shared_sites = np.where(~np.isnan(p00_i) & ~np.isnan(p00_j))[0]
    
    p00_i, p00_j = p00_i[shared_sites], p00_j[shared_sites]
    p10_i, p10_j = arr_p10[i,:][shared_sites], arr_p10[j,:][shared_sites]
    p11_i, p11_j = arr_p11[i,:][shared_sites], arr_p11[j,:][shared_sites]

     # | p00_i                   | p10_i                   | p11_i                      |
     # |-------------------------|-------------------------|----------------------------|------
    ds = (0                      + (p10_i * p00_j * d0010) + (p11_i * p00_j * d0011) +  # p00_j
         (p00_i * p10_j * d0010) + 0                       + (p11_i * p10_j * d1011) +  # p10_j
         (p00_i * p11_j * d0011) + (p10_i * p11_j * d1011) + 0                       )  # p11_j

    total_dist = np.sum(ds) / shared_sites.shape[0]
    arr_pwd[i,j], arr_pwd[j,i] = total_dist, total_dist

    return


if __name__ == "__main__":

    f = open(snakemake.log[0], 'w')
    sys.stderr = sys.stdout = f

    f.write('[{}] gmelin-larch is iteratively constructing ' \
            'the methylation phylogeny\n'.format(datetime.now()))

    mask = (np.load(snakemake.input[0], allow_pickle=True)['mask']==np.inf)
    obj = np.load(snakemake.input[1], allow_pickle=True)
    _p00, _p10, _p11, cells, sites = obj['p00'][:,mask], obj['p10'][:,mask], \
                                     obj['p11'][:,mask], obj['rows'], obj['cols']

    d0010 = float(snakemake.params.d0010) 
    d0011 = float(snakemake.params.d0011)
    d1011 = float(snakemake.params.d1011)

    f.write('[{}] distance between homozygous unmethylated alleles and ' \
            'heterozygous alleles: {}\n'.format(datetime.now(), d0010))
    f.write('[{}] distance between homozygous methylated alleles and ' \
            'heterozygous alleles: {}\n'.format(datetime.now(), d1011))
    f.write('[{}] distance between homozygous unmethylated alleles and' \
            'homozygous methylated alleles: {}\n'.format(datetime.now(), d0011))
 
    f.write('[{}] computing pairwise expected distances\n'.format(datetime.now()))

    jobs = enumerate_jobs(cells.shape[0], d0010, d0011, d1011)

    sb_pwd = allocate_shared_buffer(np.float64, (cells.shape[0], cells.shape[0]))
    arr_pwd = get_np_view_of_buffer(sb_pwd, np.float64, (cells.shape[0], cells.shape[0]))
    arr_pwd.fill(np.nan)
    np.fill_diagonal(arr_pwd, 0)

    sb_p00, sb_p10, sb_p11 = allocate_shared_buffer(_p00.dtype, _p00.shape), \
                             allocate_shared_buffer(_p10.dtype, _p10.shape), \
                             allocate_shared_buffer(_p11.dtype, _p11.shape)
    arr_p00, arr_p10, arr_p11 = get_np_view_of_buffer(sb_p00, _p00.dtype, _p00.shape), \
                                get_np_view_of_buffer(sb_p10, _p10.dtype, _p10.shape), \
                                get_np_view_of_buffer(sb_p11, _p11.dtype, _p11.shape)
    arr_p00[:], arr_p10[:], arr_p11[:] = _p00[:], _p10[:], _p11[:]
    
    mp.set_start_method('spawn')
    ps = mp.Pool(snakemake.threads, initializer=_init, \
                 initargs=(sb_p00, sb_p10, sb_p11, sb_pwd, _p00.dtype, _p00.shape, \
                          np.float64, (cells.shape[0], cells.shape[0])))
    
    try:
        results = ps.starmap(get_total_distance, jobs)
    except (KeyboardInterrupt, SystemExit):
        ps.terminate()
        ps.join()
        sys.exit(1)

    ps.close()
    ps.join()

    np.savez(snakemake.output[0], pwd=arr_pwd, rows=cells)


    f.write('[{}] DONE\n'.format(datetime.now()))
    f.close()
