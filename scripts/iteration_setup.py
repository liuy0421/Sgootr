'''
    iteration_setup.py <input.npz>
'''

import sys, signal, numpy as np, multiprocessing as mp
from datetime import datetime
from buffer_util import get_np_view_of_buffer, allocate_shared_buffer


def _init(buf_s00, buf_s10, buf_s11, dtype, shape):
    
    global sb_s00, sb_s10, sb_s11

    sb_s00 = (buf_s00, dtype, shape)
    sb_s10 = (buf_s10, dtype, shape)
    sb_s11 = (buf_s11, dtype, shape)
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def enumerate_jobs(N, M, C, p00, p10, p11, p):

    entries = np.where(~np.isnan(N))

    return [(i, j, N[i,j], M[i,j], C[i,j], p00, p10, p11, p) for \
            (i,j) in zip(entries[0], entries[1])]


def get_hom_het_likelihood(n, c, p_hom, p_het, p):
    '''
        compute likelihood for homozygous status and heterozygous status
        given n reads of the same status, copy number, prior probabilities,
        and probability of drawing the allele of the methylation status given
        a pair of alleles with different methylation statuses
    '''   
    preads_het = sum([( (gamma*p) / ((gamma*p)+(c-gamma)*(1-p)) )**n \
                      for gamma in range(1,c)]) / (c-1)
    a = p_hom + p_het * preads_het

    return (p_hom / a), (p_het * preads_het / a)
    

def compute_status_likelihood(i, j, n, m, c, p00, p10, p11, p):
    '''
        i: cell index
        j: site index
        n: number of methylated reads
        m: number of unmethylated reads
        c: copy number call at site j in cell i
        p00: prior probability of a site being homozygous unmethylated
        p10: prior probability of a site being heterozygous methylated
        p11: prior probability of a site being homozygous methylated 
    '''
    arr_s00 = get_np_view_of_buffer(*sb_s00)
    arr_s10 = get_np_view_of_buffer(*sb_s10)
    arr_s11 = get_np_view_of_buffer(*sb_s11)

    assert c > 0, 'observing reads at CNA=0 site'

    ####
    #   1. in case of chromosomal loss, only the homozygous methylation
    #      status is possible 
    ####
    if c == 1: 
        arr_s10[i,j] = 0 
        if n == 0:
            assert m > 0, 'invalid number of unmethylated reads'
            arr_s00[i,j], arr_s11[i,j] = 1, 0
        elif m == 0:
            arr_s00[i,j], arr_s11[i,j] = 0, 1
        else:
            raise ValueError('having reads of both statuses at a CNA=1 site')
        return

    ####
    #   2. if c > 1 and there are reads of both methylation statuses,
    #      they can only come from heterozygous alleles
    ####
    if n > 0 and m > 0:
        arr_s00[i,j], arr_s10[i,j], arr_s11[i,j] = 0, 1, 0
        return
    
    #### 
    #   3. if c > 1 and there are only reads from one methylation status  
    ####
    if n > 0 and m == 0: # if there are only methylated reads
        arr_s00[i,j] = 0
        arr_s11[i,j], arr_s10[i,j] = get_hom_het_likelihood(n, c, p11, p10, p)
        return
    elif n ==0 and m > 0: # if there are only unmethylated reads
        arr_s11[i,j] = 0
        arr_s00[i,j], arr_s10[i,j] = get_hom_het_likelihood(m, c, p00, p10, 1-p)
        return

    raise ValueError('at index [{},{}], we observed n={}, m={}, c={} while ' \
                     'that combination should have been impossible'.format(i,j,n,m,c))



if __name__ == "__main__":

    f = open(snakemake.log[0], 'w')
    sys.stderr = sys.stdout = f

    f.write('[{}] gmelin-larch is setting up the iterative ' \
            'procedure\n'.format(datetime.now()))

    p00, p10, p11, p = float(snakemake.params.p00), float(snakemake.params.p10), \
                       float(snakemake.params.p11), float(snakemake.params.p)
    denom = p00+p10+p11
    if denom != 1:
        f.write('[{}] P(00), P(10), and P(11) do not sum to 1. Performing ' \
                'normalization\n'.format(datetime.now(), p00, p10, p11))
        p00, p10, p11 = p00/denom, p10/denom, p11/denom
    assert (p>=0 and p<=1), 'p is not a probability'

    sct = int(snakemake.params.status_confidence_threshold)
    f.write('[{}] proceeding with parameters P(00)={}, P(10)={}, P(11)={}, p={}, ' \
            'and status_confidence_threshold={}\n'.format(datetime.now(), p00, \
                                                          p10, p11, p, sct))

    obj = np.load(snakemake.input[0], allow_pickle=True)
    _N, _M, _C, cells, sites = obj['n'], obj['m'], obj['cna'], obj['rows'], \
                               obj['cols']

    ####
    #   1. create site mask
    ####
    f.write('[{}] creating site mask\n'.format(datetime.now()))
    mask = np.full(shape=sites.shape, fill_value=np.inf)
    np.savez(snakemake.output[0], mask=mask)


    ####
    #   2. heuristically calling methylation statuses for persistence computation
    ####
    f.write('[{}] done creating mask, heuristically calling site methylation ' \
            'statuses in cells\n'.format(datetime.now()))
    _D = _N + _M
    _R = _N / _D
    _R[(_R > 0) & (_R < 1)] = .5
    _R[_D < sct] = np.nan # need at least sct reads to call status
    np.savez(snakemake.output[1], s=_R)


    ####
    #   3. compute site status likelihood
    ####
    f.write('[{}] done heuristically calling methylation statuses, computing ' \
            ' site status likelihood in cells\n'.format(datetime.now()))
    jobs = enumerate_jobs(_N, _M, _C, p00, p10, p11, p)
    s00, s10, s11 = allocate_shared_buffer(np.float64, _N.shape), \
                    allocate_shared_buffer(np.float64, _N.shape), \
                    allocate_shared_buffer(np.float64, _N.shape)
    arr_s00, arr_s10, arr_s11 = get_np_view_of_buffer(s00, np.float64, _N.shape), \
                                get_np_view_of_buffer(s10, np.float64, _N.shape), \
                                get_np_view_of_buffer(s11, np.float64, _N.shape)
    arr_s00.fill(np.nan)
    arr_s10.fill(np.nan)
    arr_s11.fill(np.nan)

    mp.set_start_method('spawn')
    ps = mp.Pool(snakemake.threads, initializer=_init, \
                 initargs=(s00, s10, s11, np.float64, _N.shape))    

    try:
        results = ps.starmap(compute_status_likelihood, jobs)
    except (KeyboardInterrupt, SystemExit):
        ps.terminate()
        ps.join()
        sys.exit(1)

    ps.close()
    ps.join()

    f.write('[{}] done computing status likelihood, writing to ' \
            'output\n'.format(datetime.now()))
    np.savez(snakemake.output[2], p00=arr_s00, p10=arr_s10, p11=arr_s11, \
             rows=cells, cols=sites)


    f.write('[{}] DONE\n'.format(datetime.now()))
    f.close()
