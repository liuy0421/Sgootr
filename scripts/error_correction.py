'''
    error_correction.py <filtered_rc_cna.npz>    
'''
import sys, numpy as np, multiprocessing as mp

f = open(snakemake.log[0],'w')
sys.stderr = sys.stdout = f


def enumerate_jobs(N, M, C, e):
   
    entries = np.where(~np.isnan(N))

    return [(i, j, N[i,j], M[i,j], C[i,j], e) for \
            (i,j) in zip(entries[0], entries[1])]


def safe_log(x): 

    return np.log(x*(10**100)+10**(-300))    


'''
    sequencing_error_correction(i,j,n,m,c,e)

    i: cell index
    j: site index
    n: number of methylated reads
    m: number of unmethylated reads
    c: copy number call at site j in cell i
    e: sequencing error parameter
'''
def seqencing_error_correction(i,j,n,m,c,e):

    if c == 0:
        return (i,j),(np.nan,np.nan), 1, 1

    L = [0 for gamma in range(c+1)]
    L[0] = safe_log((e**n)*((1-e)**m)) 
    L[c] = safe_log(((1-e)**n)*(e**m))

    for gamma in range(1,c):
        L[gamma] = safe_log(((gamma/c)**n)*((1-gamma/c)**m))

    status = L.index(max(L))

    if (status == 0 or status == c) and (n != 0 and m != 0):
        if n > m:
            return (i,j), (n+m, 0), 1, 0
        elif n < m:
            return (i,j), (0, n+m), 1, 0
        else:
            return (i,j), (np.nan, np.nan), 1, 1

    return (i,j), (n, m), 0, 0


if __name__ == "__main__":
    
    f.write('..gmelin-larch is performing read error correction') 

    obj = np.load(snakemake.input[0], allow_pickle=True)
    _N, _M, _C = obj['n'], obj['m'], obj['c']
    e = snakemake.params.e
    n_entries = np.sum(~np.isnan(_N))

    ps, jobs = mp.Pool(snakemake.threads), enumerate_jobs(_N, _M, _C, e)
    results = ps.starmap(sequencing_error_correction, jobs)
    ps.close()

    N, M = np.full(_N.shape, np.nan), np.full(_M.shape, np.nan)
    n_corrected, n_deleted = 0, 0

    for k,v,i,d in results:
        obj['n'][k], obj['m'][k] = v 
        n_corrected += i
        n_deleted += d

    print('..gmelin-larch corrected {}/{} entries, among which {} were' \
          ' deleted'.format(n_corrected, n_entries, n_deleted))
    np.savez(snakemake.output[0], n=N, m=M, rows=obj['rows'], cols=obj['cols'])                                                 

f.write('[DONE]')
f.close() 
