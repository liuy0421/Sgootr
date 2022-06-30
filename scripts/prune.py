'''
    prune.py <input.npz>
'''

import sys, numpy as np
from datetime import datetime


if __name__ == "__main__":

    f = open(snakemake.log[0], 'w')
    sys.stderr = sys.stdout = f

    f.write('[{}] gmelin-larch is iteratively constructing ' \
            'the methylation phylogeny\n'.format(datetime.now()))

    obj = np.load(snakemake.input[0], allow_pickle=True)
    _N, _M, _C, cells, sites = obj['n'], obj['m'], obj['cna'], obj['rows'], obj['cols']
    max_iter = int(snakemake.params.max_iter)
    kappa = float(snakemake.params.kappa)
    root = str(snakemake.params.root)

    assert root in cells, 'Selected root not in the input set of cells.'

    '''
        1. for each non-NA site in a cell, compute 3-tuple corresponding to 
           probability of it being homozygously methylated, heterozygous, or
           homozygously unmethylated
    '''
        #TODO
        
    ''' 
        2. allocate memory for 1) cell-by-cell pairwise distance matrix, and
           2) site persistence score array
    '''
        #TODO

    '''
        3. while round <= max_iter, 
    '''
        #TODO
        # 3a. compute pair-wise distances
        # 3b. run NJ to get (rooted) nwk
        # 3c. generate visualization w/R
        # 3d. if round > 0, compute RF distance from last tree
        # 3e. if round < max_iter, for each site, compute persistence score, 
        #     and create site mask for selected/pruned sites for the following
        #     iteration
    

    f.write('[DONE]')
    f.close()
