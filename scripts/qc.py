'''
    qc.py <methylated_rc.npz> <unmethylated_rc.npz> <[optional]cna>
'''

import sys, numpy as np
from datetime import datetime

f = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout = f

if __name__ == "__main__":

    f.write('[{}] gmelin-larch is filtering your input' \
            ' matrices\n'.format(datetime.now()))

    methylated_reads = np.load(snakemake.input[0], allow_pickle=True)
    unmethylated_reads = np.load(snakemake.input[1], allow_pickle=True)
    assert np.array_equal(methylated_reads['rows'], unmethylated_reads['rows']), \
           'Read count matrices cover different sets of cells.'
    assert np.array_equal(methylated_reads['cols'], unmethylated_reads['cols']), \
           'Read count matrices cover different sets of sites.'
    assert np.array_equal(np.isnan(methylated_reads['m']), \
                          np.isnan(unmethylated_reads['m'])), \
           'Read count matrices disagree on cells or sites covered.'

    if len(snakemake.input) == 3: # optional CNA input
        _cna = np.load(snakemake.input[2], allow_pickle=True) 
        assert np.array_equal(methylated_reads['rows'], _cna['rows']), \
               'CNA matrix covers a different set of cells.'
        assert np.array_equal(methylated_reads['rows'], _cna['rows']), \
               'CNA matrix covers a different set of sites.'
        cna = _cna['m']
        cna[np.isnan(cna)] = snakemake.params.default_cna
        cna = cna.astype(int)
    else:
        cna = np.full(methylated_reads['m'].shape, snakemake.params.default_cna, \
                      dtype=int)

    # TODO: Quality control for sites and cells

    np.savez(snakemake.output[0], m=methylated_reads['m'], \
             n=unmethylated_reads['m'], rows=methylated_reads['rows'], \
             cols=methylated_reads['cols'], cna=cna)

f.write('[{}] DONE\n'.format(datetime.now()))
f.close()
