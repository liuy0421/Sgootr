'''
    prune.py <t{i}.nwk> <heuristically_called_statuses.npz> <t{i}_site_mask.npz>
'''

import sys, signal, numpy as np, multiprocessing as mp
from datetime import datetime
import skbio
from buffer_util import get_np_view_of_buffer, allocate_shared_buffer
from scipy.spatial.distance import jensenshannon


def _init(buf_mps, dtype_mps, shape_mps, buf_sts, dtype_sts, shape_sts):

    global sb_mps, sb_sts

    sb_mps = (buf_mps, dtype_mps, shape_mps)
    sb_sts = (buf_sts, dtype_sts, shape_sts)
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def enumerate_jobs(internal_leaves, partition_validity_threshold):

    entries = range(len(internal_leaves)) 

    return [(i, internal_leaves[i], partition_validity_threshold) \
            for i in entries]


def compute_internal_persistence_score(i, leaves, partition_validity_threshold):

    arr_mps = get_np_view_of_buffer(*sb_mps)
    arr_sts = get_np_view_of_buffer(*sb_sts)    

    subtree1_indices = [int(leaf) for leaf in leaves]
    subtree1_statuses = arr_sts[subtree1_indices,:]
    subtree2_statuses = np.delete(arr_sts, subtree1_indices, axis=0)

    '''
        if either subset of the partition has fewer than 50% of its leaves with 
        heuristically called methylation status (due to lack of read coverage),
        we cannot generate valid persistence score for that site at the partition
    '''
    no_status_1 = np.sum(np.isnan(subtree1_statuses), axis=0)
    no_status_2 = np.sum(np.isnan(subtree2_statuses), axis=0)
    size1, size2 = subtree1_statuses.shape[0], subtree2_statuses.shape[0]
    valid1, valid2 = no_status_1 < partition_validity_threshold * size1, \
                     no_status_2 < partition_validity_threshold * size2 

    '''
        construct leaf status distributions for subsets induced by given partition
    '''
    _p = np.array([np.count_nonzero(subtree1_statuses==s, axis=0) for s in [0,.5,1]])
    _q = np.array([np.count_nonzero(subtree2_statuses==s, axis=0) for s in [0,.5,1]]) 
    _p[:,~valid1], _q[:,~valid2] = 1, 1 # suppress runtime warning

    arr_mps[i,:] = jensenshannon(_p, _q, axis=0)
    arr_mps[i,(~valid1 | ~valid2)] = np.nan



if __name__ == "__main__":

    # work-around for snakemake env bug
    f_nwk, f_msk, f_sts, f_out, f_sco, f_kpa, f_pvt, f_mss, f_tds, f_log = sys.argv[1:]

    f = open(f_log, 'w')
    sys.stderr = sys.stdout = f

    kappa = float(f_kpa)
    assert (kappa > 0) and (kappa < 1), 'Kappa must be a fraction.'
    partition_validity_threshold, minimum_subtree_size = float(f_pvt), float(f_mss)
    assert (partition_validity_threshold > 0) and (partition_validity_threshold < 1), \
           'partition_validity_threshold must be a fraction.'    
    assert (minimum_subtree_size > 0) and (minimum_subtree_size < 1), \
           'minimum_subtree_size must be a fraction.'


    f.write('[{}] gmelin-larch is perfroming site pruning with ' \
            'kappa={}\n'.format(datetime.now(), kappa))

    tree = skbio.TreeNode.read(f_nwk)
    msk = np.load(f_msk, allow_pickle=True)['mask']
    mask = np.isinf(msk)
    sts = np.load(f_sts, allow_pickle=True)['s']

    ####
    #   1. compute persistence score for each remaining site at each internal node
    ####     
    min_nodes = round(sts.shape[0] * minimum_subtree_size)
    max_nodes = sts.shape[0] - min_nodes
    f.write('[{}] gmelin-larch is computing persistence scores for {} ' \
            'remaining sites\n'.format(datetime.now(), np.sum(mask)))
    internal_leaves = [n.subset() for n in tree.non_tips() if \
                       len(n.subset()) >= min_nodes and \
                       len(n.subset()) <= max_nodes]
    jobs = enumerate_jobs(internal_leaves, partition_validity_threshold)
    mps = allocate_shared_buffer(np.float64, (len(internal_leaves), mask.shape[0]))
    statuses = allocate_shared_buffer(sts.dtype, sts.shape)
    arr_mps = get_np_view_of_buffer(mps, np.float64, (len(internal_leaves), mask.shape[0]))
    arr_statuses = get_np_view_of_buffer(statuses, sts.dtype, sts.shape)
    arr_mps.fill(np.nan)
    arr_statuses[:] = sts[:]
    arr_statuses[:,~mask] = np.nan
    
    mp.set_start_method('spawn')
    ps = mp.Pool(int(f_tds), initializer=_init, \
                 initargs=(mps, np.float64, (len(internal_leaves), mask.shape[0]), \
                           statuses, sts.dtype, sts.shape))

    try:
        results = ps.starmap(compute_internal_persistence_score, jobs)
    except (KeyboardInterrupt, SystemExit):
        ps.terminate()
        ps.join()
        sys.exit(1)

    ps.close()
    ps.join()

    ####
    #   2. for each site, take the max score across all internal nodes as its
    #      final tree persistence score  
    ####
    persistence_scores = np.nanmax(arr_mps, axis=0)
    no_scores = np.sum(np.isnan(persistence_scores[mask]))
    persistence_scores[~mask] = np.nan

    f.write('[{}] done computing persistence scores ({}/{} sites have no scores)' \
            ', pruning sites\n'.format(datetime.now(), no_scores, np.sum(mask)))
  
    ####
    #   3. update mask by excluding a) sites without valid persistence scores, and b)
    #      kappa of the remaining sites with the lowest persistence scores (number of
    #      sites pruned may sligtly exceed kappa in case of ties)  
    ####
    scores = persistence_scores.copy()
    scores.sort()
    cutoff = scores[int(np.sum(~np.isnan(scores))*kappa)]
    pruned = (persistence_scores <= cutoff) | \
             (np.isnan(persistence_scores) & mask)

    f.write('[{}] with a {} cutoff of {}, {} sites have persistence ' \
            'score <= to that\n'.format(datetime.now(), kappa, cutoff, \
                                      np.sum(persistence_scores <= cutoff)))

    if np.sum(mask) == mask.shape[0]: # if at iteration 0
        iteration = 0
    else:
        iteration = np.amax(msk[~np.isinf(msk)]) + 1

    msk[pruned] = iteration

    np.savez(f_out, mask=msk)
    np.savez(f_sco, scores=persistence_scores)
    f.write('[{}] {} sites remaining after pruning\n'.format(datetime.now(), \
                                                             np.sum(np.isinf(msk))))

    f.write('[{}] DONE\n'.format(datetime.now()))
    f.close()
