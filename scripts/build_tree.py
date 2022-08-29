'''
    build_tree.py <t{i}_pairwise_distances.npz> <t{i-1}.nwk> <t{i-1}_RF.txt>
'''

import numpy as np, sys, os
from datetime import datetime
import skbio


def name_internal_nodes(t):
    
    i = 0

    for node in t.traverse():
        if node.name == None:
            node.name = 'internal {}'.format(i)
            i+=1
    
    return t


def reroot(unrooted_tree, root):

    root_parent_node = unrooted_tree.find(root).parent
    tree = unrooted_tree.root_at(root_parent_node) 
    root_node = [x for x in tree.children if x.name == root][0]
    tree.remove(root_node)
   
    return name_internal_nodes(tree)


def build(dm, root, algo):

    if algo == 'N':
        unrooted_tree = skbio.tree.nj(dm)
    elif algo == 'F':
        f_out = dm[:-10] + 'fastme.nwk'
        os.system('fastme -i {} -o {} -n -s'.format(dm, f_out))
        while not os.path.exists(f_out):
            continue
        unrooted_tree = skbio.TreeNode.read(f_out)

    return reroot(unrooted_tree, root)


def compute_rf_distance(t1, t2):
     
    leaves1, leaves2 = set([x.name for x in t1.tips()]), \
                      set([x.name for x in t2.tips()])
    assert leaves1 == leaves2, 'when computing RF distance, trees contain ' \
                               'different set of leaf nodes.'
    
    subtrees1 = set(internal.subset() for internal in t1.non_tips())
    subtrees2 = set(internal.subset() for internal in t2.non_tips())

    return len(subtrees1.symmetric_difference(subtrees2))


if __name__ == "__main__":
    
    # work-around for snakemake env bug
    at_t0 = (len(sys.argv) == 7)
    if at_t0:
        f_out_nwk, f_out_rf, root, algo, f_log, f_pwd = sys.argv[1:]
    else:
        f_out_nwk, f_out_rf, root, algo, f_log, f_pwd, f_old_rf = sys.argv[1:8]
        f_old_nwks = sys.argv[8:] 

    assert (algo == 'N') or (algo == 'F'), 'Specify tree construction algorithm' \
                                           'N(eighbor-joining) or F(astME)'

    f = open(f_log, 'w')
    sys.stderr = sys.stdout = f

    obj = np.load(f_pwd, allow_pickle=True)
    pwd, cells = obj['pwd'], obj['rows']
    assert root in cells, 'Selected root not in the input set of cells.'

    ####
    #   1. construct tree with neighbor-joining (Saitou & Nei, 1987) or
    #      FastME 2.0 (Lefort et al., 2015)
    ####
    if algo == 'N': # If algorithm is Neighbor-joining
        f.write('[{}] gmelin-larch is building neighbor-joining '\
                'tree with {} as root\n'.format(datetime.now(), root))

        dm = skbio.DistanceMatrix(pwd, [str(x) for x in range(len(cells))])
        tree = build(dm, str(np.where(cells==root)[0][0]), algo)

    elif algo == 'F': # If algorithm is FastME 2.0
        f.write('[{}] gmelin-larch is building FastME 2.0 ' \
                'tree with {} as root\n'.format(datetime.now(), root))    
        f_mat = f_out_nwk[:-8] + 'fastme.mat'
        with open(f_mat, 'w') as fi:
            n_cells = pwd.shape[0]
            fi.write('{}\n'.format(n_cells))
            for i in range(n_cells):
                fi.write('{}\t{}\n'.format(i, \
                                           '\t'.join([str(x) for x in pwd[i]])))
        tree = build(f_mat, str(np.where(cells==root)[0][0]), algo)

    with open(f_out_nwk, 'w') as fi:
        fi.write(str(tree))


    ####
    #   2. compute RF distance if t{-1} exists
    ####
    if at_t0:
        with open(f_out_rf, 'w') as fi:
            fi.write('')
    else:

        rf_dists = [str(compute_rf_distance(tree, skbio.TreeNode.read(old_tree))) \
                    for old_tree in f_old_nwks]
   
        with open(f_old_rf, 'r') as fi:
            old_rf = fi.readlines()
        with open(f_out_rf, 'w') as fi:
            fi.write(''.join(old_rf))
            fi.write('{}\n'.format('\t'.join(rf_dists)))
            

    f.write('[{}] DONE\n'.format(datetime.now()))
    f.close()
