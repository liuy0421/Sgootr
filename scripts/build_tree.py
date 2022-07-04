'''
    build_tree.py <t{i}_pairwise_distances.npz> <t{i-1}.nwk> <t{i-1}_RF.txt>
'''

import numpy as np, sys
from datetime import datetime
import skbio


def name_internal_nodes(t):
    
    i = 0

    for node in t.traverse():
        if node.name == None:
            node.name = 'internal {}'.format(i)
            i+=1
    
    return t


def build_and_reroot(dm, root):

    unrooted_tree = skbio.tree.nj(dm)
    root_parent_node = unrooted_tree.find(root).parent
    tree = unrooted_tree.root_at(root_parent_node)
    root_node = [x for x in tree.children if x.name == root][0]
    tree.remove(root_node)

    return name_internal_nodes(tree)


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
    at_t0 = (len(sys.argv) == 6)
    if at_t0:
        f_pwd, f_out_nwk, f_out_rf, root, f_log = sys.argv[1:]
    else:
        f_pwd, f_old_nwk, f_old_rf, f_out_nwk, f_out_rf, root, f_log = sys.argv[1:] 

    f = open(f_log, 'w')
    sys.stderr = sys.stdout = f

    obj = np.load(f_pwd, allow_pickle=True)
    pwd, cells = obj['pwd'], obj['rows']
    assert root in cells, 'Selected root not in the input set of cells.'

    ####
    #   1. construct tree with neighbor-joining (Saitou & Nei, 1987)
    ####
    f.write('[{}] gmelin-larch is building neighbor-joining '\
            'tree with {} as root\n'.format(datetime.now(), root))

    dm = skbio.DistanceMatrix(pwd, [str(x) for x in range(len(cells))])
    tree = build_and_reroot(dm, str(np.where(cells==root)[0][0]))

    with open(f_out_nwk, 'w') as fi:
        fi.write(str(tree))

    ####
    #   2. compute RF distance if t{-1} exists
    ####
    if at_t0:
        with open(f_out_rf, 'w') as fi:
            fi.write('')
    else:

        rf_dist = compute_rf_distance(tree, skbio.TreeNode.read(f_old_nwk))
   
        with open(f_old_rf, 'r') as fi:
            old_rf = fi.readlines()
        with open(f_out_rf, 'w') as fi:
            fi.write(''.join(old_rf))
            fi.write('{}\n'.format(rf_dist))
            

    f.write('[{}] DONE\n'.format(datetime.now()))
    f.close()
