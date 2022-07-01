'''
    build_tree.py <pwd.npz>
'''

import numpy as np, sys
from datetime import datetime
import skbio


def build_and_reroot(dm, root):
    '''
        construct tree wiht neighbor-joining (Daitou & Nei, 1987)
    '''
    unrooted_tree = skbio.tree.nj(dm)
    root_parent_node = unrooted_tree.find(root).parent
    tree = unrooted_tree.root_at(root_parent_node)
    root_node = [x for x in tree.children if x.name == root][0]
    tree.remove(root_node)

    return tree


if __name__ == "__main__":
    
    # work-around for snakemake env bug
    f_in, f_out, root, f_log = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    f = open(f_log, 'w')
    sys.stderr = sys.stdout = f

    obj = np.load(f_in, allow_pickle=True)
    pwd, cells = obj['pwd'], obj['rows']
    assert root in cells, 'Selected root not in the input set of cells.'
    root = ' '.join(root.split('_'))

    f.write('[{}] gmelin-larch is building neighbor-joining tree\n'.format(datetime.now()))

    dm = skbio.DistanceMatrix(pwd, cells)
    tree = build_and_reroot(dm, root)

    with open(f_out, 'w') as fi:
        fi.write(str(tree))

    f.write('[{}] DONE\n'.format(datetime.now()))
    f.close()
