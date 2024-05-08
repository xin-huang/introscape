import numpy as np
import pandas as pd
import pybedtools
import tskit
import math
import pyranges as pr
from multiprocessing import Process, Manager

def create_final_tracts(tract_file):
    import pyranges as pr
    tract_file.columns=["Chromosome", "Start", "End", "Sample"]
    pr_tracts = pr.from_dict(tract_file)
    pr_tracts = pr_tracts.merge(strand=False, by='Sample')
    return pr_tracts


def partition_into_partitions(lst, num_partitions):
    partition_size = math.ceil(len(lst) / num_partitions)
    partitions = [lst[i:i+partition_size] for i in range(0, len(lst), partition_size)]
    return partitions

def tree_batch_processing(args):

    tree_batches, tree_object, migration_tree_object, bed_list, tgt_id = args

    ts = tree_object.value

    migration_ts = migration_tree_object.value

    for m in migration_ts.migrations():
        # Tree-sequences are sorted by the left ends of the intervals
        # Can skip those tree-sequences are not overlapped with the interval of i.

        for tree_index in tree_batches:
            t = ts.at_index(tree_index)

            if m.left >= t.interval.right: continue
            if m.right <= t.interval.left: break # [l, r)
            for n in ts.samples(tgt_id):
                if t.is_descendant(n, m.node):
                    left = m.left if m.left > t.interval.left else t.interval.left
                    right = m.right if m.right < t.interval.right else t.interval.right

                    bed_list.append([1, int(left), int(right), f'tsk_{ts.node(n).individual}_{int(n%ploidy+1)}'])

    return True
