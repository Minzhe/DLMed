###################################################################################
###                            merge.crispr.gi.py                               ###
###################################################################################
# curate weissman's cell paper of dual crispr knockout

import numpy as np
import pandas as pd
import os
import utility as util


proj_path = 'D:/projects/DLCell'
nbt_path = os.path.join(proj_path, 'data/curated/crispr.doubleKO.nbt3834.geno.pheno.norm.dense.csv')
cell_path = os.path.join(proj_path, 'data/curated/crispr.doubleKO.cell.geno.pheno.norm.dense.csv')
out_dir = os.path.join(proj_path, 'data/curated')

##############################    function    #############################
def mergeScoreTable(neg2keep, **scores):
    '''
    Merge phenotype score table from different datasets.
    '''
    names = scores.keys()
    scoretable = []
    for name in names:
        scoretable.append(scores[name] if neg2keep == name else util.rmNegative(scores[name]))
    return pd.concat(scoretable, axis=0, ignore_index=True)


##############################    main    #############################
crispr_nbt = pd.read_csv(nbt_path)
crispr_cell = pd.read_csv(cell_path)
crispr_merge = mergeScoreTable(neg2keep='cell', nbt=crispr_nbt, cell=crispr_cell)
# util.plotDist(merge=crispr_merge)
crispr_merge_sparse, duptable = util.createTrainingTable(crispr_merge, duplicate='average')
crispr_merge_sparse.to_csv(os.path.join(out_dir, 'crispr.doubleKO.merged.geno.pheno.norm.sparse.csv'), index=None)
duptable.to_csv(os.path.join(out_dir, 'crispr.doubleKO.merged.duplicate.csv'), index=None)