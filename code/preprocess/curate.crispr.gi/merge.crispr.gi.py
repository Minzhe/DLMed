###################################################################################
###                            merge.crispr.gi.py                               ###
###################################################################################
# curate weissman's cell paper of dual crispr knockout

proj_path = '/work/bioinformatics/s418336/projects/DLMed'
import numpy as np
import pandas as pd
import os
import sys
sys.path.append(os.path.join(proj_path, 'code'))
import utility.plot as p
import matplotlib.pyplot as plt



##############################    function    #############################
# def mergeScoreTable(neg2keep, **scores):
#     '''
#     Merge phenotype score table from different datasets.
#     '''
#     names = scores.keys()
#     scoretable = []
#     for name in names:
#         scoretable.append(scores[name] if neg2keep == name else util.rmNegative(scores[name]))
#     return pd.concat(scoretable, axis=0, ignore_index=True)

def create_dummies(series1, series2):
    df1 = create_dummy(series1)
    df2 = create_dummy(series2)
    add1 = list(set(df1) - set(df2)) # AARS2
    add2 = list(set(df2) - set(df1)) # None
    for x in add1:
        df2[x] = 0
    for x in add2:
        df1[x] = 0
    df1.sort_index(axis=1, inplace=True)
    df2 = df2.loc[:,df1.columns]

    return df1 + df2

def create_dummy(series):
    df = pd.get_dummies(series)
    add = list(set(series) - set(df.columns))
    assert len(add) <= 1
    if len(add) == 1:
        df[add[0]] = 0
        df[series == add[0]] = 1
    df.sort_index(axis=1, inplace=True)

    return df

##############################    phenotype score    #############################
# crispr_nbt = pd.read_csv(nbt_path)
# crispr_cell = pd.read_csv(cell_path)
# crispr_merge = mergeScoreTable(neg2keep='cell', nbt=crispr_nbt, cell=crispr_cell)
# # util.plotDist(merge=crispr_merge)
# crispr_merge_sparse, duptable = util.createTrainingTable(crispr_merge, duplicate='average')
# crispr_merge_sparse.to_csv(os.path.join(out_dir, 'crispr.doubleKO.merged.geno.pheno.norm.sparse.csv'), index=None)
# duptable.to_csv(os.path.join(out_dir, 'crispr.doubleKO.merged.duplicate.csv'), index=None)


##############################    gi score    #############################
nbt_path = os.path.join(proj_path, 'data/curated/K562/crispr.doubleKO.nbt3834.geno.giscore.dense.csv')
cell_path = os.path.join(proj_path, 'data/curated/K562/crispr.doubleKO.cell.geno.giscore.dense.csv')
dense_out_path = os.path.join(proj_path, 'data/curated/K562/crispr.doubleKO.merged.geno.giscore.dense.csv')
sparse_out_path = os.path.join(proj_path, 'data/curated/K562/crispr.doubleKO.merged.geno.giscore.sparse.csv')
gi_nbt = pd.read_csv(nbt_path)
gi_cell = pd.read_csv(cell_path)

# f, ax = plt.subplots()
# p.plot_density(title='GI score density', ax=ax, cell_gi=gi_cell['gi_score'], nbt_gi_t=gi_nbt['gi_t_score'])
# f.savefig(os.path.join(proj_path, 'data/CRISPR.GI/density.png'))

# merge
gi_nbt = gi_nbt[['Gene1', 'Gene2', 'gi_t_score']]
gi_nbt.columns = ['Gene1', 'Gene2', 'gi_score']
gi_table = pd.concat([gi_cell, gi_nbt])
gi_table = gi_table.loc[~gi_table.duplicated(subset=['Gene1', 'Gene2'], keep=False),:].sort_values(by=['Gene1', 'Gene2'])
gi_table.index = range(gi_table.shape[0])
gi_table.to_csv(dense_out_path, index=None)

# f, ax = plt.subplots()
# p.plot_density(title='GI score density', ax=ax, gi=gi_table['gi_score'])
# f.savefig(os.path.join(proj_path, 'data/CRISPR.GI/density_merged.png'))

# sparse table
gi_mat = create_dummies(gi_table['Gene1'], gi_table['Gene2'])
gi_mat['gi_score'] = gi_table['gi_score']
gi_mat.to_csv(sparse_out_path, index=None)
