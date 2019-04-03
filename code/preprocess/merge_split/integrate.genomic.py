#####################################################################################
###                            integrate.genomic.py                               ###
#####################################################################################

import os
import re
import sys
import numpy as np
from sklearn.preprocessing import scale
import pandas as pd
import pickle as pkl
proj_dir = '/work/bioinformatics/s418336/projects/DLMed/'
data_dir = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version')
sys.path.append(os.path.join(proj_dir, 'code'))
from utility import integrate as itg

genelist = 'cancergene'

###############################    cancer genes    ##################################
if genelist == 'cancergene':
    mut_path = os.path.join(data_dir, 'ccle_utsw.lung_Mutation_cancergene.csv')
    expr_path = os.path.join(data_dir, 'ccle_utsw.lung_RNAseq_cancergene.csv')
    drug_path = os.path.join(data_dir, 'ccle_utsw.lung_drug.split.csv')
    out_point_path = os.path.join(data_dir, 'ccle_utsw.lung_MutGeneExpr_cancergene_drugRes.gene_level_point.pkl')
    out_cell_path = os.path.join(data_dir, 'ccle_utsw.lung_MutGeneExpr_cancergene_drugRes.gene_level_cell.pkl')

    mut = pd.read_csv(mut_path, usecols=['Cell', 'Gene']).drop_duplicates()
    expr = pd.read_csv(expr_path, index_col=0)
    expr.columns = ['expr_'+gene for gene in expr.columns.values]
    # expr = pd.DataFrame(scale(expr), index=expr.index, columns=expr.columns)
    drug = pd.read_csv(drug_path)

    cell_id = list(drug.Cell)
    mut = mut.loc[mut.Cell.isin(cell_id)]
    expr = expr.loc[expr.index.to_series().isin(cell_id),:]

    # make mutation array
    mut['value'] = 1
    mut = mut.pivot(index='Cell', columns='Gene', values='value').fillna(0)

    # impute expr
    mean_expr = expr.mean(axis=0)
    expr = pd.concat([expr, pd.DataFrame(index=list(set(cell_id)-set(expr.index)))], sort=False).sort_index()
    expr = expr.apply(lambda x: mean_expr if all(x.isnull()) else x, axis=1)

    # concate mut, expr and cnv
    y = drug.LOG50
    drug_dummy = pd.get_dummies(drug.Drug)
    cell_mat = pd.concat([mut, expr], axis=1, sort=False)
    cell_mat = drug.Cell.apply(lambda x: cell_mat.loc[x,:])
    X = pd.concat([cell_mat, drug_dummy], axis=1, sort=False)

    # split train, test, val on point
    train_idx = drug.Train_split_point == 'Train'
    val_idx = drug.Train_split_point == 'Val'
    test_idx = drug.Train_split_point == 'Test'
    infer_idx = (drug.Train_split_point == 'Inference') & (~drug.Dropped_previously)
    infer_full_idx = drug.Train_split_point == 'Inference'

    X_train, y_train = X.loc[train_idx,:], y[train_idx]
    X_val, y_val = X.loc[val_idx,:], y[val_idx]
    X_test, y_test = X.loc[test_idx,:], y[test_idx]
    X_infer, y_infer = X.loc[infer_idx,:], y[infer_idx]
    X_infer_full, y_infer_full = X.loc[infer_full_idx,:], y[infer_full_idx]

    with open(out_point_path, 'wb') as f:
        pkl.dump({'train': (X_train, y_train), 
                'val': (X_val, y_val),
                'test': (X_test, y_test),
                'inference': (X_infer, y_infer),
                'inference_full': (X_infer_full, y_infer_full)}, file=f)
    
    # split train, test, val on cell
    train_idx = drug.Train_split_cell == 'Train'
    val_idx = drug.Train_split_cell == 'Val'
    infer_idx = drug.Train_split_point == 'Inference'

    X_train, y_train = X.loc[train_idx,:], y[train_idx]
    X_val, y_val = X.loc[val_idx,:], y[val_idx]
    X_infer, y_infer = X.loc[infer_idx,:], y[infer_idx]

    with open(out_cell_path, 'wb') as f:
        pkl.dump({'train': (X_train, y_train), 
                'val': (X_val, y_val),
                'inference': (X_infer, y_infer)}, file=f)

###############################    cancer genes    ##################################
elif genelist == 'allgene':
    mut_path = os.path.join(data_dir, 'ccle_utsw.lung_Mutation_allgene.csv')
    expr_path = os.path.join(data_dir, 'ccle_utsw.lung_RNAseq_allgene.csv')
    cnv_path = os.path.join(data_dir, 'ccle_utsw.lung_CNV_allgene.csv')
    drug_path = os.path.join(data_dir, 'ccle_utsw.lung_drug.split.csv')
    out_path = os.path.join(data_dir, 'ccle_utsw.lung_MutGeneExprCNV_allgene_drugRes.pkl')

    mut = pd.read_csv(mut_path, usecols=['Cell', 'Gene']).drop_duplicates()
    expr = pd.read_csv(expr_path, index_col=0)
    cnv = pd.read_csv(cnv_path, index_col=0)
    drug = pd.read_csv(drug_path)

    cell_id = list(set(drug.Cell) & set(mut['Cell']))
    gene_id = list(set(mut.columns) | (set(expr.columns) & set(cnv.columns)))
    
    # make mutation array
    mut['value'] = 1
    mut = mut.pivot(index='Cell', columns='Gene', values='value').fillna(0)
    mut = mut.loc[mut.index.to_series().isin(cell_id),:]
    print(mut.shape)

    # impute expr
    mean_expr = expr.mean(axis=0)
    expr = pd.concat([expr, pd.DataFrame(index=list(set(cell_id) - set(expr.index)))], sort=False).sort_index()
    expr = expr.apply(lambda x: mean_expr if all(x.isnull()) else x, axis=1)
    expr = expr.loc[expr.index.to_series().isin(cell_id),expr.columns.to_series().isin(gene_id)]
    print(expr.shape)

    # impute cnv
    mean_cnv = cnv.mean(axis=0)
    cnv = pd.concat([cnv, pd.DataFrame(index=list(set(cell_id) - set(cnv.index)))], sort=False).sort_index()
    cnv = cnv.apply(lambda x: mean_cnv if all(x.isnull()) else x, axis=1)
    cnv = cnv.loc[cnv.index.to_series().isin(cell_id),cnv.columns.to_series().isin(gene_id)]
    print(cnv.shape)

    expr.columns = ['expr_'+gene for gene in expr.columns.values]
    cnv.columns = ['cnv_'+gene for gene in cnv.columns.values]

    # concate mut, expr and cnv
    drug = drug.loc[drug.Cell.isin(cell_id),:]
    y = drug.LOG50
    drug_dummy = pd.get_dummies(drug.Drug)
    cell_mat = pd.concat([mut, expr, cnv], axis=1, sort=False)
    cell_mat = drug.Cell.apply(lambda x: cell_mat.loc[x,:])
    X = pd.concat([cell_mat, drug_dummy], axis=1, sort=False)
    print(X.head())
    print(X.shape)
    print(y.shape)

    # split train, test, val
    train_idx = drug.Train_split_point == 'Train'
    val_idx = drug.Train_split_point == 'Val'
    test_idx = drug.Train_split_point == 'Test'
    infer_idx = (drug.Train_split_point == 'Inference') & (~drug.Dropped_previously)
    infer_full_idx = drug.Train_split_point == 'Inference'

    X_train, y_train = X.loc[train_idx,:], y[train_idx]
    X_val, y_val = X.loc[val_idx,:], y[val_idx]
    X_test, y_test = X.loc[test_idx,:], y[test_idx]
    X_infer, y_infer = X.loc[infer_idx,:], y[infer_idx]
    X_infer_full, y_infer_full = X.loc[infer_full_idx,:], y[infer_full_idx]

    with open(out_path, 'wb') as f:
        pkl.dump({'train': (X_train, y_train), 
                'val': (X_val, y_val),
                'test': (X_test, y_test),
                'inference': (X_infer, y_infer),
                'inference_full': (X_infer_full, y_infer_full)}, file=f, protocol=4)