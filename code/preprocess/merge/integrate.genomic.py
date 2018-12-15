#######################################################################################
###                             integrate.genomic.py                                ###
#######################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed/'
import os
import sys
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
sys.path.append(os.path.join(proj_dir, 'code'))
from sklearn.preprocessing import scale

#######################################    main    #############################################
input_data = 'mut_expr_cnv_fix_cellline'

if input_data == 'cell_onehot':
    # cell one hot encode
    file_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_drug.csv')
    out_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_drug.ccle.onehot.pkl')

    drug_data = pd.read_csv(file_path)
    drug_data = drug_data.loc[drug_data.Source.str.contains('CCLE'),['Cell', 'Drug', 'LOG50']]

    X = pd.get_dummies(drug_data.iloc[:,[0,1]])
    y = drug_data['LOG50']

    with open(out_path, 'wb') as f:
        pkl.dump({'X': X, 'y': y}, file=f)

elif input_data == 'mutation':
    # mutation
    file_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_Mut_cancergene_drug_loci.array.pkl')
    out_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_Mut_cancergene_drug_loci.array.subset.pkl')

    with open(file_path, 'rb') as f:
        data = pkl.load(f)
    
    X, y = data['X_csr'].toarray(), data['y']
    print(X.shape)
    print(pd.Series(np.apply_along_axis(sum, 0, X)).value_counts())

elif input_data == 'expression':
    # expression and drug dummy
    file_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_MutExprCNV_cancergene_drug_loci.untruncated.pkl')
    out_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_Expr_cancergene_drug.array.pkl')

    with open(file_path, 'rb') as f:
        data = pkl.load(f)

    # recover expr
    expr = pd.DataFrame(data['expr_array'], index=data['cell_index'].values(), columns=data['gene_index'].values()).dropna(axis=0, how='all').dropna(axis=1, how='any')
    drugic50 = pd.DataFrame(data['drug_array'], columns=['Cell', 'Drug', 'LOG50', 'Source'])
    drugic50.Cell = drugic50.Cell.apply(lambda x: data['cell_index'][x])
    drugic50.Drug = drugic50.Drug.apply(lambda x: data['drug_index'][x])
    drugic50 = drugic50.loc[drugic50.Cell.isin(expr.index.values),:]
    # get dummy
    cell_dummy = drugic50.Cell.apply(lambda x: expr.loc[x,:])
    drug_dummy = pd.get_dummies(drugic50.Drug)
    X = pd.concat([cell_dummy, drug_dummy], axis=1)
    y = drugic50.LOG50
    label = drugic50.Source

    with open(out_path, 'wb') as f:
        pkl.dump({'X': X, 'y': y, 'label': label}, file=f)

elif input_data == 'mut_expr_cnv':
    # mutation (gene level), expression, cnv and drug dummmy
    mut_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_Mutation_cancergene.array.csv')
    expr_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_RNAseq_cancergene.csv')
    cnv_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_CNV_cancergene.csv')
    drug_path = os.path.join(proj_dir, 'data/curated/Lung/merged/ccle_utsw.lung_drug.csv')
    out_path = os.path.join(proj_dir, 'data/curated/Lung/merged/ccle_utsw.lung_MutExprCNV_cancergene_drug.array.pkl')

    mut = pd.read_csv(mut_path, index_col=0)
    expr = pd.read_csv(expr_path, index_col=0); expr.columns = ['expr_'+gene for gene in expr.columns.values]
    cnv = pd.read_csv(cnv_path, index_col=0); cnv.columns = ['cnv_'+gene for gene in cnv.columns.values]
    expr = pd.DataFrame(scale(expr), index=expr.index, columns=expr.columns)
    cnv = pd.DataFrame(scale(cnv), index=cnv.index, columns=cnv.columns)
    drugic50 = pd.read_csv(drug_path)

    # filter cell lines
    all_cells = list(set(mut.index) & set(drugic50['Cell']))
    drugic50 = drugic50.loc[drugic50['Cell'].isin(all_cells),:]
    mut = mut.loc[mut.index.to_series().isin(all_cells),:]
    expr = expr.loc[expr.index.to_series().isin(all_cells),:]
    mean_expr = expr.mean(axis=0)
    expr = pd.concat([expr, pd.DataFrame(index=list(set(all_cells)-set(expr.index)))]).sort_index()
    expr = expr.apply(lambda x: mean_expr if all(x.isnull()) else x, axis=1)
    cnv = cnv.loc[cnv.index.to_series().isin(all_cells),:]
    cnv = pd.concat([cnv, pd.DataFrame(index=list(set(all_cells)-set(cnv.index)))]).sort_index().fillna(0)

    # concate mut, expr and cnv
    y = np.array(drugic50.LOG50)
    drug_dummy = pd.get_dummies(drugic50.Drug)
    cell_mat = pd.concat([mut, expr, cnv], axis=1, sort=False)
    cell_mat = drugic50.Cell.apply(lambda x: cell_mat.loc[x,:])
    X = pd.concat([cell_mat, drug_dummy], axis=1, sort=False)
    
    with open(out_path, 'wb') as f:
        pkl.dump({'X': X, 'y': y}, file=f)

elif input_data == 'mut_expr_cnv_fix_cellline':
    mut_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_Mutation_cancergene.array.csv')
    expr_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_RNAseq_cancergene.csv')
    cnv_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_CNV_cancergene.csv')
    train_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_MutExprCNV_cancergene_drug_loci.untruncated.label_12_70_20_10_35324_train.csv')
    val_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_MutExprCNV_cancergene_drug_loci.untruncated.label_12_70_20_10_9120_val.csv')
    test_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_MutExprCNV_cancergene_drug_loci.untruncated.label_12_70_20_10_5608_inference.csv')
    out_path = os.path.join(proj_dir, 'data/curated/Lung/merged/ccle_utsw.lung_MutExprCNV_cancergene_drug.array.fixcell.70.20.10.pkl')

    mut = pd.read_csv(mut_path, index_col=0)
    expr = pd.read_csv(expr_path, index_col=0); expr.columns = ['expr_'+gene for gene in expr.columns.values]
    cnv = pd.read_csv(cnv_path, index_col=0); cnv.columns = ['cnv_'+gene for gene in cnv.columns.values]
    expr = pd.DataFrame(scale(expr), index=expr.index, columns=expr.columns)
    cnv = pd.DataFrame(scale(cnv), index=cnv.index, columns=cnv.columns)

    train_data = pd.read_csv(train_path, index_col=0)
    val_data = pd.read_csv(val_path, index_col=0)
    test_data = pd.read_csv(test_path, index_col=0)
    train_index = train_data.index
    val_index = val_data.index
    test_index = test_data.index
    drugic50 = pd.concat([train_data, val_data, test_data])

    cell_id = list(drugic50.cell_id)
    mut = mut.loc[mut.index.to_series().isin(cell_id)]
    expr = expr.loc[expr.index.to_series().isin(cell_id),:]
    mean_expr = expr.mean(axis=0)
    expr = pd.concat([expr, pd.DataFrame(index=list(set(cell_id)-set(expr.index)))], sort=False).sort_index()
    expr = expr.apply(lambda x: mean_expr if all(x.isnull()) else x, axis=1)
    cnv = cnv.loc[cnv.index.to_series().isin(cell_id),:]
    cnv = pd.concat([cnv, pd.DataFrame(index=list(set(cell_id)-set(cnv.index)))], sort=False).sort_index().fillna(0)

    # concate mut, expr and cnv
    y = drugic50.score
    drug_dummy = pd.get_dummies(drugic50.drug_id)
    cell_mat = pd.concat([mut, expr, cnv], axis=1, sort=False)
    cell_mat = drugic50.cell_id.apply(lambda x: cell_mat.loc[x,:])
    X = pd.concat([cell_mat, drug_dummy], axis=1, sort=False)

    X_train, y_train = X.loc[train_index,:], y[train_index]
    X_val, y_val = X.loc[val_index,:], y[val_index]
    X_test, y_test = X.loc[test_index,:], y[test_index]
    
    with open(out_path, 'wb') as f:
        pkl.dump({'train': (X_train, y_train), 
                  'val': (X_val, y_val),
                  'test': (X_test, y_test)}, file=f)
