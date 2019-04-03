#################################################################################
###                         integrate.genomic.py                              ###
#################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import re
import sys
import numpy as np
np.set_printoptions(threshold=np.nan)
from scipy import sparse
import pandas as pd
import pickle as pkl
sys.path.append(os.path.join(proj_dir, 'code'))
from utility import integrate as itg
from pprint import pprint
import matplotlib.pyplot as plt
import utility.plot as p
import utility.utility as util

drug_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.drug_response.csv')
gene_set = 'combined'

################################    main    #################################
if gene_set == 'cancer_gene':
    mut_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.mutation.cancergenes.loci.csv')
    expr_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.RPKM.cancergenes.csv')
    cnv_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.CNV.cancergenes.csv')
    out_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.MutAFExprCNV.cancergenes.drug_response.pkl')

    mut = pd.read_csv(mut_path)
    expr = pd.read_csv(expr_path, index_col=0)
    cnv = pd.read_csv(cnv_path, index_col=0)
    drug = pd.read_csv(drug_path)

    # >>>>>>>>>>> drug data <<<<<<<<<<< #
    # build cell index dict
    all_samples = set(drug['Model']) & set(mut['Sample'])
    sample_index = {sample: idx for idx, sample in enumerate(np.sort(list(all_samples)))}
    drug = drug.loc[drug['Model'].isin(sample_index.keys()),:]
    mut = mut.loc[mut['Sample'].isin(sample_index.keys()),]
    # drug index
    drug_index = {drug: idx for idx, drug in enumerate(np.unique(drug['Treatment']))}
    # response index
    response_index = {resp: idx for idx, resp in enumerate(np.unique(drug['ResponseCategory']))}
    # drug array
    drug['Model'] = drug['Model'].apply(lambda x: sample_index[x])
    drug['Treatment'] = drug['Treatment'].apply(lambda x: drug_index[x])
    drug['ResponseCategory'] = drug['ResponseCategory'].apply(lambda x: response_index[x])
    # remove models with small data points
    data_counts = drug.Model.value_counts()
    rm_idx = list(data_counts[data_counts < 10].sort_index().index)
    drug = drug.loc[~drug.Model.isin(rm_idx),:]
    drug.index = list(range(drug.shape[0]))
    # split data
    drug['Train_split'] = util.split_train_val_test(index=drug.Treatment, train_prop=0.7, val_prop=0.2)
    # plot heatmap
    # drug_mat = drug.pivot(index='Model', columns='Treatment', values='TimeToDouble').fillna(0)
    # drug_mat = np.log(drug_mat + 1)
    # f, ax = plt.subplots(figsize=(20,25))
    # p.plotHeatmap(drug_mat, mask_value=0)
    # plt.savefig(os.path.join(os.path.join(proj_dir, 'data/curated/PDX/PDX.drug_respnse.TimeToDouble.png')))

    # >>>>>>>>>>>>>> expr and cnv data <<<<<<<<<<<<<< #
    # build gene index dict
    genes = list(set(expr.columns))
    subgenes = np.unique(mut['Gene'])
    gene_index = itg.buildGeneIndex(subgenes=subgenes, genes=genes)
    # expr and cnv array
    expr_array = itg.makeExprArray(expr, cell_index=sample_index, gene_index=gene_index)
    cnv_array = itg.makeExprArray(cnv, cell_index=sample_index, gene_index=gene_index)

    # >>>>>>>>>>>>>>> mutation data <<<<<<<<<<<<<<<< #
    # make mutation table
    mut_array = itg.makeMutationArray(mut, cell_index=sample_index, gene_index=gene_index)
    print(mut_array.shape, expr_array.shape, drug.shape)
    pdx_data = {'mutation_array': mut_array,
                'expr_array': expr_array,
                'cnv_array': cnv_array,
                'drug_array': drug,
                'cell_index': {value: key for key, value in sample_index.items()},
                'gene_index': {value: key for key, value in gene_index.items()},
                'drug_index': {value: key for key, value in drug_index.items()},
                'response_index': {value: key for key, value in response_index.items()}}
    with open(out_path, 'wb') as f:
        pkl.dump(pdx_data, file=f)

elif gene_set == 'all_gene':
    mut_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.mutation.allgenes.loci.csv')
    expr_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.RPKM.allgenes.csv')
    cnv_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.CNV.allgenes.csv')
    out_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.MutAFExprCNV.allgenes.drug_response.pkl')

    mut = pd.read_csv(mut_path)
    expr = pd.read_csv(expr_path, index_col=0)
    cnv = pd.read_csv(cnv_path, index_col=0)
    drug = pd.read_csv(drug_path)

    # >>>>>>>>>>> drug data <<<<<<<<<<< #
    # build cell index dict
    all_samples = set(drug['Model']) & set(mut['Sample'])
    sample_index = {sample: idx for idx, sample in enumerate(np.sort(list(all_samples)))}
    drug = drug.loc[drug['Model'].isin(sample_index.keys()),:]
    mut = mut.loc[mut['Sample'].isin(sample_index.keys()),]
    # drug index
    drug_index = {drug: idx for idx, drug in enumerate(np.unique(drug['Treatment']))}
    # response index
    response_index = {resp: idx for idx, resp in enumerate(np.unique(drug['ResponseCategory']))}
    # drug array
    drug['Model'] = drug['Model'].apply(lambda x: sample_index[x])
    drug['Treatment'] = drug['Treatment'].apply(lambda x: drug_index[x])
    drug['ResponseCategory'] = drug['ResponseCategory'].apply(lambda x: response_index[x])
    # remove models with small data points
    data_counts = drug.Model.value_counts()
    rm_idx = list(data_counts[data_counts < 10].sort_index().index)
    drug = drug.loc[~drug.Model.isin(rm_idx),:]
    drug.index = list(range(drug.shape[0]))
    # split data
    drug['Train_split'] = util.split_train_val_test(index=drug.Treatment, train_prop=0.7, val_prop=0.2)

    # >>>>>>>>>>>>>> expr and cnv data <<<<<<<<<<<<<< #
    # build gene index dict
    genes = list(set(expr.columns) & set(cnv.columns))
    subgenes = np.unique(mut['Gene'])
    gene_index = itg.buildGeneIndex(subgenes=subgenes, genes=genes)
    # expr and cnv array
    expr_array = itg.makeExprArray(expr, cell_index=sample_index, gene_index=gene_index, impute_na=True, impute_value=0)
    cnv_array = itg.makeExprArray(cnv, cell_index=sample_index, gene_index=gene_index, impute_na=True, impute_value=2)

    # >>>>>>>>>>>>>>> mutation data <<<<<<<<<<<<<<<< #
    # make mutation table
    mut_array = itg.makeMutationArray(mut, cell_index=sample_index, gene_index=gene_index)
    print(mut_array.shape, expr_array.shape, drug.shape)
    pdx_data = {'mutation_array_allgenes': mut_array,
                'expr_array_allgenes': expr_array,
                'cnv_array_allgenes': cnv_array,
                'drug_array': drug,
                'cell_index_allgenes': {value: key for key, value in sample_index.items()},
                'gene_index_allgenes': {value: key for key, value in gene_index.items()},
                'drug_index': {value: key for key, value in drug_index.items()},
                'response_index': {value: key for key, value in response_index.items()}}
    with open(out_path, 'wb') as f:
        pkl.dump(pdx_data, file=f, protocol=4)

elif gene_set == 'combined':
    cancer_gene_data_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.MutAFExprCNV.cancergenes.drug_response.pkl')
    all_gene_data_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.MutAFExprCNV.allgenes.drug_response.pkl')
    with open(all_gene_data_path, 'rb') as f_all, open(cancer_gene_data_path, 'rb') as f_cancer:
        pdx_data_all = pkl.load(f_all)
        pdx_data_cancer = pkl.load(f_cancer)
        pdx_data_cancer['drug_array'] = pdx_data_all['drug_array']
        pdx_data_cancer['drug_index'] = pdx_data_all['drug_index']
        pdx_data_cancer['response_index'] = pdx_data_all['response_index']
    
    with open(cancer_gene_data_path, 'wb') as f:
        pkl.dump(pdx_data_cancer, file=f)




