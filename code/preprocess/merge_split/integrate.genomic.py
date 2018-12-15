#####################################################################################
###                             integrate.genomic.py                              ###
#####################################################################################

import os
import re
import sys
import numpy as np
from scipy import sparse
import pandas as pd
import pickle as pkl
proj_dir = '/work/bioinformatics/s418336/projects/DLMed/'
data_dir = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version')
sys.path.append(os.path.join(proj_dir, 'code'))
from utility import integrate as itg

mut_path = os.path.join(data_dir, 'ccle_utsw.lung_Mutation_cancergene.loci.csv')
expr_path = os.path.join(data_dir, 'ccle_utsw.lung_RNAseq_cancergene.csv')
drug_path = os.path.join(data_dir, 'ccle_utsw.lung_drug.split.csv')
out_path = os.path.join(data_dir, 'ccle_utsw.lung_MutAFExpr_cancergene_drugRes.pkl')

###############################    main    ##################################
# read data
mut = pd.read_csv(mut_path)
expr = pd.read_csv(expr_path, index_col=0)
drug = pd.read_csv(drug_path)

# >>>>>>>>>>> drug data <<<<<<<<<<< #
# build cell & drug index dict
all_cells = set(drug['Cell']) & set(mut['Cell'])
cell_index = {cell:idx for idx, cell in enumerate(np.sort(list(all_cells)))}
drug = drug.loc[drug['Cell'].isin(cell_index.keys()),:]
mut = mut.loc[mut['Cell'].isin(cell_index.keys()),]
drug_index = {drug:idx for idx, drug in enumerate(np.unique(drug['Drug']))}
# drug array
drug['Cell'] = drug['Cell'].apply(lambda x: cell_index[x])
drug['Drug'] = drug['Drug'].apply(lambda x: drug_index[x])
source_dict = {'CCLE': 1, 'UTSW':2, 'CTRP': 3}
drug['Source'] = drug['Source'].apply(lambda x: source_dict[x])

# >>>>>>>>>>>>>> expr and cnv data <<<<<<<<<<<<<< #
# build gene index dict
genes = list(set(expr.columns))
subgenes = np.unique(mut['SubGene'])
gene_index = itg.buildGeneIndex(subgenes=subgenes, genes=genes)
# expr array
expr_array = itg.makeExprArray(expr, cell_index=cell_index, gene_index=gene_index, impute_na=True)

# >>>>>>>>>>>>>>> mutation data <<<<<<<<<<<<<<<< #
# make mutation table
mut_array = itg.makeMutationArray(mut, cell_index=cell_index, gene_index=gene_index)
print(mut_array.shape, expr_array.shape, drug.shape)

# write output
merged_lung_data = {'mutation_array': mut_array,
                    'expr_array': expr_array,
                    'drug_array': drug,
                    'cell_index': {value: key for key, value in cell_index.items()},
                    'gene_index': {value: key for key, value in gene_index.items()},
                    'drug_index': {value: key for key, value in drug_index.items()}}
with open(out_path, 'wb') as f:
    pkl.dump(merged_lung_data, file=f)
