#########################################################################################
###                         integrate.multigenomic.loci.py                            ###
#########################################################################################
# integrate multi-genomic data in loci level

proj_path = '/work/bioinformatics/s418336/projects/DLMed'
import os
import re
import sys
import numpy as np
import pandas as pd
import pickle as pkl
sys.path.append(os.path.join(proj_path, 'code'))
from utility import integrate as itg
from pprint import pprint

mut_path = os.path.join(proj_path, 'data/curated/lung/utsw.mw/lung_Mutation_cancergene_loci.csv')
expr_path = os.path.join(proj_path, 'data/curated/lung/utsw.mw/lung_RNAseq_cancergene.csv')
cnv_path = os.path.join(proj_path, 'data/curated/lung/utsw.mw/lung_CNV_cancergene.csv')
drug_path = os.path.join(proj_path, 'data/curated/lung/utsw.mw/Lung_Compound_ED50.csv')
out_path = os.path.join(proj_path, 'data/curated/lung/utsw.mw/lung_MutExprCNV_cancergene_drug_loci.pkl')


###############################    main    ##################################
# read data
mut = pd.read_csv(mut_path)
expr = pd.read_csv(expr_path, index_col=0)
cnv = pd.read_csv(cnv_path, index_col=0)
drug = pd.read_csv(drug_path)

# >>>>>>>>>>> drug data <<<<<<<<<<< #
# build cell & drug index dict
all_cells = set(drug['Cell']) & set(mut['Cell'])
cell_index = {cell:idx for idx, cell in enumerate(np.sort(list(all_cells)))}
drug = drug.loc[drug['Cell'].isin(cell_index.keys()),:]
mut = mut.loc[mut['Cell'].isin(cell_index.keys()),]
drug_index = {drug:idx for idx, drug in enumerate(np.unique(drug['Compound']))}
# drug array
drug['Cell'] = drug['Cell'].apply(lambda x: cell_index[x])
drug['Compound'] = drug['Compound'].apply(lambda x: drug_index[x])
drug_array = np.array(drug[['Cell', 'Compound', 'ED50']])


# >>>>>>>>>>>>>> expr and cnv data <<<<<<<<<<<<<< #
# build gene index dict
genes = list(set(expr.columns) & set(cnv.columns))
subgenes = np.unique(mut['SubGene'])
gene_index = itg.buildGeneIndex(subgenes=subgenes, genes=genes)
# expr and cnv array
expr_array = itg.makeExprArray(expr, cell_index=cell_index, gene_index=gene_index)
cnv_array = itg.makeExprArray(cnv, cell_index=cell_index, gene_index=gene_index)

# >>>>>>>>>>>>>>> mutation data <<<<<<<<<<<<<<<< #
# make mutation table
mut_array = itg.makeMutationArray(mut, cell_index=cell_index, gene_index=gene_index)
print(mut_array.shape, expr_array.shape, cnv_array.shape, len(drug_index))
exit()

# write output
utsw_lung_data = {'mutation_array': mut_array,
                  'expr_array': expr_array,
                  'cnv_array': cnv_array,
                  'drug_array': drug_array,
                  'cell_index': {value: key for key, value in cell_index.items()},
                  'gene_index': {value: key for key, value in gene_index.items()},
                  'drug_index': {value: key for key, value in drug_index.items()}}
with open(out_path, 'wb') as f:
    pkl.dump(utsw_lung_data, file=f)
