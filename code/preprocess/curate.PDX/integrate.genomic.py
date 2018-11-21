#################################################################################
###                         integrate.genomic.py                              ###
#################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import re
import sys
import numpy as np
import pandas as pd
import pickle as pkl
sys.path.append(os.path.join(proj_dir, 'code'))
from utility import integrate as itg
from pprint import pprint

mut_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.mutation.cancergenes.loci.csv')
expr_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.RPKM.cancergenes.csv')
cnv_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.CNV.cancergenes.csv')
drug_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.drug_response.csv')
out_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.MutExprCNV.cancergenes.drug_respnse.pkl')

################################    main    #################################
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

# >>>>>>>>>>>>>> expr and cnv data <<<<<<<<<<<<<< #
# build gene index dict
genes = list(set(expr.columns) & set(cnv.columns))
subgenes = np.unique(mut['Gene'])
gene_index = itg.buildGeneIndex(subgenes=subgenes, genes=genes)
# expr and cnv array
expr_array = itg.makeExprArray(expr, cell_index=sample_index, gene_index=gene_index)
cnv_array = itg.makeExprArray(cnv, cell_index=sample_index, gene_index=gene_index)

# >>>>>>>>>>>>>>> mutation data <<<<<<<<<<<<<<<< #
# make mutation table
mut_array = itg.makeMutationArray(mut, cell_index=sample_index, gene_index=gene_index)
print(mut_array.shape, expr_array.shape, cnv_array.shape, drug.shape)
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



