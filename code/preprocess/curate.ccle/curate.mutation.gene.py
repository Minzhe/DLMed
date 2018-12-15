##################################################################################
###                         curate.mutation.gene.py                            ###
##################################################################################

import pandas as pd
import numpy as np
import os
import re

proj_path = "D:/projects/DLCell"
mutation_path = os.path.join(proj_path, "data/DepMap/CCLE_Mutation_20180718.txt")
rnaseq_lung_path = os.path.join(proj_path, "data/curated/lung_RNAseq_cancergene.csv")
geneanno_path = os.path.join(proj_path, "data/curated/cancer.gene.anno.csv")
cellline_path = os.path.join(proj_path, "data/curated/lung/celllines.lung.csv")

##############################  function  ################################
def makeMutationTable(mutations, celllist, genelist):
    '''
    Encode binary mutation table
    '''
    mut_mat = pd.DataFrame(0, index=np.sort(celllist), columns=np.sort(genelist))
    for idx, row in mutations.iterrows():
        mut_mat.loc[row['Cell'],row['Gene']] = 1

    return mut_mat


##############################  main  ################################
# lung cancer cell line
cancer_genes = pd.read_csv(geneanno_path)['Gene'].values
lung_cell = pd.read_csv(cellline_path, usecols=['Aliases'], squeeze=True).tolist()

# filter for lung cancer cell
mutation_all = pd.read_csv(mutation_path, sep='\t', usecols=['Hugo_Symbol', 'Variant_Classification', 'Tumor_Sample_Barcode'])
mutation_all.columns = ['Gene', 'MutClass', 'Cell']
mutations = mutation_all.loc[mutation_all.Cell.str.contains('LUNG'),:]
mutations['Cell'] = mutations.Cell.apply(lambda x: re.sub(r'_LUNG.*', '', x))
mutations = mutations.loc[mutations.Cell.isin(lung_cell),:]

# filter for cancer genes
mutations = mutations.loc[mutations.Gene.isin(cancer_genes),:]
mutations.sort_values(by=['Cell', 'Gene'], inplace=True)

# filter for mutation type
mutations = mutations.loc[~mutations.MutClass.isin(['Silent', 'Intron']),:]

# encode for binary mutation table
mut_mat = makeMutationTable(mutations, celllist=lung_cell, genelist=cancer_genes)
print(mut_mat.apply(sum, axis=1))
mut_mat.to_csv(os.path.join(proj_path, "data/curated/lung/lung_Mutation_cencergene.csv"))