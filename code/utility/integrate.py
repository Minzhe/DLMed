###################################################################################
###                               integrate.py                                    ###
###################################################################################

import re
import numpy as np
import pandas as pd


def makeMutationArray(df, cell_index, gene_index):
    '''
    Make cell-gene-mutation 3d array.
    '''
    mut_array = np.zeros(shape=(len(cell_index), len(gene_index), max(df['Loci'])+1), dtype=np.int8)
    for idx, row in df.iterrows():
        c_idx, g_idx, w_idx = cell_index[row['Cell']], gene_index[row['SubGene']], int(row['Loci'])
        mut_array[c_idx, g_idx, w_idx] = 1

    return mut_array

def makeExprArray(df, cell_index, gene_index):
    '''
    Make expression and copy number variation array.
    '''
    genes = np.unique(list(map(lambda x: re.sub(r'_[0-9]', '', x), gene_index.keys())))
    subgenes = list(gene_index.keys())
    max_sub = np.unique(list(map(lambda x: int(re.findall(r'_([0-9])', x)[0]), gene_index.keys())))

    mat = np.full(shape=(len(cell_index), len(gene_index)), fill_value=np.nan)
    # loop cells
    for cell in df.index:
        if cell not in cell_index.keys(): continue
        # loop genes
        for gene in df.columns:
            if gene not in genes: continue
            # loop sub genes
            for i in max_sub:
                gene_i = gene + '_' + str(i)
                if gene_i in subgenes:
                    mat[cell_index[cell], gene_index[gene_i]] = df.loc[cell, gene]

    return mat

def buildGeneIndex(subgenes, genes):
    '''
    Build subgene index with the common gene shared.
    '''
    subgenes = np.unique(subgenes)
    genes = np.unique(genes)
    genes = [gene + '_1' for gene in genes]
    all_genes = set(subgenes) | set(genes)
    gene_index = {gene:idx for idx, gene in enumerate(np.sort(list(all_genes)))}

    return gene_index