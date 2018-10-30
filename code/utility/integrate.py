###################################################################################
###                               integrate.py                                    ###
###################################################################################

import re
import numpy as np
import pandas as pd
import functools

# =========================  convert dataframe to array for training  ========================== #
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


# =========================  merge multiple data source  ========================== #
def mergeMutation(**mut):
    '''
    Merge mutation table
    '''
    mut_all = list()
    for df in mut.keys():
        if all(pd.Series(['Cell', 'Gene', 'MutStart']).isin(mut[df].columns.tolist())) is not True:
            raise KeyError('Input {} dateframe must contains column: Cell, Gene and MutStart!'.format(df))
        mut[df]['Source'] = df
        mut_all.append(mut[df][['Cell', 'Gene', 'MutStart', 'Source']])
    mut_all = pd.concat(mut_all)

    # merge
    mut_all = mut_all.groupby(by=['Cell', 'Gene', 'MutStart'], as_index=False).agg({'Source': lambda x: ','.join(x)})
    mut_all.sort_values(by=['Cell', 'Gene', 'MutStart'], inplace=True)
    
    return mut_all

def mergeExpression(merge_seq, gene_join='inner', fillna=0, **expr):
    '''
    Merge expression and cnv data.
    '''
    genes = list(map(lambda x: x.columns.tolist(), expr.values()))
    # fill all data with union gene set
    if gene_join == 'outer':
        all_genes = functools.reduce(lambda x,y: set(x) | set(y), genes)
        for name in expr.keys():
            expr[name] = pd.concat([expr[name], pd.DataFrame(columns=list(all_genes-set(expr[name].columns)))], axis=1, sort=True).fillna(0)
    # select common genes
    elif gene_join == 'inner':
        all_genes = functools.reduce(lambda x,y: set(x) & set(y), genes)
        for name in expr.keys():
            expr[name] = expr[name].loc[:,expr[name].columns.to_series().isin(list(all_genes))]
    # merge according to the importance sequence
    expr_all = pd.DataFrame()
    for name in merge_seq:
        tmp_expr = expr[name].loc[list(set(expr[name].index) - set(expr_all.index)),:]
        expr_all = pd.concat([expr_all, tmp_expr], sort=True)
    expr_all.sort_index(inplace=True)
    return expr_all

def mergeDrug(duplicate='keep', **drug):
    '''
    Merge drug sensitivity data.
    '''
    for name in drug.keys():
        drug[name]['Source'] = name
    drug_all = pd.concat(drug.values())
    drug_all.sort_values(by=['Cell', 'Drug', 'Source'], inplace=True)
    drug_all.index = list(range(drug_all.shape[0]))
    dup_idx = drug_all.duplicated(subset=['Cell', 'Drug'], keep=False)
    if duplicate == 'keep':
        return drug_all, drug_all.loc[dup_idx,:]

    
    