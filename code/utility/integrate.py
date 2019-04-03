###################################################################################
###                               integrate.py                                    ###
###################################################################################

import re
import numpy as np
# np.set_printoptions(threshold=np.nan)
import pandas as pd
import functools

# =========================  convert dataframe to array for training  ========================== #
def makeMutationArray(df, cell_index, gene_index):
    '''
    Make cell-gene-mutation 3d array.
    '''
    mut_array = np.zeros(shape=(len(cell_index), len(gene_index), max(df['Loci'])+1), dtype=np.float)
    for idx, row in df.iterrows():
        c_idx, g_idx, w_idx = cell_index[row[0]], gene_index[row[1]], int(row[2])
        mut_array[c_idx, g_idx, w_idx] = round(row[3], 3)

    return mut_array

def makeExprArray(df, cell_index, gene_index, impute_na=True, impute_value=0):
    '''
    Make expression and copy number variation array.
    '''
    mat = df.copy()
    mat.columns = mat.columns.to_series().apply(lambda x: x + '_0')
    for sub_gene in gene_index.keys():
        if sub_gene not in list(mat.columns):
            gene = re.sub(r'_[0-9]+', '', sub_gene)
            if gene + '_0' in list(mat.columns):
                mat[sub_gene] = mat[gene + '_0']
            else:
                mat[sub_gene] = np.nan

    for cell in cell_index.keys():
        if cell not in list(mat.index):
            mat.loc[cell,:] = np.nan
    
    mat = mat.loc[:,mat.columns.to_series().isin(list(gene_index.keys()))]
    mat.columns = list(map(lambda x: gene_index[x], list(mat.columns)))
    mat.sort_index(axis=1, inplace=True)

    mat = mat.loc[mat.index.to_series().isin(list(cell_index.keys())),:]
    mat.index = list(map(lambda x: cell_index[x], list(mat.index)))
    mat.sort_index(axis=0, inplace=True)

    mat = np.array(mat)
    # check nan value
    r_cell_index = {value: key for key, value in cell_index.items()}
    r_gene_index = {value: key for key, value in gene_index.items()}
    r_idx = np.where(np.apply_along_axis(all, 1, np.isnan(mat)))[0]
    c_idx = np.where(np.apply_along_axis(all, 0, np.isnan(mat)))[0]
    if impute_na:
        if len(r_idx) != 0:
            print('Warning: {} row - {} is all nan values.'.format(r_idx, [r_cell_index[x] for x in r_idx]))
            col_mean = np.apply_along_axis(np.nanmean, 0, mat)
            mat[r_idx,:] = col_mean
        if len(c_idx) != 0:
            print('Warning: {} column - {} is all nan values.'.format(c_idx, [r_gene_index[x] for x in c_idx]))
            mat[:,c_idx] = impute_value

    return mat

def buildGeneIndex(subgenes, genes, join='inner'):
    '''
    Build subgene index with the common gene shared.
    '''
    subgenes = np.unique(subgenes)
    genes = np.unique(genes)
    genes = [gene + '_0' for gene in genes]
    if join == 'inner':
        all_genes = set(subgenes) & set(genes)
    elif join == 'outer':
        all_genes = set(subgenes) | set(genes)
    elif join == 'left':
        all_genes = set(subgenes)
    elif join == 'right':
        all_genes = set(genes)
    gene_index = {gene:idx for idx, gene in enumerate(np.sort(list(all_genes)))}

    return gene_index


# =========================  merge multiple data source  ========================== #
def mergeMutation(**mut):
    '''
    Merge mutation table
    '''
    def mergeAF(af):
        if len(af) == sum(af == 0.5):
            return 0.5
        else:
            return af[af != 0.5].mean()
    
    mut_all = list()
    for df in mut.keys():
        if all(pd.Series(['Cell', 'Gene', 'MutStart', 'AF']).isin(mut[df].columns.tolist())) is not True:
            raise KeyError('Input {} dateframe must contains column: Cell, Gene, MutStart and AF!'.format(df))
        mut[df]['Source'] = df
        mut_all.append(mut[df][['Cell', 'Gene', 'MutStart', 'AF', 'Source']])
    mut_all = pd.concat(mut_all)

    # merge
    mut_all = mut_all.groupby(by=['Cell', 'Gene', 'MutStart'], as_index=False).agg({'AF': mergeAF, 'Source': lambda x: ','.join(x)})
    mut_all.sort_values(by=['Cell', 'Gene', 'MutStart'], inplace=True)
    
    return mut_all

def mergeExpression(keep_dup, gene_join='inner', fillna=None, **expr):
    '''
    Merge expression and cnv data.
    '''
    genes = list(map(lambda x: x.columns.tolist(), expr.values()))
    # fill all data with union gene set
    if gene_join == 'outer':
        all_genes = functools.reduce(lambda x,y: set(x) | set(y), genes)
        for name in expr.keys():
            expr[name] = pd.concat([expr[name], pd.DataFrame(columns=list(all_genes-set(expr[name].columns)))], axis=1, sort=True).fillna(fillna)
    # select common genes
    elif gene_join == 'inner':
        all_genes = functools.reduce(lambda x,y: set(x) & set(y), genes)
        for name in expr.keys():
            expr[name] = expr[name].loc[:,expr[name].columns.to_series().isin(list(all_genes))]
    # merge according to the importance sequence
    expr_all = pd.DataFrame()
    for name in keep_dup:
        tmp_expr = expr[name].loc[list(set(expr[name].index) - set(expr_all.index)),:]
        expr_all = pd.concat([expr_all, tmp_expr], sort=True)
    expr_all.sort_index(inplace=True)
    return expr_all

def mergeDrug(keep_dup='all', **drug):
    '''
    Merge drug sensitivity data.
    '''
    def prior_select(data, prior):
        for item in prior:
            if item in list(data):
                return data == item
        return pd.Series([False] * len(data))
    # formate data
    for name in drug.keys():
        drug[name]['Source'] = name
    drug_all = pd.concat(drug.values())
    drug_all.sort_values(by=['Cell', 'Drug', 'Source'], inplace=True)
    drug_all.index = list(range(drug_all.shape[0]))
    # duplicate
    if keep_dup == 'all':
        pass
    elif all(pd.Series(keep_dup).isin(drug.keys())):
        idx = drug_all.groupby(by=['Cell', 'Drug'], sort=False, as_index=False)['Source'].apply(lambda x: prior_select(x, keep_dup))
        drug_all = drug_all.loc[list(idx),:]
    else:
        raise ValueError('Unrecognizable keep duplicate actions.')
    drug_all.index = list(range(drug_all.shape[0]))

    return drug_all

    
    