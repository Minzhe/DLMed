#########################################################################################
###                         integrate.multigenomic.gene.py                            ###
#########################################################################################
# integrate multi-genomic data in gene level

import os
import functools
import pandas as pd
import numpy as np
from sklearn.preprocessing import scale


proj_path = 'D:/projects/DLCell'
lung_expr_path = os.path.join(proj_path, 'data/curated/lung/lung_RNAseq_cancergene.csv')
lung_mutation_path = os.path.join(proj_path, 'data/curated/lung/lung_Mutation_cencergene.csv')
lung_cnv_path = os.path.join(proj_path, 'data/curated/lung/lung_CNV_cancergene.csv')
lung_drug_path = os.path.join(proj_path, 'data/curated/lung/ccle.lung.gdsc.csv')


###############################    function    #################################
def concatMD(**genomedata):
    '''
    Concatenate genomic data to multiple dimension array.
    '''
    # filter common genes and cell
    datalist = list(genomedata.values())

    celllist = list(map(lambda x: x.index.tolist(), datalist))
    genelist = list(map(lambda x: x.columns.tolist(), datalist))
    cells = functools.reduce(lambda x,y: set(x) & set(y), celllist)
    genes = functools.reduce(lambda x,y: set(x) & set(y), genelist)

    datalist = list(map(lambda x: np.array(x.loc[cells,genes]), datalist))

    # create multi-dimensional array
    return np.array(datalist)

def concat1D(**genomedata):
    '''
    Concatenate genomic data to one dimension array.
    '''
    for name in genomedata.keys():
        genomedata[name].columns = list(map(lambda x: '_'.join([name, x]), genomedata[name].columns.tolist()))
    data_mat = pd.concat(list(genomedata.values()), join='inner', axis=1)
    data_mat.reset_index(level=0, inplace=True)
    data_mat.columns.values[0] = 'Cell'

    return data_mat

def concatGenomicDrug(genomedata, drugdata):
    '''
    Concat genomic data with drug IC50 data.
    '''
    genome_drug = genome_data.merge(drugdata, left_on='Cell', right_on='Aliases', how='inner')
    genome_drug.drop(['Aliases'], inplace=True, axis=1)
    dummy_drug = pd.get_dummies(genome_drug.Drug)
    dummy_drug.columns = list(map(lambda x: 'Drug_' + x, dummy_drug.columns.values))
    genome_drug = pd.concat([genome_drug.iloc[:,:-2], dummy_drug, genome_drug.iloc[:,-1]], axis=1)
    return genome_drug


###############################    main    #################################
expr = pd.read_csv(lung_expr_path, index_col=0)
expr = expr.apply(scale, axis=0)
mutation = pd.read_csv(lung_mutation_path, index_col=0)
cnv = pd.read_csv(lung_cnv_path, index_col=0)
cnv = cnv.apply(scale, axis=0)

# common cell and gene
# genome_data = concat1D(expr=expr, mutation=mutation, cnv=cnv)
genome_data = concat1D(cnv=cnv)

# drug data
drug = pd.read_csv(lung_drug_path, usecols=['Aliases', 'Drug', 'IC50'])

# concat genomic and drug data
data_mat = concatGenomicDrug(genome_data, drug)
# data_mat.to_csv(os.path.join(proj_path, 'data/curated/lung/lung.genomic.drug.csv'), index=None)
data_mat.to_csv(os.path.join(proj_path, 'data/curated/lung/lung.cnv.drug.csv'), index=None)


