###########################################################################
###                         curate.RNASeq.py                            ###
###########################################################################

import numpy as np
import pandas as pd
import os
import re


############################   function   ##########################
def filterRNASeq(rnaseq, tissue, genelist):
    '''
    Select specific tissue cell lines.
    '''
    ### filter tissue
    tissue = tissue.upper()
    expr = rnaseq.loc[:,np.insert(rnaseq.columns[2:].str.contains(tissue), 0, [True, True])]
    expr.columns = list(map(lambda x: re.sub(r'_{}.*'.format(tissue), '', x), expr.columns.values))
    expr.columns.values[:2] = ['ENSG', 'Gene']
    expr = expr.iloc[:,1:]
    expr.Gene = expr.Gene.astype(str)
    print('  Number of cell lines: {}'.format(len(expr.columns)-1))

    ### filter cancer genes
    expr = expr.loc[expr.Gene.isin(genelist),:]
    expr.sort_values(by=['Gene'], inplace=True)
    expr.set_index('Gene', inplace=True)

    return expr.transpose()


proj_path = "D:/projects/DLCell"
rnaseq_path = os.path.join(proj_path, "data/DepMap/CCLE_RNAseq_RPKM_20180718.gct")
geneanno_path = os.path.join(proj_path, "data/curated/cancer.gene.anno.csv")


############################   main   ################################
# read cancer gene list
anno = pd.read_csv(geneanno_path)['Gene'].values

# read gene expression data
ccle_rnaseq = pd.read_csv(rnaseq_path, sep='\t', header=0, skiprows=2)

# filter for lung cancer cell lines
print('Filtering lung cancer cell line ...')
lung_expr = filterRNASeq(rnaseq=ccle_rnaseq, tissue='LUNG', genelist=anno)
lung_expr.to_csv(os.path.join(proj_path, 'data/curated/lung_RNAseq_cancergene.csv'))

# # filter for breast cancer cell lines
# print('Filtering breast cancer cell lines ...')
# bc_expr = filterRNASeq(rnaseq=ccle_rnaseq, tissue='BREAST', genelist=anno)
# bc_expr.to_csv(os.path.join(proj_path, "data/curated/breast_RNAseq_cancergene.csv"), index=None)

# # filter for kidney cancer cell lines
# print('Filtering kidney cancer cell lines ...')
# kidney_expr = filterRNASeq(rnaseq=ccle_rnaseq, tissue='KIDNEY', genelist=anno)
# kidney_expr.to_csv(os.path.join(proj_path, "data/curated/kidney_RNAseq_cancergene.csv"), index=None)

# # filter for leukemia cancer cell lines
# print('Filtering leukemia cancer cell lines ...')
# leukemia_expr = filterRNASeq(rnaseq=ccle_rnaseq, tissue='HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', genelist=anno)
# leukemia_expr.to_csv(os.path.join(proj_path, "data/curated/leukemia_RNAseq_cancergene.csv"), index=None)