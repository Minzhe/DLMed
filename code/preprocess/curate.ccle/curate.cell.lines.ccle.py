##################################################################################################
###                               curate.cell.lines.ccle.py                                    ###
##################################################################################################
# This script is to summarize the cell line information contained in different datasets.

import os
import re
import numpy as np
import pandas as pd
import utility as util

proj_dir = "D:/projects/DLCell"
cellline_path = os.path.join(proj_dir, 'data/DepMap/DepMap.celllines.2018.q3.csv')
expr_path = os.path.join(proj_dir, 'data/DepMap/CCLE_RNAseq_RPKM_20180718.gct')
mutation_path = os.path.join(proj_dir, 'data/DepMap/CCLE_Mutation_20180718.txt')
cnv_path = os.path.join(proj_dir, 'data/DepMap/DepMap.CNV.2018.q3.csv')
rppa_path = os.path.join(proj_dir, 'data/DepMap/CCLE_RPPA_20180123.csv')
rnai_path = os.path.join(proj_dir, 'data/DepMap/DepMap.RNAi.score.csv')
crispr_path = os.path.join(proj_dir, 'data/DepMap/DepMap.CRISPRi.score.csv')
gdsc_path = os.path.join(proj_dir, 'data/Cosmic/GDSC.fitted_dose_response.v17.3.csv')
out_path = os.path.join(proj_dir, 'data/curated/celllines.csv')
lung_path = os.path.join(proj_dir, 'data/curated/celllines.lung.csv')

##################################    function    ######################################
def getBroadID(series):
    '''
    Get broad id of cell line.
    '''
    return [re.findall(r'\((ACH-.*)\)', sample)[0] for sample in series]


def getDrugNo(cell, drugtable):
    '''
    Get drug number of a cell line
    '''
    if cell in drugtable.Cell.values:
        return drugtable.loc[drugtable['Cell'] == cell,'Drug'].values[0]
    else:
        return 0



##################################    main    ######################################
### cell lines
cell_line = pd.read_csv(cellline_path, usecols=[0,1,2,5,6,7])
cell_line['Aliases'] = cell_line['CCLE_Name'].apply(lambda x: x.split('_')[0])

### gene expression
expr_ids = getBroadID(pd.read_csv(expr_path, sep='\t', header=0, skiprows=2, nrows=1).columns[2:].values)
cell_line['EXPR'] = cell_line['Broad_ID'].apply(lambda x: 1 if x in expr_ids else 0)

### mutation
mutation_ids = np.unique(pd.read_csv(mutation_path, sep='\t', usecols=['Broad_ID'], squeeze=True))
cell_line['MUTATION'] = cell_line['Broad_ID'].apply(lambda x: 1 if x in mutation_ids else 0)

### CNV
cnv_ids = pd.read_csv(cnv_path, usecols=[0]).iloc[:,0].tolist()
cell_line['CNV'] = cell_line['Broad_ID'].apply(lambda x: 1 if x in cnv_ids else 0)

## RPPA
rppa_ids = pd.read_csv(rppa_path, usecols=[0]).iloc[:,0].tolist()
cell_line['RPPA'] = cell_line['CCLE_Name'].apply(lambda x: 1 if x in rppa_ids else 0)

### RNAi
rnai_ids = pd.read_csv(rnai_path, nrows=1, header=0).columns[1:].values
cell_line['RNAi'] = cell_line['CCLE_Name'].apply(lambda x: 1 if x in rnai_ids else 0)

### CRISPRi
crispr_ids = pd.read_csv(crispr_path, usecols=[0]).iloc[:,0].tolist()
cell_line['CRISPRi'] = cell_line['Broad_ID'].apply(lambda x: 1 if x in crispr_ids else 0)

### GDSC
gdsc_data = pd.read_csv(gdsc_path, usecols=['CELL_LINE_NAME', 'DRUG_NAME'])
gdsc_data.columns = ['Cell', 'Drug']
gdsc_data['Cell'] = gdsc_data['Cell'].apply(util.cleanAlias)
gdsc_data = gdsc_data.groupby(by=['Cell'], as_index=False).count()
cell_line['GDSC'] = cell_line['Aliases'].apply(lambda x: getDrugNo(x, gdsc_data))

### sort
cell_line.sort_values(by=['Primary Disease', 'Aliases'], inplace=True)
cell_line = cell_line[['Broad_ID', 'Aliases', 'Primary Disease', 'EXPR', 'MUTATION', 'CNV', 'RPPA', 'RNAi', 'CRISPRi', 'GDSC']]
cell_line.to_csv(out_path, index=None)

### lung
cell_line_lung = cell_line.loc[cell_line['Primary Disease'].str.contains('Lung'),:]
cell_line_lung.to_csv(lung_path, index=None)


