##################################################################################################
###                               curate.cell.lines.utsw.py                                    ###
##################################################################################################
# This script is to summarize the cell line information contained in different datasets.

import os
import re
import numpy as np
import pandas as pd

proj_dir = "D:/projects/DLCell"
cellline_path = os.path.join(proj_dir, 'data/UTSW_John/Cell.Line.Database.v1.72.csv')
ccle_path = os.path.join(proj_dir, 'data/curated/celllines.lung.csv')
expr_path = os.path.join(proj_dir, 'data/UTSW_John/UTSW Cell Lines.csv')
mutation_path = os.path.join(proj_dir, 'data/UTSW_John/All Mutations (CPRIT and COSMIC-84).csv')
cnv_path = os.path.join(proj_dir, 'data/UTSW_John/All CNVs by Genes.csv')
methy_path = os.path.join(proj_dir, 'data/UTSW_John/450K Methylation Profiling (CLD) 2-1-12.csv')
drug_path = os.path.join(proj_dir, 'data/UTSW_John/table6_chemical_screen_data.csv')
out_path = os.path.join(proj_dir, 'data/curated/utsw.john.cell.lines.csv')


##############################    function     ################################
def cleanCellLine(cell):
    '''
    Function to clean cell line name.
    '''
    return re.sub(r'[-\.\s]', '', cell).upper()


##############################    main     ################################
cell_line = pd.read_csv(cellline_path, index_col=0, usecols=['Sort', 'Cell Line', 'Tumor Type'])
cell_line.columns = ['Cell', 'Type']
cell_line.Cell = cell_line.Cell.apply(cleanCellLine)
cell_line.sort_values(by=['Type', 'Cell'], inplace=True)

### lung
cell_line = cell_line.loc[cell_line['Type'].str.contains('Lung', na=False),:]

### ccle
ccle_lung = pd.read_csv(ccle_path, usecols=['Aliases'], squeeze=True).tolist()
cell_line['CCLE'] = cell_line['Cell'].apply(lambda x: 1 if x in ccle_lung else 0)

### expr
expr_cell = pd.read_csv(expr_path, nrows=1).columns[33:].to_series().apply(cleanCellLine).tolist()
cell_line['EXPR'] = cell_line['Cell'].apply(lambda x: 1 if x in expr_cell else 0)

### mutation
mutation_cell = pd.read_csv(mutation_path, nrows=1).columns[3:].to_series().apply(cleanCellLine).tolist()
cell_line['Mutation'] = cell_line['Cell'].apply(lambda x: 1 if x in mutation_cell else 0)

### cnv
cnv_cell = pd.read_csv(cnv_path, nrows=1).columns[5:].to_series().apply(cleanCellLine).tolist()
cnv_cell = set([cell.split('_')[0] for cell in cnv_cell])
cell_line['CNV'] = cell_line['Cell'].apply(lambda x: 1 if x in cnv_cell else 0)

### methylation
methy_cell = pd.read_csv(methy_path, nrows=1).columns[7:-14].to_series().apply(cleanCellLine).tolist()
cell_line['METHY'] = cell_line['Cell'].apply(lambda x: 1 if x in methy_cell else 0)

### drug
drug_cell = pd.read_csv(drug_path, nrows=1).columns[7:].to_series().apply(cleanCellLine).tolist()
cell_line['DRUG'] = cell_line['Cell'].apply(lambda x: 1 if x in drug_cell else 0)
print(cell_line.iloc[:,2:].apply(sum, axis=0))

cell_line.to_csv(out_path, index=None)