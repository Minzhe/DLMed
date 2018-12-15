####################################################################################################
###                               curate.cell.lines.encode.py                                    ###
####################################################################################################
# This script is to summarize the cell line information contained in different datasets.

import os
import re
import numpy as np
import pandas as pd

proj_dir = "D:/projects/DLCell"
metadata_path = os.path.join(proj_dir, 'data/ENCODE/metadata.tsv')
lungcell_path = os.path.join(proj_dir, 'data/curated/celllines.lung.csv')
out_path = os.path.join(proj_dir, 'data/curated/encode.celllines.lung.csv')


###################################     main     #################################
# metadata = pd.read_csv(metadata_path, sep='\t')
metadata = pd.read_table(metadata_path, sep='\t', usecols=['Assay', 'Biosample term name', 'Biosample type', 'Biosample organism', 'Experiment target', 'Library made from'])
metadata.columns = ['Assay', 'Sample', 'Type', 'Organism', 'Target', 'Library']

# human
metadata = metadata.loc[metadata['Organism'] == 'Homo sapiens',:]

# lung cell line
lungcell = metadata.loc[metadata['Type'] == 'cell line',:]
lungcell['Sample'] = lungcell['Sample'].apply(lambda x: re.sub(r'[-\.\s]', '', x).upper())
lung_ids = pd.read_csv(lungcell_path, usecols=['Aliases'], squeeze=True)
lungcell = lungcell.loc[lungcell['Sample'].isin(lung_ids),:]

# lung tissue
lungtissue = metadata.loc[metadata['Sample'].str.upper().str.contains('LUNG'),:]

# count unique and sort
lungdata = pd.concat([lungcell, lungtissue])
lungdata.drop(columns=['Organism', 'Library'], inplace=True)
lungdata.sort_values(by=['Sample', 'Assay', 'Target'], inplace=True)
lungdata = lungdata.loc[~lungdata.duplicated(),:]

# count
# lungdata['Num'] = 1
# lungdata.fillna('', inplace=True)
# lungdata = lungdata.groupby(by=['Assay', 'Sample', 'Target'], axis=0, sort=False, as_index=False).sum()

# output
# lungdata = lungdata[['Assay', 'Sample', 'Target']]
print(lungdata)
exit()
lungdata.to_csv(out_path, index=False)

# 
