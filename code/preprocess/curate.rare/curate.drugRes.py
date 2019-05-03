##############################################################################
###                          curate.drugRes.py                             ###
##############################################################################
proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import numpy as np
import pandas as pd
import pickle as pkl
import sys
sys.path.append(os.path.join(proj_dir, 'code'))
import utility.utility as util
import utility.plot as p
import matplotlib.pyplot as plt
import seaborn as sns

##############################   function   ########################
def plot_heatmap(df, path):
    data = df.pivot(index='CELL', columns='DRUG', values='LOGIC50')
    vmax = np.nanmax(data.abs())
    plt.subplots(figsize=(20, 16))
    p.plotHeatmap(data, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    plt.savefig(path)

def label_tissue(cell_id, rare, lung):
    tissue = ''
    if cell_id in rare: tissue += '|Rare'
    if cell_id in lung: tissue += '|Lung'
    return tissue.strip('|')

def indexArray(df):
    cells = df[['CELL', 'TISSUE']].drop_duplicates().set_index('CELL').sort_index().squeeze()
    cell_index = dict()
    for i, cell in enumerate(cells.index):
        cell_index[i] = (cell, cells[cell])
        df['CELL'][df['CELL'] == cell] = i
    drug_index = dict()
    for i, drug in enumerate(np.sort(df['DRUG'].unique())):
        drug_index[i] = drug
        df['DRUG'][df['DRUG'] == drug] = i
    return cell_index, drug_index, np.array(df[['CELL', 'DRUG', 'LOGIC50']])


###############################    main   ###############################
rare_path = os.path.join(proj_dir, 'data/curated/Rare/cellline.rare.csv')
ccle_path = os.path.join(proj_dir, 'data/DepMap/DepMap.celllines.2018.q3.csv')
utsw_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_drug.split.csv')
gdsc_path = os.path.join(proj_dir, 'data/COSMIC/GDSC.fitted_dose_response.v17.3.csv')
out_path = os.path.join(proj_dir, 'data/curated/Rare/rare_lung.drugRes.csv')
img_path = os.path.join(proj_dir, 'data/curated/Rare/img.png')
pkl_path = os.path.join(proj_dir, 'data/curated/Rare/rare_lung.drugRes.one_hot.pkl')

# # get lung and rare cell line
# rare = pd.read_csv(rare_path, usecols=['COSMIC_ID'], squeeze=True)
# ccle = pd.read_csv(ccle_path, usecols=['CCLE_Name', 'COSMIC_ID'])
# lung = ccle.loc[ccle['CCLE_Name'].apply(lambda x: True if 'LUNG' in x else False),'COSMIC_ID'].dropna().apply(int)
# cell = list(rare) + list(lung)

# # get drug ic50
# gdsc = pd.read_csv(gdsc_path, usecols=['COSMIC_ID', 'CELL_LINE_NAME', 'DRUG_ID', 'DRUG_NAME', 'PUTATIVE_TARGET', 'LN_IC50'])
# gdsc = gdsc.loc[gdsc['COSMIC_ID'].isin(cell),:]
# gdsc['Tissue'] = gdsc['COSMIC_ID'].apply(lambda x: label_tissue(x, list(rare), list(lung)))
# gdsc['CELL_LINE_NAME'] = gdsc['CELL_LINE_NAME'].apply(util.cleanCellLine)
# gdsc['DRUG_NAME'] = gdsc[['DRUG_ID', 'DRUG_NAME']].apply(lambda x: util.cleanDrugName(x[1]) + '.' + str(x[0]), axis=1)
# gdsc = gdsc.drop(['COSMIC_ID', 'DRUG_ID'], axis=1)
# gdsc.columns = ['CELL', 'DRUG', 'TARGET', 'LOGIC50', 'TISSUE']
# tissue = gdsc[['CELL', 'TISSUE']].drop_duplicates().set_index('CELL').sort_index().squeeze()

# # get utsw lung drug
# utsw = pd.read_csv(utsw_path, usecols=[0,1,2,3])
# utsw = utsw.loc[utsw['Source'] == 'UTSW',:].drop(['Source'], axis=1)
# utsw.columns = ['CELL', 'DRUG', 'LOGIC50']
# utsw['TARGET'] = np.nan
# utsw['TISSUE'] = utsw['CELL'].apply(lambda x: tissue[x] if x in tissue.index else 'Lung')

# # merge
# data = pd.concat([gdsc, utsw])
# data.to_csv(out_path, index=None)

# index table
data = pd.read_csv(out_path)
gene_index, drug_index, data = indexArray(data)
with open(pkl_path, 'wb') as f:
    pkl.dump({'gene_index':gene_index, 'drug_index':drug_index, 'logic50': data}, file=f)
# plot_heatmap(data, img_path)
