#####################################################################################
###                             integrate.dummy.py                                ###
#####################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed/'
import os
import sys
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
sys.path.append(os.path.join(proj_dir, 'code'))
import utility.utility as util
import utility.plot as p
cmap4 = ["#95a5a6", "#3498db", "#2ecc71", "#e74c3c"]

#########################  main  #########################
file_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_drug.csv')
map1_path = os.path.join(proj_dir, 'data/curated/Lung/merged/cell.drug.matrix.cluster.ic50.sub.png')
map2_path = os.path.join(proj_dir, 'data/curated/Lung/merged/cell.drug.matrix.cluster.source.sub.png')
map3_path = os.path.join(proj_dir, 'data/curated/Lung/merged/cell.drug.matrix.source.sub.png')

drug_data = pd.read_csv(file_path)
drug_data = drug_data.loc[drug_data['Source'] != 'CTRP',:]

# >>>>>>>>> plot heatmap <<<<<<<<<<<< #
# IC50
max_value = drug_data.LOG50.max() + 1
drug_logic50 = drug_data[['Cell', 'Drug', 'LOG50']].pivot(index='Cell', columns='Drug').fillna(max_value)
drug_logic50.columns = drug_logic50.columns.droplevel(0)
drug_logic50
# source
drug_source = drug_data[['Cell', 'Drug', 'Source']]
source_dict = {'CCLE': 1, 'UTSW': 2, 'CTRP': 3}
drug_source.Source = drug_source.Source.apply(lambda x: source_dict[x])
drug_source = drug_source.pivot(index='Cell', columns='Drug').fillna(0)
drug_source.columns = drug_source.columns.droplevel(0)
#####  cluster by IC50  #####
f, ax1 = plt.subplots(figsize=(32,30))
drug_logic50 = p.plotHeatmap(drug_logic50, cmap='mako', mask_value=max_value, ax=ax1, cluster=True, return_df=True, cbar=False)
f.savefig(map1_path)
f, ax2 = plt.subplots(figsize=(32,30))
drug_source = drug_source.loc[drug_logic50.index, drug_logic50.columns]
p.plotHeatmap(drug_source, cmap=cmap4, mask_value=0, ax=ax2, cluster=False, return_df=False, cbar=False)
f.savefig(map2_path)
#####  cluster by source  #####
drug_source = util.clusterMatrix(drug_source, method='categorical', value_order=[1,2])
f, ax3 = plt.subplots(figsize=(32,30))
p.plotHeatmap(drug_source, cmap=cmap4, mask_value=0, ax=ax3, cluster=False, return_df=False, cbar=False)
f.savefig(map3_path)