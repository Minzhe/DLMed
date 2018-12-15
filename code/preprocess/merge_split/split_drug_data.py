#################################################################################
###                               split_data.py                               ###
#################################################################################
# hold 10% of cell line in UTSW for test

proj_dir = '/work/bioinformatics/s418336/projects/DLMed/'
import os, sys
import numpy as np
import pandas as pd
sys.path.append(os.path.join(proj_dir, 'code'))
import matplotlib.pyplot as plt
import utility.plot as p
import utility.utility as util


#########################        main        #########################
out_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_drug.split.csv')
fig1_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_drug.split_point.png')
fig2_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_drug.split_cell.png')
# all cell lines
train_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_MutExprCNV_cancergene_drug_loci.untruncated.label_12_70_20_10_35324_train.csv')
val_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_MutExprCNV_cancergene_drug_loci.untruncated.label_12_70_20_10_9120_val.csv')
test_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_MutExprCNV_cancergene_drug_loci.untruncated.label_12_70_20_10_5608_inference.csv')
train_cell = set(pd.read_csv(train_path, usecols=['cell_id'], squeeze=True))
val_cell = set(pd.read_csv(val_path, usecols=['cell_id'], squeeze=True))
test_cell = set(pd.read_csv(test_path, usecols=['cell_id'], squeeze=True))
train_drug = set(pd.read_csv(train_path, usecols=['drug_id'], squeeze=True))
val_drug = set(pd.read_csv(val_path, usecols=['drug_id'], squeeze=True))
test_drug = set(pd.read_csv(test_path, usecols=['drug_id'], squeeze=True))
all_cell = train_cell | val_cell | test_cell
all_drug = train_drug | val_drug | test_drug

# cell line source
ccle_mut_path = os.path.join(proj_dir, 'data/curated/Lung/ccle/ccle.lung_Mutation_cancergene.csv')
utsw_mut_path = os.path.join(proj_dir, 'data/curated/Lung/utsw.mw/utsw.lung_Mutation_cancergene.csv')
ccle_cell = set(pd.read_csv(ccle_mut_path, usecols=['Cell'], squeeze=True))
utsw_cell = set(pd.read_csv(utsw_mut_path, usecols=['Cell'], squeeze=True))
# utsw_only_cell = utsw_cell - ccle_cell
# print(utsw_only_cell)

# read drug ic50 data
drug_path = os.path.join(proj_dir, 'data/curated/Lung/merged/ccle_utsw.lung_drug.csv')
drug = pd.read_csv(drug_path)
drug['Dropped_previously'] = ~(drug.Cell.isin(list(all_cell)) & drug.Drug.isin(list(all_drug)))

# split on data points
inference_idx = drug.Cell.isin(list(test_cell))
drug['Train_split_point'] = util.split_train_val_test(drug.Drug, train_prop=0.7, val_prop=0.2, masked=(inference_idx | drug.Dropped_previously))
drug.loc[inference_idx, 'Train_split_point'] = 'Inference'
# drug = drug.loc[~drug.Dropped_previously,:]
print(drug.Train_split_point.value_counts())
# train: 31099
# val: 8855
# test: 4490
# inference: 5662

# split on cell line
drug['Train_split_cell'] = np.nan
drug.loc[drug.Cell.isin(list(train_cell)),'Train_split_cell'] = 'Train'
drug.loc[drug.Cell.isin(list(val_cell)),'Train_split_cell'] = 'Val'
drug.loc[drug.Dropped_previously,'Train_split_cell'] = np.nan
drug.loc[drug.Cell.isin(list(test_cell)),'Train_split_cell'] = 'Inference'
print(drug.Train_split_cell.value_counts())
# train: 35324
# val: 9120
# inference: 5662

drug.to_csv(out_path, index=None)

# plot heatmap
heat_mat = drug.pivot(index='Cell', columns='Drug', values='Train_split_point').fillna(0)
heat_mat.replace('Train', 5, inplace=True)
heat_mat.replace('Val', 4, inplace=True)
heat_mat.replace('Test', 3, inplace=True)
heat_mat.replace('Inference', 2, inplace=True)
f, ax = plt.subplots(figsize=(32,30))
p.plotHeatmap(heat_mat, cmap='mako', mask_value=0, ax=ax, cbar=False)
f.savefig(fig1_path)


heat_mat = drug.pivot(index='Cell', columns='Drug', values='Train_split_cell').dropna(how='all', axis=0).dropna(how='all', axis=1).fillna(0)
heat_mat.replace('Train', 1, inplace=True)
heat_mat.replace('Val', 2, inplace=True)
heat_mat.replace('Inference', 3, inplace=True)
f, ax = plt.subplots(figsize=(32,30))
p.plotHeatmap(heat_mat, cmap=['#F2F2F6', '#FC8D62', '#CB1B4F', '#401B44'], mask_value=0, ax=ax, cbar=False)
f.savefig(fig2_path)