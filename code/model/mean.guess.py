#############################################################################
###                            mean.guess.py                              ###
#############################################################################

import os
import pandas as pd
import pickle as pkl
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
plt.style.use('seaborn')
import seaborn as sns

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
train_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_MutExprCNV_cancergene_drug_loci.untruncated.label_12_70_20_10_35324_train.csv')
val_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_MutExprCNV_cancergene_drug_loci.untruncated.label_12_70_20_10_9120_val.csv')
test_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_MutExprCNV_cancergene_drug_loci.untruncated.label_12_70_20_10_5608_inference.csv')
mean_path = os.path.join(proj_dir, 'result/genomic.drugsens/mean.guess.png')
median_path = os.path.join(proj_dir, 'result/genomic.drugsens/median.guess.png')

train_data = pd.read_csv(train_path, index_col=0)
val_data = pd.read_csv(val_path, index_col=0)
test_data = pd.read_csv(test_path, index_col=0)

train_ave = train_data.groupby(by=['drug_id'])['score'].mean()
train_median = train_data.groupby(by=['drug_id'])['score'].median()

# predict
val_data['mean_pred'] = val_data['drug_id'].apply(lambda x: train_ave[x])
val_data['median_pred'] = val_data['drug_id'].apply(lambda x: train_median[x])
test_data['mean_pred'] = test_data['drug_id'].apply(lambda x: train_ave[x])
test_data['median_pred'] = test_data['drug_id'].apply(lambda x: train_median[x])

r2_val_mean = r2_score(val_data['score'], val_data['mean_pred'])
r2_val_median = r2_score(val_data['score'], val_data['median_pred'])
r2_test_mean = r2_score(test_data['score'], test_data['mean_pred'])
r2_test_median = r2_score(test_data['score'], test_data['median_pred'])

f, (ax1, ax2) = plt.subplots(1, 2, figsize=(18,8))
ax1 = sns.regplot(val_data['mean_pred'], val_data['score'], scatter_kws={'alpha':0.5}, ax=ax1)
ax2 = sns.regplot(test_data['mean_pred'], test_data['score'], scatter_kws={'alpha':0.5}, ax=ax2)
ax1.set_title('Mean guess (val)')
ax2.set_title('Mean guess (test)')
f.savefig(mean_path)

f, (ax1, ax2) = plt.subplots(1, 2, figsize=(18,8))
ax1 = sns.regplot(val_data['median_pred'], val_data['score'], scatter_kws={'alpha':0.5}, ax=ax1)
ax2 = sns.regplot(test_data['median_pred'], test_data['score'], scatter_kws={'alpha':0.5}, ax=ax2)
ax1.set_title('Median guess (val)')
ax2.set_title('Median guess (test)')
f.savefig(median_path)


