####################################################################################################
###                                   gene-drug.association.py                                   ###
####################################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import sys
sys.path.append(os.path.join(proj_dir, 'code'))
import glob
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.distributions.empirical_distribution import ECDF
import pickle as pkl
from utility import utility as util
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('seaborn')
from utility import plot as p

###############################   function  ##################################
def agg_gene(df):
    data = df.copy()
    data.index = data.index.to_series().apply(lambda x: x.split('_')[0])
    data = data.groupby(by=data.index).agg(lambda x: max(x, key=abs))
    return data

def get_top_n(df, target, cutoff):
    data = df.copy()
    data.index.name = 'Gene'
    data.reset_index(inplace=True)
    # data frame
    data = pd.melt(data, id_vars=['Gene'], value_vars=data.columns[1:], var_name='Drug', value_name='LogFC')
    data = data.loc[data['LogFC'].sort_values(ascending=False).index,:]
    ecdf = ECDF(abs(data['LogFC']))
    data['p_val'] = 1 - ecdf(abs(data['LogFC']))
    # filter
    data = data.loc[(abs(data['LogFC']) > cutoff) & (data['Drug'].apply(lambda x: not x.startswith('SW'))),:]
    data['Target'] = data['Drug'].apply(lambda x: target[x])
    data = data[['Gene', 'Drug', 'Target', 'LogFC', 'p_val']]
    data = data.loc[abs(data['LogFC']).sort_values(ascending=False).index,:]
    return data

def search_gene(gene, gene_list):
    res = [True for x in gene_list if x in gene]
    return len(res) > 0


###############################   main  ##################################
# ### sclc & nsclc
# diff1_path = os.path.join(proj_dir, 'result/simulation/random_single_mut_simulation.gene-drug.association_median.NSCLC.csv')
# diff2_path = os.path.join(proj_dir, 'result/simulation/random_single_mut_simulation.gene-drug.association_median.SCLC.csv')
# fig_path = os.path.join(proj_dir, 'result/simulation/fig/gene-drug.association.heatmap.median.subtype_diff.png')
# target_path = os.path.join(proj_dir, 'data/curated/Lung/drug_target.csv')

# drug_sens_change_sclc = pd.read_csv(diff1_path, index_col=0)
# drug_sens_change_nsclc = pd.read_csv(diff2_path, index_col=0)
# drug_sens_change = abs(drug_sens_change_nsclc - drug_sens_change_sclc)
# # aggregate
# drug_sens_change.index = drug_sens_change.index.to_series().apply(lambda x: x.split('_')[0])
# drug_sens_change = drug_sens_change.groupby(by=drug_sens_change.index).agg(lambda x: max(x, key=abs))

# # # plot heatmap
# # drug_sens_change = drug_sens_change.loc[drug_sens_change.apply(lambda x: not all(x < 0.1), axis=1), drug_sens_change.apply(lambda x: not all(x < 0.1), axis=0)]
# # f, ax = plt.subplots(figsize=(80,40))
# # p.plotHeatmap(abs(drug_sens_change), cmap='YlGnBu', cluster=True, ax=ax)
# # f.savefig(fig_path, transparent=True)
# # exit()

# # get top n largest
# target = pd.read_csv(target_path, index_col=['Drug'], squeeze=True)
# drug_sens_change_sclc = agg_gene(drug_sens_change_sclc)
# topn = get_top_n(abs(drug_sens_change_sclc), target=target, cutoff=0.1)
# topn.to_csv(os.path.join(proj_dir, 'result/simulation/random_single_mut_simulation.gene-drug.association_median.high_sclc.csv'), index=False)
# drug_sens_change_nsclc = agg_gene(drug_sens_change_nsclc)
# topn = get_top_n(abs(drug_sens_change_nsclc), target=target, cutoff=0.1)
# topn.to_csv(os.path.join(proj_dir, 'result/simulation/random_single_mut_simulation.gene-drug.association_median.high_nsclc.csv'), index=False)

### single
diff_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis/random_single_mut_simulation.gene-drug.association_median.csv')
target_path = os.path.join(proj_dir, 'data/curated/Lung/drug_target.csv')
# shuffle data
drug_sens_change = pd.read_csv(diff_path, index_col=0)
# n_row, n_col = drug_sens_change.shape
# shuffle_change = np.random.choice(drug_sens_change.values.reshape(-1), n_row*n_col, replace=False).reshape(n_row, n_col)
# shuffle_change = pd.DataFrame(shuffle_change, index=drug_sens_change.index, columns=drug_sens_change.columns)
drug_sens_change = agg_gene(drug_sens_change)
# shuffle_change = agg_gene(shuffle_change)
cutoff = np.percentile(abs(drug_sens_change).values.reshape(-1), 95)
print('Cut off of p-value 0.05 is: {}'.format(cutoff))

# traget
target = pd.read_csv(target_path, index_col=['Drug'], squeeze=True)
topn = get_top_n(drug_sens_change, target, cutoff=0.1)

druggable_gene = ['ABL', 'AKT', 'ATM', 'ATR', 'BCL', 'BCL', 'BRCA', 'BRD', 'CASP', 'CDK', 'EML', 'ERBB', 'FAS', 'FGFR', 'FLT', 'HSP', 'IGF',\
                  'JAK', 'KIT', 'KRAS', 'MAP', 'MDM', 'MET', 'MTOR', 'MYC', 'PDGF', 'PIK', 'PPAR', 'PRK', 'SMO', 'SRC', 'STAT', 'TGF', 'TOP1',\
                  'TP53', 'TRAF']
topn = topn.loc[topn['Gene'].apply(lambda x: search_gene(x, druggable_gene)),:]
topn = topn.loc[topn['Target'].notnull() & (topn['Target'] != 'Chemo'),:]
topn.to_csv(os.path.join(proj_dir, 'result/simulation/random_single_mut_simulation.gene-drug.association_median.high.csv'), index=False)