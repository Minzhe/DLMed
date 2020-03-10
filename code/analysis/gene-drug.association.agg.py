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
def aggGene(df):
    data = df.copy()
    data.index = data.index.to_series().apply(lambda x: x.split('_')[0])
    data = data.groupby(by=data.index).agg(lambda x: max(x, key=abs))
    return data

def getTopN(df, target, cutoff):
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

def searchGene(gene, gene_list):
    res = [True for x in gene_list if x in gene]
    return len(res) > 0

def selectTarget(df, target):
    druggable_gene = ['ABL', 'AKT', 'ATM', 'ATR', 'BCL', 'BCL', 'BRCA', 'BRD', 'CASP', 'CDK', 'EML', 'ERBB', 'FAS', 'FGFR', 'FLT', 'HSP', 'IGF', 'JAK', 
                      'KIT', 'KRAS', 'MAP', 'MDM', 'MET', 'MTOR', 'MYC', 'PDGF', 'PIK', 'PPAR', 'PRK', 'SMO', 'SRC', 'STAT', 'TGF', 'TOP1', 'TP53', 'TRAF']
    topn = df.loc[df['Gene'].apply(lambda x: searchGene(x, druggable_gene)),:]
    topn = topn.loc[topn['Target'].notnull() & (topn['Target'] != 'Chemo'),:]

def plotHeatMap(df, path, cutoff=None):
    if cutoff is not None:
        df = df.loc[df.apply(lambda x: not all(x < cutoff), axis=1), df.apply(lambda x: not all(x < cutoff), axis=0)]
    f, ax = plt.subplots(figsize=(80,40))
    p.plotHeatmap(abs(df), cmap='YlGnBu', cluster=False, ax=ax)
    f.savefig(path)

def plotDensity(df, path):
    arr = df.values.reshape(-1)
    f, ax = plt.subplots(figsize=(12,8))
    p.plot_density(title='Distribution of gene drug effect', ax=ax, effect=arr)
    f.savefig(path)



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

# # get top n largest
# target = pd.read_csv(target_path, index_col=['Drug'], squeeze=True)
# drug_sens_change_sclc = agg_gene(drug_sens_change_sclc)
# topn = get_top_n(abs(drug_sens_change_sclc), target=target, cutoff=0.1)
# topn.to_csv(os.path.join(proj_dir, 'result/simulation/random_single_mut_simulation.gene-drug.association_median.high_sclc.csv'), index=False)
# drug_sens_change_nsclc = agg_gene(drug_sens_change_nsclc)
# topn = get_top_n(abs(drug_sens_change_nsclc), target=target, cutoff=0.1)
# topn.to_csv(os.path.join(proj_dir, 'result/simulation/random_single_mut_simulation.gene-drug.association_median.high_nsclc.csv'), index=False)

### single
model = 'gene_cnn'
diff_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis.{}'.format(model),     '{}.single_mut_simu.gene-drug.associ.median.csv'.format(model))
agg_path  = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis.{}'.format(model),     '{}.single_mut_simu.gene-drug.associ.median.gene.csv'.format(model))
top_path  = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis.{}'.format(model),     '{}.single_mut_simu.gene-drug.associ.median.gene_top.csv'.format(model))
heat_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis.{}/img'.format(model), '{}.single_mut_simu.gene-drug.associ.heat.png'.format(model))
dens_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis.{}/img'.format(model), '{}.single_mut_simu.gene-drug.associ.dist.png'.format(model))
target_path = os.path.join(proj_dir, 'data/curated/Lung/drug_target.csv')

drug_sens_change = pd.read_csv(diff_path, index_col=0)
drug_sens_change = aggGene(drug_sens_change)
drug_sens_change.to_csv(agg_path)

# cutoff = np.percentile(drug_sens_change.abs().values.reshape(-1), 95)
# print('Cut off of p-value 0.05 is: {}'.format(cutoff))
# plotHeatMap(drug_sens_change, heat_path)
# plotDensity(drug_sens_change, dens_path)

# target = pd.read_csv(target_path, index_col=['Drug'], squeeze=True)
# topn = getTopN(drug_sens_change, target, cutoff=cutoff)
# topn.to_csv(top_path, index=False)