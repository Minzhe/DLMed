####################################################################################################
###                                   gene-drug.association.py                                   ###
####################################################################################################
proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import pandas as pd
import numpy as np
from scipy.stats import norm

######################################      function     ###################################
def readLociEff(full_path, gene_path, actual_path, gene_list=None, p_eff_cut=None, p_gene_cut=None, p_drug_cut=None):
    # read simulation result
    res_full, full_cut = read_loci_eff(full_path, p_eff_cut, p_gene_cut, p_drug_cut)
    res_gene, gene_cut = read_loci_eff(gene_path, p_eff_cut, p_gene_cut, p_drug_cut)
    data = res_full.merge(res_gene, on=['Gene', 'Drug'], suffixes=['-full.model', '-gene.model'])
    # subset on gene list
    genes = data['Gene'].apply(lambda x: x.split('_')[0])
    if gene_list is not None:
        data = data.loc[genes.isin(gene_list),:]
    # subset on cut off
    if full_cut is not None and gene_cut is not None:
        data = data.loc[(data['LogIC50Diff-full.model'].abs() > full_cut) & (data['LogIC50Diff-gene.model'].abs() > gene_cut),:]
    # read actual value
    data = append_obs(data, actual_path)
    return data

### >>>>>>>>>>>>>>>>>>>    function    <<<<<<<<<<<<<<<<<<< ###
def read_loci_eff(path, p_eff_cut=None, p_gene_cut=None, p_drug_cut=None):
    data = pd.read_csv(path, index_col=0)
    # general cut off
    cutoff = calculate_percentile(data, p_eff_cut) if p_eff_cut is not None else None
    # calculate p value based on z score
    data = meltSubset(data, p_gene_cut=p_gene_cut, p_drug_cut=p_drug_cut)
    return data, cutoff

def meltSubset(df, p_gene_cut=None, p_drug_cut=None):
    p_gene = df.apply(calculate_p, axis=1)
    p_drug = df.apply(calculate_p, axis=0)
    # melt data frame
    data = pd.melt(df.reset_index(), id_vars=['index'], var_name='Drug', value_name='LogIC50Diff')
    p_gene = pd.melt(p_gene.reset_index(), id_vars=['index'], var_name='Drug', value_name='neg.log.p.gene')
    p_drug = pd.melt(p_drug.reset_index(), id_vars=['index'], var_name='Drug', value_name='neg.log.p.drug')
    # merge
    data = data.merge(p_gene, on=['index', 'Drug'])
    data = data.merge(p_drug, on=['index', 'Drug']).rename(columns={'index': 'Gene'})#.sort_values(by=['Gene', 'p.val'])
    # subset
    data['neg.log.p.gene'] = calculate_neg_log_p(data['neg.log.p.gene'])
    data['neg.log.p.drug'] = calculate_neg_log_p(data['neg.log.p.drug'])
    if p_gene_cut is not None:
        data = data.loc[data['neg.log.p.gene'] > calculate_neg_log_p(p_gene_cut),:]
    if p_drug_cut is not None:
        data = data.loc[data['neg.log.p.drug'] > calculate_neg_log_p(p_drug_cut),:]
    data.index = list(range(data.shape[0]))
    return data

def append_obs(data, obs_path):
    print('Appending actual observation.')
    # append freq
    obs = pd.read_csv(obs_path, index_col=0, usecols=[0,1], squeeze=True)
    freq = dict()
    for gene in obs.keys(): freq[gene] = obs[gene]
    data['Freq'] = data.groupby(by=['Gene'], as_index=True, sort=False)['Gene'].apply(lambda x: append_obs_freq(x, freq))
    # append observation
    obs = pd.read_csv(obs_path, index_col=0).drop(['Freq'], axis=1)
    data = append_obs_eff(data, obs)
    return data

def calculate_percentile(data, p):
    cutoff = np.percentile(data.abs().values.reshape(-1), 100-p*100)
    print('Cut off for effective interaction: {}, pval {}'.format(cutoff, p))
    return cutoff

def calculate_p(x):
    z = (x - np.mean(x)) / np.std(x)
    p = norm.cdf(-abs(z))*2
    return p

def calculate_neg_log_p(x, keep=4):
    return round(-np.log(x), keep)

def get_max_eff(df):
    return df.loc[df['LogIC50Diff'].abs().idxmax(),:]

def append_obs_freq(s, freq):
    count = freq.get(list(set(s))[0], 0)
    return pd.Series([count] * len(s), index=s.index)

def append_obs_eff(data, obs):
    obs.index.name = 'Gene'
    obs = obs.reset_index()
    obs = obs.melt(id_vars=['Gene'], var_name='Drug', value_name='LogIC50Diff-obs')
    data = data.merge(obs, on=['Gene', 'Drug'], how='left')
    return data


#######################################     main    ######################################
# gene_list = ['TP53', 'KRAS', 'KEAP1', 'STK11', 'EGFR', 'NF1', 'BRAF', 'SETD2', 'RBM10', 'MET', 'ARID1A', 'PIK3CA', 'SMARCA4', 'RB1', 'CDKN2A', 'U2AF1', 'ERBB2']
gene_list = None
# data_path   = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis.{}/{}.single_mut_simu.gene-drug.associ.median.csv')
# full_path   = data_path.format('full_cnn', 'full_cnn')
# gene_path   = data_path.format('gene_cnn', 'gene_cnn')
# actual_loci_path = os.path.join(proj_dir, 'result/simulation/single_mut_pooled_analysis/single_mut_pooled.gene-drug.association.actual.csv')
# out_loci_path    = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis/single_mut_simu.gene-drug.associ.median.loci.pval.csv')
# data = readLociEff(full_path, gene_path, actual_loci_path, gene_list=gene_list, p_eff_cut=None, p_gene_cut=None, p_drug_cut=None)
# data.to_csv(out_loci_path, index=None)

data_path   = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis.{}/{}.single_mut_simu.gene-drug.associ.median.gene.csv')
full_path   = data_path.format('full_cnn', 'full_cnn')
gene_path   = data_path.format('gene_cnn', 'gene_cnn')
actual_gene_path = os.path.join(proj_dir, 'result/simulation/single_mut_pooled_analysis/single_mut_pooled.gene-drug.association.actual.gene.csv')
out_gene_path    = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis/single_mut_simu.gene-drug.associ.median.gene.pval.csv')
data = readLociEff(full_path, gene_path, actual_gene_path, gene_list=gene_list, p_eff_cut=None, p_gene_cut=None, p_drug_cut=None)
data.to_csv(out_gene_path, index=None)