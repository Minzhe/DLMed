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
import multiprocessing as mp
from statsmodels.distributions.empirical_distribution import ECDF
# pd.set_option('display.max_rows', None)
import matplotlib.pyplot as plt
import seaborn as sns

###############################   function  ##################################
### >>>>>>>>>>>>>>>>>>>>>    aggregation    <<<<<<<<<<<<<<<<<<<<<<< ###
def agg_gene(df, agg):
    data = df.copy()
    data.index = data.index.to_series().apply(lambda x: x.split('_')[0])
    if agg == 'max':
        data = data.groupby(by=data.index).agg(lambda x: max(x, key=abs))
    elif agg == 'median':
        data = data.groupby(by=data.index).agg(np.median)
    else:
        raise ValueError('Unrecognized agg function.')
    return data

def calculate_p(df, drug):
    sw_df = df[drug]
    df = abs(df)
    # empirical distribution
    p_val = list()
    for idx, row in df.iterrows():
        ecdf = ECDF(list(row))
        p = 1 - ecdf(abs(sw_df.loc[idx,:]))
        p_val.append(p)
    p_val = pd.DataFrame(p_val, index=sw_df.index, columns=[item+'_p' for item in drug])
    sw_df = round(pd.concat([sw_df, p_val], axis=1).sort_index(axis=1), 6)
    return sw_df

def summarize_drug(drugs, out_dir, **data):
    for drug in drugs:
        out_path = os.path.join(out_dir, 'gene-drug.association.{}.csv'.format(drug))
        drug_df = []
        for name in data.keys():
            df = data[name]
            colnames = df.columns.tolist()
            colnames = [col for col in colnames if drug in col]
            df = df[colnames]
            df.columns = [col+'_'+name for col in colnames]
            drug_df.append(df)
        drug_df = pd.concat(drug_df, axis=1)
        drug_df = drug_df.loc[drug_df.iloc[:,1] < 0.05,:]
        drug_df = drug_df.loc[drug_df.iloc[:,1].sort_values().index,:]
        drug_df.to_csv(out_path)

### >>>>>>>>>>>>>>>>>>>>>    heatmap    <<<<<<<<<<<<<<<<<<<<<<< ### 
def plot_simu_heatmap(diff_path, mut_loci_path, mut_gene_path, gene, out_dir, drug):
    def loci(x): return (x[0]-1)*96+x[1]
    fig_path = os.path.join(out_dir, '{}.All.png'.format(gene))
    # location and gene level mut frequency
    mut_loci = pd.read_csv(mut_loci_path, index_col=0).sum()
    mut_gene = pd.read_csv(mut_gene_path, index_col=0).sum()
    gene_freq = mut_gene[gene]
    # simu drug response
    diff = pd.read_csv(diff_path, index_col=0)
    diff = diff.loc[diff.index.to_series().apply(lambda x: x.split('_')[0]) == gene,:]
    diff['Freq'] = mut_loci[diff.index]
    diff['Loci'] = diff.index.to_series().apply(lambda x: loci(list(map(int, x.split('_')[1:]))))
    diff = diff.sort_values(by=['Loci'])
    diff.index = diff.apply(lambda x: '{}.{}'.format(int(x['Loci']), int(x['Freq'])), axis=1)
    diff = diff.drop(['Freq', 'Loci'], axis=1)
    # heatmap
    plt.subplots(figsize=(60,24))
    sns.heatmap(diff, vmin=-diff.abs().values.max(), vmax=diff.abs().values.max(), cmap='RdBu_r')
    plt.title('Gene {} mutation response map (# of cell line with mutation: {})'.format(gene, gene_freq))
    plt.savefig(fig_path)


###############################   main  ##################################
##############   gene drug interaction with three SW drugs  ################
'''
diff_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis/random_single_mut_simulation.gene-drug.association_median.csv')
# out_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis/random_single_mut_simulation.gene-drug.association_median.SW_drug.csv')
out_dir = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis')
drug = ['SW036310', 'SW022906', 'SW069087']

diff = pd.read_csv(diff_path, index_col=0)
diff_max = agg_gene(diff, agg='max')
diff_median = agg_gene(diff, agg='median')
diff_max = calculate_p(diff_max, drug)
diff_median = calculate_p(diff_median, drug)
summarize_drug(drug, out_dir, max=diff_max, median=diff_median)
'''

##############   gene drug interaction with three SW drugs  ################
diff_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis/random_single_mut_simulation.gene-drug.association_median.csv')
mut_loci_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_Mut_cancergene_loci_table.csv')
mut_gene_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_Mut_cancergene_table.csv')
out_dir = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis/img/')
drug = ['SW036310', 'SW022906', 'SW069087', 'SW197335']
plot_simu_heatmap(diff_path, mut_loci_path, mut_gene_path, 'RBM10', out_dir, drug)
# plot_simu_heatmap(diff_path, mut_loci_path, mut_gene_path, 'VHL', out_dir, drug)
# plot_simu_heatmap(diff_path, mut_loci_path, mut_gene_path, 'HIF1A', out_dir, drug)
# plot_simu_heatmap(diff_path, mut_loci_path, mut_gene_path, 'KRAS', out_dir, drug)
exit()

"""
Three drugs:
TP53_1: 673, TP53_1: 674
TP53_1_69: 7578208
TP53_1_32:
TP53_1_5:

SW036310
VHL_1: 694
VHL_1_3: 

SW069087
HIF1A_1: 296
HIF1A_1_8: 62207796
HIF1A_1_5: 62207705
"""

# mut_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_Mutation_cancergene.csv')
# mut = pd.read_csv(mut_path)
# mut = mut.loc[mut.Gene == 'HIF1A',:].sort_values(by=['MutStart'])
# print(np.unique(mut['MutStart']))

# sw022906_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis_SW/SW022906.csv')
# sw036310_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis_SW/SW036310.csv')
# sw069087_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis_SW/SW069087.csv')
# sw022906 = pd.read_csv(sw022906_path, index_col=0).loc['TP53_1_69',]
# sw036310 = pd.read_csv(sw036310_path, index_col=0).loc['TP53_1_69',]
# sw069087 = pd.read_csv(sw069087_path, index_col=0).loc['TP53_1_69',]
# print(sw022906.sort_values().tail(8))
# print(sw036310.sort_values().tail(8))
# print(sw069087.sort_values().tail(8))


##############   raw cell line data  ################
'''
file_path = os.path.join(proj_dir, 'result/simulation/simulation.random_single_mut/single_mutation_diff_logic50/*.csv')
file_path = glob.glob(file_path)
sw036310, sw022906, sw069087 = dict(), dict(), dict()
for path in file_path:
    print(path)
    cell = path.split('/')[-1].split('_')[-1].strip('.csv')
    data = pd.read_csv(path, index_col=0)   
    sw036310[cell] = data['SW036310']
    sw022906[cell] = data['SW022906']
    sw069087[cell] = data['SW069087']
sw036310 = pd.DataFrame(sw036310, index=data.index)
sw022906 = pd.DataFrame(sw022906, index=data.index)
sw069087 = pd.DataFrame(sw069087, index=data.index)
sw036310.to_csv(os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis_SW/SW036310.csv'))
sw022906.to_csv(os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis_SW/SW022906.csv'))
sw069087.to_csv(os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis_SW/SW069087.csv'))
'''
