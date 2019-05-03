####################################################################################################################
###                                   combination.screen.kras.downstream.py                                      ###
####################################################################################################################
proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import sys
sys.path.append(os.path.join(proj_dir, 'code'))
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

###########################################    function   ##########################################
def plot_loci_importance(eff, path, rank='inverse_rank'):
    eff['Gene'] = eff.index.to_series().apply(lambda x: '_'.join(x.split('_')[:2]))
    eff['Loci'] = eff.index.to_series().apply(lambda x: int(x.split('_')[-1]))
    eff = eff.set_index(keys=['Gene', 'Loci']).abs()
    # get rank of importance within each gene
    if rank == 'rank':
        eff = eff.groupby(by=eff.index.get_level_values(0)).apply(lambda x: (x.rank()-1)/len(x))
    elif rank == 'inverse_rank':
        eff = eff.groupby(by=eff.index.get_level_values(0)).apply(lambda x: 1/x.rank(ascending=False))
    # average importance of each location across gene
    eff = eff.groupby(by=eff.index.get_level_values(1)).mean()
    plt.subplots(figsize=(48,24))
    sns.heatmap(eff, cmap='Blues')
    plt.savefig(path)

def plot_loci_gene(eff, loci, out_dir):
    path = os.path.join(out_dir, 'loci_importance.{}.png'.format(loci))
    eff['Loci'] = eff.index.to_series().apply(lambda x: int(x.split('_')[-1]))
    eff = eff.loc[eff['Loci'] == loci,:].drop(['Loci'], axis=1).mean(axis=1).reset_index()
    eff.columns = ['Gene', 'Eff']
    f, ax = plt.subplots(figsize=(45,12))
    ax = sns.barplot(x='Gene', y='Eff', data=eff)
    ax.set_xticklabels(eff['Gene'], rotation=90)
    f.savefig(path)


###########################################    main   ##########################################
model = 'full_cnn'
eff_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis.{}/single_mut_simulation.gene-drug.association.median.csv'.format(model))
img_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis.{}/img/loci_importance.rank.png'.format(model))
loci_dir = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis.{}/img'.format(model))
eff = pd.read_csv(eff_path, index_col=0)

plot_loci_importance(eff, img_path, rank='rank')
# plot_loci_gene(eff, 14, loci_dir)