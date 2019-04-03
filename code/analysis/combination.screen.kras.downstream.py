####################################################################################################################
###                                   combination.screen.kras.downstream.py                                      ###
####################################################################################################################
proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import sys
sys.path.append(os.path.join(proj_dir, 'code'))
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns
import utility.plot as p

################################    main   ################################
model = 'full_cnn'
eff_path = os.path.join(proj_dir, 'result/simulation/simulation.combination_screen_kras.{}.analysis/kras.combination.eff.csv'.format(model))
anno_path = os.path.join(proj_dir, 'result/simulation/simulation.combination_screen_kras.{}.analysis/ic50_cnn_model_alleles.pkl'.format(model))
out_path = os.path.join(proj_dir, 'result/simulation/simulation.combination_screen_kras.{}.analysis/kras.combination.top.anno.{}.csv'.format(model, model))
fig_path = os.path.join(proj_dir, 'result/simulation/simulation.combination_screen_kras.{}.analysis/kras.combination.heatmap.1.{}.png'.format(model, model))
target_path = os.path.join(proj_dir, 'data/curated/Lung/drug_target.csv')

eff = pd.read_csv(eff_path)
f = open(anno_path, 'rb')
anno = pkl.load(f)
gene_index, drug_index = anno['gene_index'], anno['drug_index']
eff.columns = ['drug', 'gene', 'coef', 'neg.log.p']
target = pd.read_csv(target_path, index_col=0, squeeze=True)
# annotation
eff['drug'] = eff['drug'].apply(lambda x: drug_index[x])
eff['gene'] = eff['gene'].apply(lambda x: gene_index[x])
eff['target'] = eff['drug'].apply(lambda x: target[x])
eff = eff.loc[abs(eff['coef']).sort_values(ascending=False).index,:]
# cutoff = np.percentile(abs(eff['coef']), 95)
# eff = eff.loc[abs(eff['coef']) > cutoff,:]
eff = eff[['drug', 'target', 'gene', 'coef', 'neg.log.p']]
# eff.to_csv(out_path, index=False)

##################   heatmap  ####################
eff = eff.pivot(index='drug', columns='gene', values='coef')
vmax = eff.abs().values.max()
f, ax = plt.subplots(figsize=(60,35))
p.plotHeatmap(eff, cmap='RdBu_r', ax=ax, cluster=True, vmin=0-vmax, vmax=vmax)
# sns.clustermap(eff, vmin=0-vmax, vmax=vmax, cmap='RdBu_r')
f.savefig(fig_path)
