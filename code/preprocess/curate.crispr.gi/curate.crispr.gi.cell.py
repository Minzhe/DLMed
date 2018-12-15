###################################################################################
###                         curate.crispr.gi.cell.py                            ###
###################################################################################
# curate weissman's cell paper of dual crispr knockout

import numpy as np
import pandas as pd
import os
import re
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
import utility as util


proj_path = 'D:/projects/DLCell'
pheno_path = os.path.join(proj_path, 'data/CRISPRi/crispr.GI.cell.phenotype.txt')
out_dir = os.path.join(proj_path, 'data/curated')

##############################    function    #############################
def calPhenoScore(gene1, gene2, scoretable):
    '''
    Average phenotype score values for combination.
    '''
    scores = []
    print(gene1, gene2)
    idx = ((scoretable.iloc[:,0] == gene1) & (scoretable.iloc[:,1] == gene2)) | ((scoretable.iloc[:,1] == gene1) & (scoretable.iloc[:,0] == gene2))
    score = scoretable.loc[idx,].iloc[:,2]
    if len(score) != 0:
        return np.mean(scores)
    else:
        print('Warning: combination phenotype not found!')
        return np.nan


# def calsgRNACorr(pheno, col, fileout, figout):
#     '''
#     Calculate AB pair and BA pair score.
#     '''
#     scoretable['Gene1'] = scoretable.sgRNAs.apply(lambda x: x.split('++')[0])
#     scoretable['Gene2'] = scoretable.sgRNAs.apply(lambda x: x.split('++')[1])
#     scoretable = scoretable[['Gene1', 'Gene2', col]]
#     allsgs = scoretable.Gene1.append(scoretable.Gene2).unique()

#     pheno_sg = pd.DataFrame(list(combinations(allsgs, 2)), columns=['sgRNA1', 'sgRNA2'])
#     pheno_sg['score1'] = pheno_sg['score2'] = None
#     for idx, row in pheno_sg.iterrows():
#         print(idx)
#         sg1, sg2 = pheno_sg['sgRNA1'][idx], pheno_sg['sgRNA2'][idx]
#         pheno_sg.iloc[idx, 2] = scoretable.loc[(scoretable['Gene1'] == sg1) & (scoretable['Gene2'] == sg2),col].values
#         pheno_sg.iloc[idx, 3] = scoretable.loc[(scoretable['Gene1'] == sg2) & (scoretable['Gene2'] == sg1),col].values


##############################    main    #############################
# read table
pheno = pd.read_csv(pheno_path, sep='\t', skiprows=4, header=None, names=['sgRNAs', 'Jurkat_rep1', 'Jurkat_rep2', 'Jurkat_ave', 'K562_rep1', 'K562_rep2', 'K562_ave'])

# gene level
pheno['Gene1'] = pheno.sgRNAs.apply(lambda x: x.split('++')[0].split('_')[0])
pheno['Gene2'] = pheno.sgRNAs.apply(lambda x: x.split('++')[1].split('_')[0])
pheno = pheno[['Gene1', 'Gene2', 'Jurkat_rep1', 'Jurkat_rep2', 'Jurkat_ave', 'K562_rep1', 'K562_rep2', 'K562_ave']]

# ===================== K562 ======================= #
pheno_K562 = pd.DataFrame({'Gene1': pheno['Gene1'], 'Gene2': pheno['Gene2'], 'score': pheno['K562_ave']})
pheno_K562[['Gene1', 'Gene2']] = pheno[['Gene1', 'Gene2']].apply(np.sort, axis=1)
pheno_K562 = pheno_K562.loc[pheno_K562.score.notnull(),:]
pheno_K562 = pheno_K562.groupby(['Gene1', 'Gene2'], as_index=False).mean()

# compare single gene condition
# ---------- plot ---------- #
# plotNegative(phenotable=pheno_K562, path=os.path.join(proj_path, 'data/curated/plot/density.single.crispr.png'))

# filter out same gene combination
print('1')
pheno_K562 = util.rmSameTarget(pheno_K562)
print('2')
pheno_K562_norm = util.quantile_normalize(pheno_K562, center='negative')
print('3')
pheno_K562_norm.to_csv(os.path.join(out_dir, 'crispr.doubleKO.cell.geno.pheno.norm.dense.csv'), index=None)

# convert to trainging format
# geno_pheno_K562 = makeTrainingTable(pheno_K562)
# geno_pheno_K562_dense = makeTrainingTableDense(geno_pheno_K562)
# geno_pheno_K562.to_csv(os.path.join(out_path, 'crispr.doubleKO.geno.pheno.sparse.csv'), index=None)
# geno_pheno_K562_dense.to_csv(os.path.join(out_path, 'crispr.doubleKO.geno.pheno.dense.csv'), index=None)