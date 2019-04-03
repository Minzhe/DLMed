#######################################################################################
###                         curate.crispr.gi.nbt.3834.py                            ###
#######################################################################################
# curate nbt.3834 paper of dual crispr knockout

import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import quantile_transform, scale
import matplotlib.pyplot as plt
import seaborn as sns

proj_path = '/work/bioinformatics/s418336/projects/DLMed'
data_path = os.path.join(proj_path, 'data/CRISPR.GI/crispr.GI.nbt3834.phenotype.xlsx')

#################################    function    ###################################
def getSingleGenePheno(scoretable):
    '''
    Get single gene phentype score.
    '''
    gene_pheno = pd.DataFrame([tuple(row) for row in scoretable], columns=['Gene1', 'score'])
    gene_pheno = gene_pheno.groupby(['Gene1'], as_index=False).mean()
    gene_pheno['Gene2'] = 'negative'
    return gene_pheno

def getDualGenePheno(scoretable):
    '''
    Get dual gene phenotype score.
    '''
    idx = scoretable.apply(lambda x: x[0] == x[1], axis=1)
    return scoretable.loc[~idx,:]

def concatSingleDualGenePheno(singletable, dualtable):
    '''
    Concat single and dual gene knock out phenotype score.
    '''
    geno_pheno = singletable.append(dualtable, ignore_index=True)
    geno_pheno.loc[len(geno_pheno)] = ['negative', 'negative', 0.0]
    geno_pheno.sort_values(by=['Gene1', 'Gene2'], inplace=True)
    return geno_pheno



###########################################    phenotype    #############################################
# out_dir = os.path.join(proj_path, 'data/curated')
# out2_path = os.path.join(proj_path, 'data/curated/crispr.doubleKO.nbt3834.geno.pheno.norm.dense.csv')
# pheno = pd.read_excel(pheno_path, usecols=[0,1,2,3], names = ['Gene_Pair', 'Gene1_score', 'Gene2_score', 'score'])
# pheno['Gene1'] = pheno['Gene_Pair'].apply(lambda x: x.split('__')[0])
# pheno['Gene2'] = pheno['Gene_Pair'].apply(lambda x: x.split('__')[1])
# pheno = pheno[['Gene1', 'Gene2', 'Gene1_score', 'Gene2_score', 'score']]
# allgenes = np.sort(pheno.Gene1.append(pheno.Gene2).unique())

# single_gene_pheno = getSingleGenePheno(np.concatenate((pheno[['Gene1', 'Gene1_score']].values, pheno[['Gene2', 'Gene2_score']].values)))
# dual_gene_pheno = getDualGenePheno(pheno[['Gene1', 'Gene2', 'score']])
# geno_pheno = concatSingleDualGenePheno(single_gene_pheno, dual_gene_pheno)
# geno_pheno_norm = util.quantile_normalize(geno_pheno, center='negative')
# geno_pheno_norm.to_csv(os.path.join(out_dir, 'crispr.doubleKO.nbt3834.geno.pheno.norm.dense.csv'), index=None)

# geno_pheno_norm = concatSingleDualGenePheno(single_gene_pheno, dual_gene_pheno, quantile_norm=True)
# geno_pheno_norm.to_csv(os.path.join(out_dir, 'crispr.doubleKO.nbt3834.geno.pheno.norm.sparse.csv'), index=None)
# geno_pheno_norm = makeTrainingTable(geno_pheno_norm)
# geno_pheno_norm.to_csv(os.path.join(out_dir, 'crispr.doubleKO.nbt3834.geno.pheno.norm.dense.csv'), index=None)

###########################################    phenotype    #############################################
dense_out = os.path.join(proj_path, 'data/curated/K562/crispr.doubleKO.nbt3834.geno.giscore.dense.csv')
data = pd.read_excel(data_path, usecols=[0,6,7], names=['Gene_Pair', 'gi_t_score', 'gi_m_score'])
data['Gene1'] = data['Gene_Pair'].apply(lambda x: x.split('__')[0])
data['Gene2'] = data['Gene_Pair'].apply(lambda x: x.split('__')[1])
data = data[['Gene1', 'Gene2', 'gi_t_score', 'gi_m_score']]
data = data.loc[data.Gene1 != data.Gene2,:]
data.to_csv(dense_out, index=None)