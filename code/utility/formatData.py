######################################################################################
###                               formatData.py                                    ###
######################################################################################

import re
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.preprocessing import OneHotEncoder

############################   drug sensitivity   ############################
def buildTrainingArray(mut_array, drug_array):
    '''
    Build SVM traning data from mutation and drug sentivity data.
    '''
    # drug
    drug_onehot = OneHotEncoder(sparse=False)
    drugs = drug_onehot.fit_transform(np.reshape(drug_array[:,1], (-1,1)))
    # cell
    cells = np.reshape(mut_array, (mut_array.shape[0],-1))
    # merge
    X = np.concatenate([np.array(list(map(lambda x: cells[int(x)], drug_array[:,0]))), drugs], axis=1)
    y = drug_array[:,2]
    # remove 0 column
    idx = np.apply_along_axis(sum, 0, X) == 0
    X = X[:,~idx]
    return  csr_matrix(X), y


############################   fitness   ############################
def createTrainingTable(scoretable, duplicate='remove'):
    '''
    Convert double knock out score table to training table suitable for machine learning model. 
    Input columns: Gene1, Gene2, score
    '''
    scoretable, duptable = rmDuplicate(scoretable, duplicate=duplicate)
    allgenes = scoretable.Gene1.append(scoretable.Gene2).unique()
    allgenes = allgenes[allgenes != 'negative']
    geno_pheno = pd.DataFrame(0, index=scoretable.index, columns=np.concatenate([allgenes, ['score']]))
    for idx, row in scoretable.iterrows():
        print(idx)
        geno_pheno.loc[idx, 'score'] = row['score']
        if row['Gene1'] != 'negative':
            geno_pheno.loc[idx, row['Gene1']] = 1
        if row['Gene2'] != 'negative':
            geno_pheno.loc[idx, row['Gene2']] = 1
    return geno_pheno, duptable


def sparse2dense(sparsetable):
    '''
    Convert sparse one hot encoded table to dense indexed table.
    '''
    def getGenes(idx, genes):
        targets = tuple(genes[idx == 1])
        if len(targets) == 2:
            return targets
        elif len(targets) == 1:
            return targets + (None,)
        elif len(targets) == 0:
            return (None, None)
        else:
            raise ValueError('Target gene numbers not correct: {}.'.format(targets))
    genes = sparsetable.columns[:-1].to_series()
    targets = []
    for idx, row in sparsetable.iterrows():
        targets.append(getGenes(row[:-1], genes))
    denseTable = pd.DataFrame(targets, columns=['Gene1', 'Gene2'])
    denseTable['score'] = sparsetable['score'].tolist()
    return denseTable


def rmDuplicate(scoretable, duplicate='remove'):
    '''
    Check if table contains duplicates.
    '''
    scoretable[['Gene1', 'Gene2']] = scoretable[['Gene1', 'Gene2']].apply(np.sort, axis=1)
    scoretable.sort_values(by=['Gene1', 'Gene2'], inplace=True)
    scoretable.reset_index(drop=True, inplace=True)
    idx = scoretable.duplicated(subset=['Gene1', 'Gene2'], keep=False)
    if sum(idx) != 0:
        print('Warning: duplicated row found.')
        print(scoretable.loc[idx,:])
        if duplicate == 'remove':
            return scoretable.loc[~idx,:], scoretable.loc[idx,:]
        elif duplicate == 'average':
            scoretable_ave = scoretable.groupby(by=['Gene1', 'Gene2'], sort=False, as_index=False).mean()
            return scoretable_ave, scoretable.loc[idx,:]
    else:
        print('No duplicates found.')
        return scoretable, None


def rmSameTarget(phenotable):
    '''
    Remove entry with two same genes.
    '''
    same_idx = (phenotable.Gene1 == phenotable.Gene2) & (phenotable.Gene1 != 'negative')
    return phenotable.loc[~same_idx,:]


def rmNegative(scoretable):
    '''
    Remove dual negative sample.
    '''
    idx = (scoretable['Gene1'] == 'negative') & (scoretable['Gene2'] == 'negative')
    return scoretable.loc[~idx,:]