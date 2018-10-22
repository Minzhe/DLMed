###################################################################################
###                               utility.py                                    ###
###################################################################################

import re
import numpy as np
import pandas as pd
from sklearn.preprocessing import quantile_transform, scale


# =========================  normalize  ============================= #
def quantile_normalize(scoretable, center):
    '''
    quantile normalize to normal distribution and center to negative value
    '''
    scoretable.score = quantile_transform(scoretable.score.values.reshape(-1,1), output_distribution='normal', random_state=1234)
    scoretable.score = scale(scoretable.score.values, with_mean=False, with_std=True)
    neg = scoretable.loc[(scoretable.iloc[:,0] == 'negative') & (scoretable.iloc[:,1] == 'negative'),'score'].values
    if center == 'negative':
        scoretable.score = scoretable.score - neg
    elif center == '0':
        scoretable.score = scoretable.score - scoretable.score.mean()
    return scoretable


# =========================  gene  ============================= #
def cleanAlias(name):
    '''
    Get alias of cell line name.
    '''
    alias = re.sub(r'\(.*\)', '', name)
    alias = re.sub(r'[-_\[\]]', '', alias)
    return alias.upper()


class geneEncoder(object):
    '''
    Convert between different cell line naming system.
    '''
    def __init__(self, broadid, alias):
        self.broad_dict = {ach: name for (ach, name) in zip(broadid, alias)}
        # self.alias_dict = {name: ach for (ach, name) in zip(broadid, alias)}
    
    def broad2name(self, broadids):
        return [self.broad_dict.get(ach, 'NA') for ach in broadids]


def geneLocEncoder(genes, locs, maxwin=96):
    '''
    Encode gene locations.
    @return: dict[gene, startsite] = subgene, win
    '''
    mut_loc = pd.DataFrame({'gene': genes, 'loc':locs})
    mut_loc.sort_values(by=['gene', 'loc'], inplace=True)
    mut_loc.drop_duplicates(inplace=True)
    mut_loc['win'] = 0

    # add window index
    genes = np.unique(mut_loc.gene)
    for gene in genes:
        mut_loc.loc[mut_loc['gene'] == gene,'win'] = list(range(1, sum(mut_loc['gene'] == gene)+1))
    
    # subdivide window
    max_sub = max(mut_loc.win) // maxwin + 1
    gene_loc_encoder = dict()
    for i in range(1, max_sub+1):
        tmp_mat = mut_loc.loc[(mut_loc.win <= i*maxwin) & (mut_loc.win > (i-1)*maxwin),:]
        # create mapping dict
        for idx, row in tmp_mat.iterrows():
            tmp_gene = row['gene'] + '_' + str(i)
            tmp_win = row['win'] % maxwin
            gene_loc_encoder[(row['gene'], row['loc'])] = (tmp_gene, tmp_win)

    return gene_loc_encoder


def divideGene(mut, gene_loc_encoder):
    '''
    Divide gene with windows.
    '''
    if sum(pd.Series(['Cell', 'Gene', 'MutStart']).isin(mut.columns.tolist())) != 3:
        raise ValueError('Parameter mut should contains columns "Cell", "Gene", "MutStart".')
    mut = mut[['Cell', 'Gene', 'MutStart']]
    mut.sort_values(by=['Cell', 'Gene', 'MutStart'], inplace=True)
    mut.reset_index(drop=True, inplace=True)
    mut.columns = ['Cell', 'SubGene', 'Loci']
    
    for idx, row in mut.iterrows():
        if idx % 100 == 1: print(idx)
        mut.loc[idx,['SubGene','Loci']] = gene_loc_encoder[(row['SubGene'], row['Loci'])]
    mut.sort_values(by=['Cell', 'SubGene', 'Loci'], inplace=True)
    
    return mut



# =========================  cell line  ============================= #
def cleanCellLine(name):
    '''
    Clean cell line name.
    '''
    name = re.sub(r'[-_\[\]]', '', name.upper())
    if re.match(r'^H[0-9]+$', name):
        name = 'NCI' + name
    return name


