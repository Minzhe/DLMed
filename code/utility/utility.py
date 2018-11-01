###################################################################################
###                               utility.py                                    ###
###################################################################################

import re
import numpy as np
import pandas as pd
from scipy.stats import shapiro, normaltest
from sklearn.preprocessing import quantile_transform, scale


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  normalize  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
def quantile_normalize(df, target='auto', ifcenter=True, ifscale=True):
    '''
    Quantile normalize data frame to target distribution and centered. 
    If set target to auto, then regular quantile normalization, otherwise use sklearn version of quantile normalization.
    '''
    scale(df, with_mean=ifcenter, with_std=ifscale, copy=False)
    # if target == 'auto':
    #     rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    #     df = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
        
    return df

def safe_scale(df, max_std=None):
    '''
    Scale and reduce outlier value.
    '''
    scale(df, with_mean=True, with_std=True, copy=False)
    if max_std is not None:
        df[df > max_std] = max_std
        df[df < -max_std] = -max_std
        scale(df, with_mean=True, with_std=True, copy=False)
    return df

def detect_unnormal(df, test, cutoff):
    '''
    Detect distribution that do not look like bell shape and unimodal.
    '''
    pval = []
    for col in df:
        if test == 'shapiro':
            stat, p = shapiro(df[col])
        elif test == 'normaltest':
            stat, p = normaltest(df[col])
        pval.append(p)
    pval = np.array(pval)
    return pval, list(df.columns[pval < cutoff])


# >>>>>>>>>>>>>>>>>>>>>>>>>>>  gene  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
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
    
    mut['SubGene'], mut['Loci'] = zip(*list(mut.apply(lambda x: gene_loc_encoder[(x['SubGene'], x['Loci'])], axis=1)))
    mut.sort_values(by=['Cell', 'SubGene', 'Loci'], inplace=True)
    
    return mut



# >>>>>>>>>>>>>>>>>>>>>>>>>>  cell line  <<<<<<<<<<<<<<<<<<<<<<<<<< #
def cleanCellLine(name):
    '''
    Clean cell line name.
    '''
    name = re.sub(r'_[ABCD]$', '', name.upper())
    name = re.sub(r'[-_\[\]]', '', name)
    if re.match(r'^H[0-9]+$', name):
        name = 'NCI' + name
    return name

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  drug  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
def cleanDrugName(name):
    '''
    Clean drug name.
    '''
    return re.sub(r'[-,\s]', '', name).upper()

def cleanDrugData(df, max_value=50, min_value=1e-5, duplicate='mean'):
    '''
    Clean drug sensitivity data.
    '''
    df.columns = ['Cell', 'Drug', 'LOG50']
    df.Drug = df.Drug.apply(cleanDrugName)
    df['LOG50'][df['LOG50'] > max_value] = max_value
    df['LOG50'][df['LOG50'] < min_value] = min_value
    df['LOG50'] = np.log(df['LOG50'])
    if duplicate == 'mean':
        df = df.groupby(by=['Cell', 'Drug'], as_index=False).mean()
    df.index = list(range(df.shape[0]))
    return df

