###################################################################################
###                               utility.py                                    ###
###################################################################################

import re
import numpy as np
import pandas as pd
from scipy.spatial import distance as dist
from scipy.cluster import hierarchy as h
from scipy.stats import shapiro, normaltest
from sklearn.preprocessing import quantile_transform, scale
import matplotlib.pyplot as plt


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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  split data  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
def split_train_val_test(index, train_prop=0.7, val_prop=0.2, masked=None):
    labels = index.copy()
    if masked is not None:
        if len(masked) != len(index):
            raise ValueError('maksed should have same length as index, or None.')
        else:
            labels[masked.astype(bool)] = np.nan
    items = list(set(labels[labels.notnull()]))
    for item in items:
        train_idx, val_idx, test_idx = split_train_val_test_once(np.where(labels == item)[0], train_prop, val_prop)
        labels[train_idx] = 'Train'
        labels[val_idx] = 'Val'
        labels[test_idx] = 'Test'
    return labels

def split_train_val_test_once(sequence, train, val):
    seq = sequence.copy()
    np.random.shuffle(seq)
    train_idx = seq[:int(round(len(seq)*train))]
    val_idx = seq[int(round(len(seq)*train)):int(round(len(seq)*(train+val)))]
    test_idx = seq[int(round(len(seq)*(train+val))):]
    if len(train_idx) == 0 or len(val_idx) == 0 or len(test_idx) == 0:
        raise ValueError('Not enough long sequence to split')
    return train_idx, val_idx, test_idx

def split_train_val(index, train_prop=0.75, masked=None):
    if masked is not None:
        if len(masked) != len(index):
            raise ValueError('masked should have same length as index, or None.')
        else:
            items = np.array(list(set(index[~masked])))
    else:
        items = np.array(list(set(index)))
    np.random.shuffle(items)
    train = items[:int(round(len(items)*train_prop))]
    val = items[int(round(len(items)*train_prop)):]
    train_idx = index.isin(list(train))
    val_idx = index.isin(list(val))

    labels = index.copy()
    labels[train_idx] = 'Train'
    labels[val_idx] = 'Val'
    labels[~(train_idx | val_idx)] = np.nan
    return labels

# >>>>>>>>>>>>>>>>>>>>>>>>>  data matrix  <<<<<<<<<<<<<<<<<<<<<<<<<< #
def sort_by_index(df, index, ascending, axis=0):
    if axis == 0:
        df['index'] = index
        df.sort_values(by=['index'], ascending=ascending, axis=0, inplace=True)
        del df['index']
        return df
    elif axis == 1:
        df.loc['index'] = index
        df.sort_values(by=['index'], ascending=ascending, axis=1, inplace=True)
        df.drop(index=['index'], inplace=True)
    return df

def clusterMatrix(df, method, value_order=None):
    '''
    Reordering matrix by hierarchical clustering.
    '''
    if method == 'hclust':
        # col
        corr = np.asarray(df.corr())
        col_linkage = h.linkage(dist.pdist(corr), method='average')
        col_dn = np.array(h.dendrogram(col_linkage, distance_sort='ascending', color_threshold=0.01)['ivl'], dtype='int16')
        # row
        corr = np.asarray(df.T.corr())
        row_linkage = h.linkage(dist.pdist(corr), method='average')
        row_dn = np.array(h.dendrogram(row_linkage, distance_sort='ascending', color_threshold=0.01)['ivl'], dtype='int16')
        df = df.iloc[row_dn, col_dn]
    elif method == 'categorical' and value_order is not None:
        df.sort_index(inplace=True)
        df = df.reindex(sorted(df.columns), axis='columns')
        # sort row
        for value in list(reversed(value_order)):
            index = df.apply(lambda x: 1 if value in x.values else 0, axis=1)
            df = sort_by_index(df, index=index, ascending=False, axis=0)
        # sort col
        for value in list(reversed(value_order)):
            index = df.apply(lambda x: 1 if value in x.values else 0, axis=0)
            df = sort_by_index(df, index=index, ascending=False, axis=1)
    return df

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
    if sum(pd.Series(['Cell', 'Gene', 'MutStart', 'AF']).isin(mut.columns.tolist())) != 4:
        raise ValueError('Parameter mut should contains columns "Cell", "Gene", "MutStart" and "AF".')
    mut = mut[['Cell', 'Gene', 'MutStart', 'AF']]
    mut.sort_values(by=['Cell', 'Gene', 'MutStart'], inplace=True)
    mut.reset_index(drop=True, inplace=True)
    mut.columns = ['Cell', 'SubGene', 'Loci', 'AF']
    
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

def cleanDrugData(df, max_value=None, min_value=None, duplicate='mean'):
    '''
    Clean drug sensitivity data.
    '''
    df.columns = ['Cell', 'Drug', 'LOG50']
    df.Drug = df.Drug.apply(cleanDrugName)
    if max_value is not None:
        df = df.loc[df['LOG50'] < max_value,:]
    if min_value is not None:
        df = df.loc[df['LOG50'] > min_value,:]
    df['LOG50'] = np.log(df['LOG50'])
    if duplicate == 'mean':
        df = df.groupby(by=['Cell', 'Drug'], as_index=False).mean()
    df.index = list(range(df.shape[0]))
    return df

