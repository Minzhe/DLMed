####################################################################################################
###                                   combination.effect.py                                      ###
####################################################################################################
proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import pandas as pd
import numpy as np
import statsmodels.api as sm
import multiprocessing as mp

#####################################   function  ###################################
def find_combination(mut):
    '''
    Find combination of genes that have cells have no, single and double mutation.
    '''
    genes = list(set(np.unique(mut['Gene'])))
    mut_cell = {gene: set(mut['Cell'][mut['Gene'] == gene]) for gene in genes}
    # find combination
    comb_genes = list()
    for i in range(len(genes)):
        gene_i = genes[i]
        cell_i = mut_cell[gene_i]
        if len(cell_i) == 0:
            continue
        for j in range(i+1, len(genes)):
            gene_j = genes[j]
            cell_j = mut_cell[gene_j]
            if len(cell_j) == 0:
                continue
            tmp_comb = mut_matrix(cell_i, cell_j, set(mut['Cell']))
            if all(tmp_comb):
                tmp_comb = (gene_i, gene_j), tmp_comb
                comb_genes.append(sum(tmp_comb, ()))
    comb_genes = pd.DataFrame(comb_genes)
    comb_genes.columns = ['Gene1', 'Gene2', 'Mut0', 'Mut1', 'Mut2', 'Mut12', 'Cell0', 'Cell1', 'Cell2', 'Cell12']
    return comb_genes.sort_values(by=['Gene1', 'Gene2'])

def mut_matrix(set1, set2, total):
    no_mut = total - (set1 | set2)
    mut1 = set1 - (set1 & set2)
    mut2 = set2 - (set1 & set2)
    mut12 = set1 & set2
    return len(no_mut), len(mut1), len(mut2), len(mut12), no_mut, mut1, mut2, mut12


def cal_interaction_effect(comb, ic50, n=8):
    '''
    Calculate the double gene mutation interaction effect.
    '''
    comb_gene_eff = list()
    for i in range(n+1):
        length = comb.shape[0] // n
        s_idx, e_idx = i*length, min(i*length+length, comb.shape[0]+1)
        if (s_idx + 1) >= e_idx:
            continue
        tmp_comb = comb.iloc[s_idx:e_idx,:]
        print('{} - {}: '.format(s_idx, e_idx))
        manager = mp.Manager()
        return_dict = manager.dict()
        jobs = []
        for idx, row in tmp_comb.iterrows():
            gene1, gene2 = row['Gene1'], row['Gene2']
            print('.', end='')
            cell0, cell1, cell2, cell12 = row['Cell0'], row['Cell1'], row['Cell2'], row['Cell12']
            p = mp.Process(target=cal_pair_interaction, args=(gene1, gene2, cell0, cell1, cell2, cell12, ic50, return_dict,))
            jobs.append(p)
            p.start()
        for proc in jobs:
            proc.join()
        tmp_comb_eff = pd.concat(list(return_dict.values()), axis=0)
        comb_gene_eff.append(tmp_comb_eff)
        print('\n')

    comb_gene_eff = pd.concat(comb_gene_eff)

    return comb_gene_eff

def cal_pair_interaction(gene1, gene2, cell0, cell1, cell2, cell12, ic50, return_dict):
    '''
    Calculate gene pair interaction effect of all drugs
    '''
    inter_eff = list()
    for drug in list(ic50.columns):
        eff0 = get_ic50(ic50, drug, list(cell0), 'raw')
        eff1 = get_ic50(ic50, drug, list(cell1), 'raw')
        eff2 = get_ic50(ic50, drug, list(cell2), 'raw')
        eff12 = get_ic50(ic50, drug, list(cell12), 'raw')
        if all((len(eff0), len(eff1), len(eff2), len(eff12))):
            coef, pval = test_interaction(eff0, eff1, eff2, eff12)
            inter_eff.append([gene1, gene2, drug, len(eff0), len(eff1), len(eff2), len(eff12), coef, pval])
    if len(inter_eff) > 0:
        inter_eff = pd.DataFrame(inter_eff)
        inter_eff.columns = ['Gene1', 'Gene2', 'Drug', 'Mut0', 'Mut1', 'Mut2', 'Mut12', 'Coef', 'Pval']
        return_dict[(gene1, gene2)] = inter_eff

def get_ic50(ic50, drug, cells, value='raw'):
    val = ic50.loc[cells, drug]
    val = list(val[val.notnull()])
    if value == 'raw':
        return val
    if value == 'mean':
        return np.mean(val) if len(val) > 0 else np.nan
    if value == 'median':
        return np.median(val) if len(val) > 0 else np.nan

def test_interaction(y0, y1, y2, y12):
    '''
    Test whether there is an interaction effect between factor x1 and x2.
    '''
    y = np.array(list(y0) + list(y1) + list(y2) + list(y12))
    x1 = np.array([0]*len(y0) + [1]*len(y1) + [0]*len(y2) + [1]*len(y12))
    x2 = np.array([0]*len(y0) + [0]*len(y1) + [1]*len(y2) + [1]*len(y12))
    x12 = np.array([0]*len(y0) + [0]*len(y1) + [0]*len(y2) + [1]*len(y12))
    X = np.array([x1, x2, x12]).T
    ols = sm.OLS(y, X).fit()
    return round(ols.params[2], 6), round(-np.log(ols.pvalues[2]), 6)

def get_common_comb(comb1, comb2, suffix1, suffix2, how, mincell=1):
    comb1 = comb1.loc[(comb1['Mut0'] >= mincell) & (comb1['Mut1'] >= mincell) & (comb1['Mut2'] >= mincell) & (comb1['Mut12'] >= mincell),:]
    comb2 = comb2.loc[(comb2['Mut0'] >= mincell) & (comb2['Mut1'] >= mincell) & (comb2['Mut2'] >= mincell) & (comb2['Mut12'] >= mincell),:]
    return comb1.merge(comb2, left_on=['Gene1', 'Gene2', 'Drug'], right_on=['Gene1', 'Gene2', 'Drug'], how=how, suffixes=(suffix1, suffix2))




#######################################   main  ##########################################
# cell_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/cell_inference.csv')
# mut_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_Mutation_cancergene.csv')
# drug_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_drug.split.csv')
# out_path = os.path.join(proj_dir, 'result/combination_screen/inference.interaction_effect.csv')

# cells = pd.read_csv(cell_path, squeeze=True)
# mut = pd.read_csv(mut_path, usecols=['Cell', 'Gene'])
# mut = mut.loc[mut['Cell'].isin(cells),:]
# genes = list(set(np.unique(mut['Gene'])))
# ic50 = pd.read_csv(drug_path, usecols=['Cell', 'Drug', 'LOG50'])
# ic50 = ic50.pivot(index='Cell', columns='Drug', values='LOG50')
# ic50 = ic50.loc[ic50.index.to_series().isin(list(cells)),:]

# comb_genes = find_combination(mut)
# target_genes = ['TP53', 'KRAS', 'ALK', 'EGFR', 'BRAF', 'MET', 'PIK3CA']
# comb_genes = comb_genes.loc[comb_genes['Gene1'].isin(target_genes) | comb_genes['Gene2'].isin(target_genes),:]
# print(comb_genes.shape)
# comb_genes_drug = cal_interaction_effect(comb_genes, ic50, n=1)
# comb_genes_drug.to_csv(out_path, index=False)

diff_trainval_path = os.path.join(proj_dir, 'result/combination_screen/interaction_effect.train_val.csv')
diff_test_path = os.path.join(proj_dir, 'result/combination_screen/interaction_effect.inference.csv')
diff_trainval = pd.read_csv(diff_trainval_path)
diff_test = pd.read_csv(diff_test_path)
diff_merge = get_common_comb(diff_test, diff_trainval, suffix1='_infer', suffix2='_train', how='inner', mincell=1)
diff_merge = diff_merge.sort_values(by=['Coef_infer'])
diff_merge.to_csv(os.path.join(proj_dir, 'result/combination_screen/inference_interaction_effect.train_val.csv'), index=False)
print(diff_merge)