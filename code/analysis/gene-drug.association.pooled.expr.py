###############################################################################################################
###                                   gene-drug.association.pooled.mut.py                                   ###
###############################################################################################################
proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import sys
import pandas as pd
import numpy as np
import pickle as pkl
from scipy.stats import ttest_ind
from statsmodels.distributions.empirical_distribution import ECDF
import multiprocessing as mp

###################################    function   ######################################
def pool_gene_drug_expr(expr, resp, out_dir):
    genes = list(expr.columns)
    drugs = np.unique(resp['Drug'])
    groups = divide_cells(expr)
    for i in range(0, len(genes), 16):
        tmp_genes = genes[i:i+16]
        jobs = []
        for gene in tmp_genes:
            print('Analyzing gene {} ...'.format(gene))
            p = mp.Process(target=pool_gene_drug_expr_onegene, args=(expr, gene, drugs, groups, out_dir,))
            jobs.append(p)
            p.start()
        for proc in jobs:
            proc.join()

def pool_gene_drug_expr_onegene(expr, gene, drugs, groups, out_dir):
    diff = list()
    low, high = groups[gene]
    for comp in drugs:
        tmp_resp = resp.loc[resp['Drug'] == comp,:]
        resp_low = tmp_resp.loc[tmp_resp['Cell'].isin(low), 'LOG50']
        resp_high = tmp_resp.loc[tmp_resp['Cell'].isin(high), 'LOG50']
        mean_diff = round(np.mean(resp_high) - np.mean(resp_low), 5)
        p_val = round(ttest_ind(resp_low, resp_high)[1], 5)
        diff.append((gene, comp, mean_diff, p_val))
    diff = pd.DataFrame(diff, columns=['GENE', 'COMP', 'OBSLOGIC50DIFF', 'PVAL'])
    out_path = os.path.join(out_dir, '{}.csv'.format(gene))
    diff.to_csv(out_path, index=None)

def divide_cells(expr):
    genes = list(expr.columns)
    groups = dict()
    for gene in genes:
        lower = np.percentile(expr[gene], 25)
        upper = np.percentile(expr[gene], 75)
        groups[gene] = list(expr.index[expr[gene] <= lower]), list(expr.index[expr[gene] >= upper])
    return groups

def read_simu_expr(path):
    data = pd.read_csv(path)
    data['LOGIC50DIFF'] = data['High_resp'] - data['Low_resp']
    data = data[['Gene', 'Drug', 'LOGIC50DIFF']].rename(columns={'Gene': 'GENE', 'Drug': 'COMP'})
    data = data.groupby(by=['GENE'], as_index=False).apply(lambda x: calculate_empirical_p(x, 'P.GENE'))
    data = data.groupby(by=['COMP'], as_index=False).apply(lambda x: calculate_empirical_p(x, 'P.COMP'))
    return data

def calculate_empirical_p(df, name):
    s = abs(df['LOGIC50DIFF'])
    ecdf = ECDF(s)
    df[name] = np.round(0-np.log(1-ecdf(s)), 4)
    return df


#####################################    main   #########################################
expr_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_RNAseq_cancergene.csv')
resp_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_drug.split.csv')
out_dir = os.path.join(proj_dir, 'result/simulation/single_expr_pooled_analysis/expr_simu_gene')
out_path = os.path.join(proj_dir, 'result/simulation/single_expr_pooled_analysis/single_expr_simu.gene_drug_association.csv')
# expr = pd.read_csv(expr_path, index_col=0)
# resp = pd.read_csv(resp_path, usecols=['Cell', 'Drug', 'LOG50'])

# pool_gene_drug_expr(expr, resp, out_dir)
# res = pd.concat([pd.read_csv(os.path.join(out_dir, path)) for path in os.listdir(out_dir)])
# res.to_csv(out_path, index=None)

# compare prediction with observation
res = pd.read_csv(out_path)
pred = read_simu_expr(os.path.join(proj_dir, 'result/simulation/simulation.expression/simu_expr.all_genes.fill.csv'))
pred = pred.merge(res, on=['GENE', 'COMP'])
out1_path = os.path.join(proj_dir, 'result/simulation/single_expr_pooled_analysis/single_expr_simu.gene_drug_association.pred_vs_obs.csv')
pred.to_csv(out1_path, index=False)
out2_path = os.path.join(proj_dir, 'result/simulation/single_expr_pooled_analysis/single_expr_simu.gene_drug_association.pred_vs_obs.top.csv')
pred = pred.loc[(pred['P.GENE'] > 3) & (pred['P.COMP'] > 3),:].sort_values(by=['P.GENE', 'P.COMP'], ascending=False).sort_values(by=['PVAL'], ascending=True)
pred = pred.loc[pred['LOGIC50DIFF'] * pred['OBSLOGIC50DIFF'] > 0,:]
pred.to_csv(out2_path, index=False)