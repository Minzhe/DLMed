####################################################################################
###                                expr.simu.py                                  ###
####################################################################################
proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import numpy as np
import pandas as pd
import pickle as pkl

############################     function     ##############################
def simulate_expr(expr, gene, resp):
    print('Preparing for gene {} ...'.format(gene))
    low_simu = np.percentile(expr[gene], 5, interpolation='lower')
    high_simu = np.percentile(expr[gene], 95, interpolation='higher')
    cells = expr.index[(expr[gene] < np.percentile(expr[gene], 33.3)) | (expr[gene] > np.percentile(expr[gene], 66.7))]
    df = resp.loc[resp['Cell'].isin(list(cells)),['Cell', 'Drug', 'LOG50']].rename(columns={'LOG50': 'Resp'})
    df['Gene'] = gene
    df['Expr'] = df['Cell'].apply(lambda x: expr.loc[x,gene])
    df = df[['Cell', 'Drug', 'Gene', 'Expr', 'Resp']]
    df['Low_expr'] = low_simu
    df['Low_resp'] = np.nan
    df['High_expr'] = high_simu
    df['High_resp'] = np.nan
    return df

def simulate_expr_sep(expr, gene, resp, out_dir):
    print('Preparing for gene {} ...'.format(gene))
    low_simu = np.percentile(expr[gene], 5, interpolation='lower')
    high_simu = np.percentile(expr[gene], 95, interpolation='higher')
    df = resp.loc[resp['Cell'].isin(list(expr.index)),['Cell', 'Drug', 'LOG50']].rename(columns={'LOG50': 'Resp'})
    df['Gene'] = gene
    df['Expr'] = df['Cell'].apply(lambda x: expr.loc[x,gene])
    df = df[['Cell', 'Drug', 'Gene', 'Expr', 'Resp']]
    df['Low_expr'] = low_simu
    df['Low_resp'] = np.nan
    df['High_expr'] = high_simu
    df['High_resp'] = np.nan
    out_path = os.path.join(out_dir, '{}.csv'.format(gene))
    df.to_csv(out_path, index=None)


############################       main      #################################
expr_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_RNAseq_cancergene.csv')
resp_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_drug.split.csv')
out_path = os.path.join(proj_dir, 'result/simulation/simulation.expression/simu_expr.driver_genes.csv')
out_dir = os.path.join(proj_dir, 'result/simulation/simulation.expression/simu_expr.gene')
genes = ['TP53', 'KRAS', 'KEAP1', 'STK11', 'EGFR', 'NF1', 'BRAF', 'SETD2', 'RBM10', 'MET', 'ARID1A', 'PIK3CA', 'SMARCA4', 'RB1', 'CDKN2A', 'U2AF1', 'ERBB2']

# driver
# expr = pd.read_csv(expr_path, index_col=0)[genes]
# resp = pd.read_csv(resp_path)
# simu = pd.concat([simulate_expr(expr, gene, resp) for gene in genes])
# simu.to_csv(out_path, index=None)

# all cancer gene
# expr = pd.read_csv(expr_path, index_col=0)
# genes = list(expr.columns)
# resp = pd.read_csv(resp_path)
# for gene in genes:
#     simulate_expr_sep(expr, gene, resp, out_dir)
out_dir = os.path.join(proj_dir, 'result/simulation/simulation.expression/simu_expr.gene_fill')
out_path = os.path.join(proj_dir, 'result/simulation/simulation.expression/simu_expr.all_genes.fill.csv')
simu = pd.concat([pd.read_csv(os.path.join(out_dir, path), usecols=['Drug', 'Gene', 'Low_resp', 'High_resp']) for path in os.listdir(out_dir)])
simu = simu.groupby(by=['Drug', 'Gene'], as_index=False).mean()
simu.to_csv(out_path, index=None)