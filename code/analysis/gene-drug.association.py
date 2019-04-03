####################################################################################################
###                                   gene-drug.association.py                                   ###
####################################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import sys
sys.path.append(os.path.join(proj_dir, 'code'))
import glob
import pandas as pd
import numpy as np
from scipy import stats
import pickle as pkl
from utility import utility as util
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('seaborn')
from utility import plot as p

##################################  function  ######################################
def sign_max_df(df):
    s = df.apply(lambda x: sign_max(x), axis=0)
    return s

def sign_max(x):
    x = np.array(x)
    return x[np.argmax(abs(x))]

def read_simu(paths, loci2gene, cell_index, drug_index, mut_table, subset, val='median'):
    '''
    Load simulation result
    '''
    diff_logic50 = list()
    for path in list(paths):
        with open(path, 'rb') as f:
            data = pkl.load(f)
            perb_ic50, base_ic50 = data['res'][0], np.array(data['base'][0]['y_hat_me_cnn'])
            cell, drugs = cell_index[data['cellline_id']], [drug_index[id_] for id_ in data['drug_ids']]
            # check subtype
            if subset is not None and cell not in subset:
                continue
            print('Loading {}'.format(path))
            tmp_diff = perb_ic50 - base_ic50
        # check original mutation status
        mut_idx = np.where(np.array(mut_table.loc[cell,:]) == 1)[0]
        tmp_diff[mut_idx,:] = tmp_diff[mut_idx,:] * (-1)
        diff_logic50.append(tmp_diff)
    # summarize
    diff_logic50 = np.array(diff_logic50)
    print(diff_logic50.shape)
    if val == 'median':
        diff_logic50 = np.apply_along_axis(np.median, 0, diff_logic50)
    if val == 'mean':
        diff_logic50 = np.apply_along_axis(np.mean, 0, diff_logic50)
    if val == 'max':
        diff_logic50 = np.apply_along_axis(np.max, 0, diff_logic50)
    if val == '75percentile':
        diff_logic50 = np.apply_along_axis(lambda x: np.percentile(abs(x), 75), 0, diff_logic50)
    if val == 't.test.p':
        diff_logic50 = np.apply_along_axis(lambda x: stats.ttest_1samp(x, 0)[1], 0, diff_logic50)
    # collpase to genes
    print(diff_logic50.shape)
    diff_logic50 = pd.DataFrame(diff_logic50, index=loci2gene, columns=drugs)
    return diff_logic50

##################################  main  ######################################
geneset_path = os.path.join(proj_dir, 'data/curated/geneset.pkl')
meta_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/cell.line.meta.csv')
simu_path = glob.glob(os.path.join(proj_dir, 'result/simulation/random_single_mut_simulation/ic50_cnn_model_cell_*_drug_*.pkl'))
loci_anno_path = os.path.join(proj_dir, 'result/simulation/random_single_mut_simulation/ic50_cnn_model_alleles.pkl')
mut_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_MutAFExpr_cancergene_drugRes.pkl')
mut_out = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_Mut_cancergene_loci_table.csv')

# location annotation
with open(loci_anno_path, 'rb') as f:
    data = pkl.load(f)
    loci, gene_index, drug_index, cell_index = data['alleles'], data['gene_index'], data['drug_index'], data['cell_index']
    loci2gene = [gene_index[gene]+'_'+str(loc) for gene, loc in loci]

# cell line mutation table
# with open(mut_path, 'rb') as f:
#     mut = pkl.load(f)['mutation_array']
# cell_mut = np.zeros(shape=(len(cell_index), len(loci2gene)))
# for i in range(len(cell_index)):
#     print(i)
#     for j in range(len(loci)):
#         if mut[i,loci[j][0],loci[j][1]] != 0:
#             cell_mut[i,j] = 1
# cell_mut = pd.DataFrame(cell_mut, index=[cell_index[i] for i in range(len(cell_index))], columns=loci2gene)
# cell_mut.to_csv(mut_out)
cell_mut = pd.read_csv(mut_out, index_col=0)

# cell line info
info = pd.read_csv(meta_path)
sclc = list(info['CCLE_Name'][info['Primary Disease'] == 'Lung SCLC'].apply(lambda x: x.split('_')[0]))
nsclc = list(info['CCLE_Name'][info['Primary Disease'] == 'Lung NSCLC'].apply(lambda x: x.split('_')[0]))

# simulated logic50 change
val = '75percentile'
diff_logic50 = read_simu(simu_path, loci2gene=loci2gene, cell_index=cell_index, drug_index=drug_index, mut_table=cell_mut, subset=None, val=val)
simu_out_path = os.path.join(proj_dir, 'result/simulation/random_single_mut_simulation.gene-drug.association_{}.abs.csv'.format(val))
diff_logic50.to_csv(simu_out_path)

