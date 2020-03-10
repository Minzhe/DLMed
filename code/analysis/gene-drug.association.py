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
# from utility import utility as util
# import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('seaborn')
# from utility import plot as p

##################################  function  ######################################
def sign_max_df(df):
    s = df.apply(lambda x: sign_max(x), axis=0)
    return s

def sign_max(x):
    x = np.array(x)
    return x[np.argmax(abs(x))]

def aggMatrix(fun, mat):
    print('Aggregating matrix ...')
    if fun == 'median':
        return np.apply_along_axis(np.median, 0, mat)
    elif fun == 'mean':
        return np.apply_along_axis(np.mean, 0, mat)
    elif fun == 'max':
        return np.apply_along_axis(np.max, 0, mat)
    elif fun == '75.percentile':
        return np.apply_along_axis(lambda x: np.percentile(abs(x), 75), 0, mat)
    elif fun == 't.test.p':
        return np.apply_along_axis(lambda x: stats.ttest_1samp(x, 0)[1], 0, mat)
    else:
        raise ValueError('Unrecognized value aggregation.')

def readSubType(path):
    '''
    Read cell line meta information of lung cancer subtypes
    '''
    info = pd.read_csv(path)
    sclc = list(info['CCLE_Name'][info['Primary Disease'] == 'Lung SCLC'].apply(lambda x: x.split('_')[0]))
    nsclc = list(info['CCLE_Name'][info['Primary Disease'] == 'Lung NSCLC'].apply(lambda x: x.split('_')[0]))
    return sclc, nsclc

def readLociAnno(path):
    '''
    Read location annotation.
    '''
    # location annotation
    with open(loci_anno_path, 'rb') as f:
        data = pkl.load(f)
        loci, gene_index, drug_index, cell_index = data['alleles'], data['gene_index'], data['drug_index'], data['cell_index']
        loci2gene = [gene_index[gene] + '_' + str(loc) for gene, loc in loci]
    return gene_index, drug_index, cell_index, loci2gene

def readSimu(simu_paths, loci_anno_path, mut_path, subset, val='median'):
    '''
    Load simulation result, and aggregate multiple cell lines.
    '''
    # location annotation
    _, drug_index, cell_index, loci2gene = readLociAnno(loci_anno_path)
    mut_table = pd.read_csv(mut_path, index_col=0)
    # simulation
    diff_logic50 = list()
    for path in list(simu_paths):
        with open(path, 'rb') as f:
            data = pkl.load(f)
            perb_ic50, base_ic50 = data['res'][0], np.array(data['base'][0]['y_hat_me_cnn'])
            cell, drugs = cell_index[data['cellline_id']], [drug_index[id_] for id_ in data['drug_ids']]
            # check subset
            if subset is not None and cell not in subset: continue
            print('Loading {}'.format(path))
            tmp_diff = perb_ic50 - base_ic50
        # check original mutation status
        mut_idx = np.where(np.array(mut_table.loc[cell,:]) == 1)[0]
        tmp_diff[mut_idx,:] = tmp_diff[mut_idx,:] * (-1)
        diff_logic50.append(tmp_diff)
    # summarize, collapse to genes
    diff_logic50 = aggMatrix(fun='median', mat=np.array(diff_logic50))
    return pd.DataFrame(diff_logic50, index=loci2gene, columns=drugs)

def sparseMutationArrayLoci(mut_path, loci_anno_path, out_path):
    '''
    Summarize the mutation information for each cell line in loci level.
    '''
    # location annotation
    with open(loci_anno_path, 'rb') as f:
        data = pkl.load(f)
        loci, gene_index, drug_index, cell_index = data['alleles'], data['gene_index'], data['drug_index'], data['cell_index']
        loci2gene = [gene_index[gene]+'_'+str(loc) for gene, loc in loci]
    # mutation array
    with open(mut_path, 'rb') as f:
        mut = pkl.load(f)['mutation_array']
    cell_mut = np.zeros(shape=(len(cell_index), len(loci2gene)))
    for i in range(len(cell_index)):
        for j in range(len(loci)):
            if mut[i,loci[j][0],loci[j][1]] != 0:
                cell_mut[i,j] = 1
    cell_mut = pd.DataFrame(cell_mut, index=[cell_index[i] for i in range(len(cell_index))], columns=loci2gene)
    cell_mut.to_csv(out_path)

def SparseMutationArrayGene(mut_path, out_path):
    '''
    Summarize the mutation information for each cell line in gene level.
    '''
    mut_gene = pd.read_csv(mut_path, index_col=0).T
    mut_gene.index = mut_gene.index.to_series().apply(lambda x: x.split('_')[0])
    mut_gene = mut_gene.groupby(mut_gene.index).agg(lambda x: int(any(x))).T
    mut_gene.to_csv(out_path)


##################################  main  ######################################
model = 'gene_cnn'
geneset_path  = os.path.join(proj_dir, 'data/curated/geneset.pkl')
meta_path     = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version', 'cell.line.meta.csv')
mut_path      = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version', 'ccle_utsw.lung_MutAFExpr_cancergene_drugRes.pkl')
mut_loci_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version', 'ccle_utsw.lung_Mut_cancergene_loci_table.csv')
mut_gene_out  = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version', 'ccle_utsw.lung_Mut_cancergene_table.csv')
simu_path      = os.path.join(proj_dir, 'result/simulation/simulation.random_single_mut.{}'.format(model), 'ic50_cnn_model_cell_*_drug_*.pkl')
loci_anno_path = os.path.join(proj_dir, 'result/simulation/simulation.random_single_mut.{}'.format(model), 'ic50_cnn_model_alleles.pkl')

# simulated logic50 change
val = 'median'
simu_path = sorted(glob.glob(simu_path))
diff_logic50 = readSimu(simu_path, loci_anno_path, mut_loci_path, subset=None, val=val)
simu_out_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis.{}/{}.single_mut_simu.gene-drug.associ.{}.csv'.format(model, model, val))
diff_logic50.to_csv(simu_out_path)

