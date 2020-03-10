####################################################################################################
###                                         read_simu.py                                         ###
####################################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import sys
sys.path.append(os.path.join(proj_dir, 'code'))
import glob
import pandas as pd
import numpy as np
import pickle as pkl
import multiprocessing as mp
from utility import utility as util


##################################  function  ######################################
def read_simu(path, loci2gene, cell_index, drug_index):
    '''
    Load simulation result
    '''
    print('Loading {}'.format(path))
    with open(path, 'rb') as f:
        data = pkl.load(f)
        perb_ic50, base_ic50 = data['res'][0], np.array(data['base'][0]['y_hat_me_cnn'])
        cell, drugs = cell_index[data['cellline_id']], [drug_index[id_] for id_ in data['drug_ids']]
    diff_logic50 = round(pd.DataFrame(perb_ic50-base_ic50, index=loci2gene, columns=drugs), 5)
    out_path = os.path.join(os.path.dirname(path), 'single_mutation_diff_logic50', 'diff_logic50_{}.csv'.format(cell))
    diff_logic50.to_csv(out_path)


##################################  main  ######################################
simu_path = glob.glob(os.path.join(proj_dir, 'result/simulation/simulation.random_single_mut/ic50_cnn_model_cell_*_drug_*.pkl'))
loci_anno_path = os.path.join(proj_dir, 'result/simulation/simulation.random_single_mut/ic50_cnn_model_alleles.pkl')

# location annotation
with open(loci_anno_path, 'rb') as f:
    data = pkl.load(f)
    loci, gene_index, drug_index, cell_index = data['alleles'], data['gene_index'], data['drug_index'], data['cell_index']
    loci2gene = [gene_index[l]+'_'+str(pos) for l, pos in loci]

for i in range(0, len(simu_path), 16):
    jobs = []
    for path in simu_path[i:i+16]:
        p = mp.Process(target=read_simu, args=(path, loci2gene, cell_index, drug_index,))
        jobs.append(p)
        p.start()
    for proc in jobs:
        proc.join()
