#####################################################################################
###                            integrate.for.svr.py                               ###
#####################################################################################

proj_path = 'D:/projects/DLMed'
import os
import sys
import numpy as np
import pickle as pkl
sys.path.append(os.path.join(proj_path, 'code'))
from utility.formatData import buildTrainingArray

file_path = os.path.join(proj_path, 'data/curated/Lung/merged/merged.lung_MutExprCNV_cancergene_drug_loci.pkl')
out_path = os.path.join(proj_path, 'data/curated/Lung/merged/merged.lung_MutExprCNV_cancergene_drug_loci.array.pkl')

#########################  main  #########################
with open(file_path, 'rb') as f:
    data = pkl.load(f)

mut_array = data['mutation_array']
drug_array = data['drug_array']
cell_index = data['cell_index']
gene_index = data['gene_index']
drug_index = data['drug_index']

X, y = buildTrainingArray(mut_array, drug_array)
print(X.shape)
print(y.shape)
with open(out_path, 'wb') as f:
    pkl.dump({'X_csr':X, 'y':y}, file=f)
