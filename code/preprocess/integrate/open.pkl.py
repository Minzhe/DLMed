import os
import numpy as np
import pickle as pkl
np.set_printoptions(threshold=np.nan)

proj_path = 'D:/projects/DLCell'
file_path = os.path.join(proj_path, 'data/curated/lung/loci.level/lung_MutExprCNV_cancergene_drug_loci.pkl')

with open(file_path, 'rb') as f:
    data = pkl.load(f)

mut_array = data['mutation_array']
expr_array = data['expr_array']
cnv_array = data['cnv_array']
drug_array = data['drug_array']
cell_index = data['cell_index']
gene_index = data['gene_index']
drug_index = data['drug_index']

print(expr_array.shape)
print(np.apply_along_axis(lambda x: sum(np.isnan(x)), 0, expr_array))
print(np.apply_along_axis(lambda x: sum(np.isnan(x)), 0, cnv_array))
print(drug_index)