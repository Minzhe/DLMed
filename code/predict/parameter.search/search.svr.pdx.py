###################################################################################
###                            search.svr.pdx.py                                ###
###################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import pickle as pkl
import sys
sys.path.append(os.path.join(proj_dir, 'code'))
from model.param_grid_search import SVR_grid_search


data_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.MutExpr.cancergenes.drug_response.pkl')
model_prefix = os.path.join(proj_dir, 'code/predict/model/svr/SVR.pdx.MutAFExpr')
out_path = os.path.join(proj_dir, 'result/genomic.drugsens/pdx.MutAFExpr.svr.result.csv')

with open(data_path, 'rb') as f:
    data = pkl.load(f)
print(data.keys())

# hyperparameters
C = [2000, 3000, 4000, 6000]
epsilon = [0.1, 0.2, 0.4, 0.8]
param = {'C': C, 'epsilon': epsilon}