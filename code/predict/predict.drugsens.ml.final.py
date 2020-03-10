#######################################################################################
###                         predict.drugsens.ml.final.py                            ###
#######################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import sys
import numpy as np
import pandas as pd
import pickle as pkl
from sklearn.externals import joblib
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.svm import SVR
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor

###############################  main  ##############################
data_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_MutGeneExpr_cancergene_drugRes.gene_level_cell.pkl')
with open(data_path, 'rb') as f:
    data = pkl.load(f)
X_train, y_train = data['train'][0], data['train'][1]
X_val, y_val = data['val'][0], data['val'][1]
X_inf, y_inf = data['inference'][0], data['inference'][1]

### svr model

