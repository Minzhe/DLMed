###################################################################################
###                          predict.drugsens.ml.py                             ###
###################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import sys
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
sys.path.append(os.path.join(proj_dir, 'code'))
import utility.plot as p
from sklearn.externals import joblib
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.svm import SVR
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor

#############################    main   ################################
data_path = os.path.join(proj_dir, 'data/curated/Lung/merged/merged.lung_Expr_cancergene_drug.array.pkl')
model_path = os.path.join(proj_dir, 'code/predict/model/gbr')
fig_path = os.path.join(proj_dir, 'result/genomic.drugsens/gbr.png')

with open(data_path, 'rb') as f:
    data = pkl.load(f)
X, y, label = data['X'], data['y'], data['label']
cmap = list(map(lambda x: {1:'r', 2:'b', 3:'g'}[x], label))

# model
# svr = SVR(verbose=True)
# rfr = RandomForestRegressor()
gbr = joblib.load(os.path.join(model_path, 'GBR.lr.0.5.n_estimator.500.max_feature.log2.joblib'))

# train
print('Predicting with GBR model ...')
y_gbr = gbr.predict(X)
print(r2_score(y[label == 1], y_gbr[label == 1]))
print(r2_score(y[label == 2], y_gbr[label == 2]))
print(r2_score(y[label == 3], y_gbr[label == 3]))

# plot
# f, ax = plt.subplots(figsize=(18,15))
# p.plotScatter(y, y_gbr, ax=ax, title='GBR', xlabel='Truth', ylabel='Prediction', alpha=0.05, cmap=cmap)
# plt.savefig(fig_path)

