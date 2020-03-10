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
# import seaborn as sns; sns.set()
sys.path.append(os.path.join(proj_dir, 'code'))
# import utility.plot as p
from sklearn.externals import joblib
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.svm import SVR
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor

def plot_correlation(model, X_train, y_train, X_val, y_val, X_test, y_test, name, path):
    print('Predicting training ...')
    pred_train = model.predict(X_train)
    print('Predicting validation ...')
    pred_val = model.predict(X_val)
    print('Predicting testing ...')
    pred_test = model.predict(X_test)
    r2 = r2_score(y_train, pred_train), r2_score(y_val, pred_val), r2_score(y_test, pred_test)
    # plot 
    lims = [-11.0, 11.0]
    loc = [-10, 10]
    s = 1
    f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True, figsize=(27,8))
    ax1.scatter(pred_train, y_train, s=s)
    ax2.scatter(pred_val, y_val, s=s)
    ax3.scatter(pred_test, y_test, s=s)
    f.suptitle('Prediction of {} Model'.format(name), fontsize=15)
    ax1.set_title('Training')
    ax2.set_title('Validation')
    ax3.set_title('Testing')
    for i, ax in enumerate([ax1, ax2, ax3]):
        ax.plot(lims, lims, 'r--', alpha=0.75, zorder=0)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.text(-10, 9.5, r'$r^2$: {}'.format(round(r2[i], 3)), fontsize=12)
        ax.set(xlabel='Predicted', ylabel='Truth')
    f.savefig(path)

#############################    main   ################################
data_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_MutGeneExpr_cancergene_drugRes.gene_level_cell.pkl')
model_path = os.path.join(proj_dir, 'code/predict/model_repo')
fig_svr_path = os.path.join(proj_dir, 'result/genomic.drugsens/fig/predict_cell.svr.png')
fig_rfr_path = os.path.join(proj_dir, 'result/genomic.drugsens/fig/predict_cell.rfr.png')
fig_gbr_path = os.path.join(proj_dir, 'result/genomic.drugsens/fig/predict_cell.gbr.png')


with open(data_path, 'rb') as f:
    data = pkl.load(f)
(X_train, y_train), (X_val, y_val), (X_test, y_test) = data['train'], data['val'], data['inference']

# model
# svr = SVR(verbose=True)
print('Loading SVR ...')
svr = joblib.load(os.path.join(model_path, 'svr/SVR.ccle_utsw.cell.MutExpr.kernel.linear.gamma.auto.C.200.epsilon.0.1.joblib'))
print('Loading RFR ...')
rfr = joblib.load(os.path.join(model_path, 'rfr/RFR.ccle_utsw.cell.MutExpr.n_estimator.400.max_feature.auto.criterion.mse.max_depth.None.joblib'))
print('Loading GBR ...')
gbr = joblib.load(os.path.join(model_path, 'gbr/GBR.ccle_utsw.cell.MutExpr.loss.ls.lr.0.2.n_estimator.1000.max_feature.auto.subsample1.0.joblib'))

plot_correlation(svr, X_train, y_train, X_val, y_val, X_test, y_test, 'Support Vector Machine', fig_svr_path)
plot_correlation(rfr, X_train, y_train, X_val, y_val, X_test, y_test, 'Random Forest', fig_rfr_path)
plot_correlation(gbr, X_train, y_train, X_val, y_val, X_test, y_test, 'Gradient Boosting', fig_gbr_path)
