##############################################################################
###                         run.predict.svr.py                             ###
##############################################################################

import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.svm import SVR

if __name__ == '__main__':
    proj_path = '/work/bioinformatics/s418336/projects/DLCell'
    # proj_path = 'D:/projects/DLCell'
    data_path = os.path.join(proj_path, 'data/curated/crispr.doubleKO.geno.pheno.csv')
    out_path = os.path.join(proj_path, 'result/crispr.DKO/svr.result.csv')

    geno_pheno = pd.read_csv(data_path)
    X, y = geno_pheno.iloc[:, 0:-1], geno_pheno.iloc[:,-1]

    ### SVR
    svr = SVR(verbose=True)
    param = {'C':[0.1, 1, 10], 'epsilon':[0.01, 0.1, 0.2]}
    clf = GridSearchCV(estimator=svr, param_grid=param, scoring={'r2':'r2', 'mse':'neg_mean_squared_error'}, cv=4, n_jobs=-1, refit=False, verbose=10, return_train_score=True)
    clf.fit(X, y)

    res = pd.DataFrame(clf.cv_results_)
    res.to_csv(out_path, index=None)
