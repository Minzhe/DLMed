##############################################################################
###                           search.gbr.py                                ###
##############################################################################

import os
import numpy as np
import pandas as pd
import pickle as pkl
import datetime as dt
from sklearn.externals import joblib
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.ensemble import GradientBoostingRegressor
import multiprocessing as mp

#############################  function  ##############################
def GBR_model(X_train, X_test, y_train, y_test, save_path, return_dict, loss, learning_rate, n_estimators, max_features, subsample=1, random_state=1234, verbose=1):
    '''
    Gradient Boosting Regressor model.
    '''
    if os.path.isdir(os.path.dirname(save_path)):
        model_path = '{}.loss.{}.lr.{}.n_estimator.{}.max_feature.{}.subsample{}.joblib'.format(save_path, loss, learning_rate, n_estimators, max_features, subsample)
    else:
        raise ValueError('Folder not exist!')
    # model
    print('Building model with GBR: loss={}, lr={}, n_estimator={}, max_feature={}, subsample={}'.format(loss, learning_rate, n_estimators, max_features, subsample))
    start_time = dt.datetime.now()
    gbr = GradientBoostingRegressor(loss=loss, 
                                    learning_rate=learning_rate,
                                    n_estimators=n_estimators,
                                    max_features=max_features,
                                    subsample=subsample,
                                    random_state=random_state,
                                    verbose=verbose)
    gbr.fit(X_train, y_train)
    joblib.dump(gbr, model_path)
    # predict
    pred_train = gbr.predict(X_train)
    pred_test = gbr.predict(X_test)
    r2_train = r2_score(y_train, pred_train)
    r2_test = r2_score(y_test, pred_test)
    # result
    end_time = dt.datetime.now()
    print('\n' + '*'*60 + '\nModel {} finished! \nStart time: {}, End time: {}, Time used: {}'.format(os.path.basename(model_path), start_time, end_time, end_time-start_time))
    print('Traing r2: {}, Testing r2: {}\n'.format(r2_train, r2_test) + '*'*60 + '\n')
    return_dict[(loss, learning_rate, n_estimators, max_features, subsample)] = loss, learning_rate, n_estimators, max_features, subsample, r2_train, r2_test

#############################  main  ##############################
if __name__ == '__main__':
    proj_path = '/work/bioinformatics/s418336/projects/DLMed'
    # proj_path = 'D:/projects/DLMed'
    data_path = os.path.join(proj_path, 'data/curated/Lung/merged/ccle_utsw.lung_MutExprCNV_cancergene_drug.array.pkl')
    model_path = os.path.join(proj_path, 'code/predict/model/gbr/GBR.ccle_utsw.MutExprCNV')
    out_path = os.path.join(proj_path, 'result/genomic.drugsens/ccle_utsw.MutExprCNV.gbr.result.3.csv')

    # read data
    with open(data_path, 'rb') as f:
        data = pkl.load(f)
    X, y = data['X'], data['y']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1234)

    # hyperparameters
    lr = [0.1, 0.2, 0.3, 0.5, 0.9]
    n_estimators = [300, 400, 600, 900]
    subsample = [1.0]
    max_features = ['auto']
    loss = ['ls', 'lad', 'huber', 'quantile']

    # multiprocess training
    manager = mp.Manager()
    return_dict = manager.dict()
    jobs = []
    for loss_ in loss:
        for lr_ in lr:
            for n_ in n_estimators:
                for f_ in max_features:
                    for sub_ in subsample:
                        p = mp.Process(target=GBR_model, args=(X_train, X_test, y_train, y_test, model_path, return_dict, loss_, lr_, n_, f_, sub_,))
                        jobs.append(p)
                        p.start()
    
    for proc in jobs:
        proc.join()

    result = pd.DataFrame(list(return_dict.values()), columns=['loss', 'learning_rate', 'n_estimators', 'max_features', 'subsample', 'r2_train', 'r2_test'])
    result['rank_test'] = result.r2_test.rank(ascending=False)
    result.to_csv(out_path, index=None)

