###################################################################################
###                            search.rfr.cell.py                               ###
###################################################################################

import os
import numpy as np
import pandas as pd
import pickle as pkl
import multiprocessing as mp
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.ensemble import RandomForestRegressor
from sklearn.externals import joblib
import datetime as dt

#############################  function  ##############################
def RFR_model(X_train, X_val, X_test, y_train, y_val, y_test, 
              n_estimators, max_features, criterion, max_depth, 
              save_path, return_dict, random_state=1234, verbose=1):
    '''
    Gradient Boosting Regressor model.
    '''
    if os.path.isdir(os.path.dirname(save_path)):
        model_path = '{}.n_estimator.{}.max_feature.{}.criterion.{}.max_depth.{}.joblib'.format(save_path, n_estimators, max_features, criterion, max_depth)
    else:
        raise ValueError('Folder not exist!')
    # model
    print('Building model with RFR: n_estimator={}, max_feature={}, criterion={}, max_depth={}'.format(n_estimators, max_features, criterion, max_depth))
    start_time = dt.datetime.now()
    rfr = RandomForestRegressor(n_estimators=n_estimators,
                                max_features=max_features,
                                criterion=criterion,
                                max_depth=max_depth,
                                random_state=random_state,
                                verbose=verbose)
    rfr.fit(X_train, y_train)
    joblib.dump(rfr, model_path)

    # predict
    pred_train = rfr.predict(X_train)
    pred_val = rfr.predict(X_val)
    pred_test = rfr.predict(X_test)
    r2_train = r2_score(y_train, pred_train)
    r2_val = r2_score(y_val, pred_val)
    r2_test = r2_score(y_test, pred_test)

    # result
    end_time = dt.datetime.now()
    print('\n' + '*'*60 + '\nModel {} finished! \nStart time: {}, End time: {}, Time used: {}'.format(os.path.basename(model_path), start_time, end_time, end_time-start_time))
    print('Traing r2: {}, Validation r2: {}, Testing r2: {}\n'.format(r2_train, r2_val, r2_test) + '*'*60 + '\n')
    return_dict[(n_estimators, max_features, criterion, max_depth)] = n_estimators, max_features, criterion, max_depth, r2_train, r2_val, r2_test

#####################################    main    ######################################
if __name__ == '__main__':
    proj_path = '/work/bioinformatics/s418336/projects/DLMed'
    # proj_path = 'D:/projects/DLMed'
    data_path = os.path.join(proj_path, 'data/curated/Lung/merged/ccle_utsw.lung_MutExprCNV_cancergene_drug.array.fixcell.70.20.10.pkl')
    model_path = os.path.join(proj_path, 'code/predict/model/rfr/RFR.ccle_utsw.cell.MutExprCNV')
    out_path = os.path.join(proj_path, 'result/genomic.drugsens/ccle_utsw.cell.MutExprCNV.rfr.result.csv')

    with open(data_path, 'rb') as f:
        data = pkl.load(f)
    (X_train, y_train), (X_val, y_val), (X_test, y_test) = data['train'], data['val'], data['test']

    # hyperparameters
    n_estimators = [20, 40, 60, 80, 100, 200, 400]
    max_features = ['auto']
    criterion = ['mse']
    max_depth = [None]

    # multiprocess training
    manager = mp.Manager()
    return_dict = manager.dict()
    jobs = []
    for n_ in n_estimators:
        for f_ in max_features:
            for c_ in criterion:
                for d_ in max_depth:
                    p = mp.Process(target=RFR_model, args=(X_train, X_val, X_test, y_train, y_val, y_test, n_, f_, c_, d_, model_path, return_dict,))
                    jobs.append(p)
                    p.start()
    
    for proc in jobs:
        proc.join()

    result = pd.DataFrame(list(return_dict.values()), columns=['n_estimators', 'max_features', 'criterion', 'max_depth', 'r2_train', 'r2_val', 'r2_test'])
    result['rank_test'] = result.r2_test.rank(ascending=False)
    result.to_csv(out_path, index=None)