##############################################################################
###                            search.svr.py                               ###
##############################################################################

import os
import numpy as np
import pandas as pd
import pickle as pkl
import multiprocessing as mp
import datetime as dt
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.svm import SVR
from sklearn.externals import joblib

#############################  function  ##############################
def SVR_model(X_train, X_test, y_train, y_test, save_path, return_dict, C, epsilon, random_state=1234, verbose=1):
    '''
    Gradient Boosting Regressor model.
    '''
    if os.path.isdir(os.path.dirname(save_path)):
        model_path = '{}.C.{}.epsilon.{}.joblib'.format(save_path, C, epsilon)
    else:
        raise ValueError('Folder not exist!')
    # model
    print('Building model with SVR: C={}, epsilon={}'.format(C, epsilon))
    start_time = dt.datetime.now()
    svr = SVR(gamma='scale', C=C, epsilon=epsilon, verbose=verbose, cache_size=12000)
    svr.fit(X_train, y_train)
    joblib.dump(svr, model_path)
    # predict
    pred_train = svr.predict(X_train)
    pred_test = svr.predict(X_test)
    r2_train = r2_score(y_train, pred_train)
    r2_test = r2_score(y_test, pred_test)
    # result
    end_time = dt.datetime.now()
    print('\n' + '*'*60 + '\nModel {} finished! \nStart time: {}, End time: {}, Time used: {}'.format(os.path.basename(model_path), start_time, end_time, end_time-start_time))
    print('Traing r2: {}, Testing r2: {}\n'.format(r2_train, r2_test) + '*'*60 + '\n')
    return_dict[(C, epsilon)] = C, epsilon, r2_train, r2_test

#####################################    main    ######################################
if __name__ == '__main__':
    proj_path = '/work/bioinformatics/s418336/projects/DLMed'
    # proj_path = 'D:/projects/DLMed'
    data_path = os.path.join(proj_path, 'data/curated/Lung/merged/ccle_utsw.lung_MutExprCNV_cancergene_drug.array.pkl')
    model_path = os.path.join(proj_path, 'code/predict/model/svr/SVR.ccle_utsw.MutExprCNV')
    out_path = os.path.join(proj_path, 'result/genomic.drugsens/ccle_utsw.MutExprCNV.svr.result.3.csv')

    with open(data_path, 'rb') as f:
        data = pkl.load(f)
    X, y = data['X'], data['y']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1234)

    # hyperparameters
    C = [20, 30, 40, 50]
    epsilon = [1, 2]
    param = {'C': C, 'epsilon': epsilon}

    # multiprocess training
    manager = mp.Manager()
    return_dict = manager.dict()
    jobs = []
    for c_ in C:
        for e_ in epsilon:
            p = mp.Process(target=SVR_model, args=(X_train, X_test, y_train, y_test, model_path, return_dict, c_, e_,))
            jobs.append(p)
            p.start()
    
    for proc in jobs:
        proc.join()

    result = pd.DataFrame(list(return_dict.values()), columns=['C', 'eplison', 'r2_train', 'r2_test'])
    result['rank_test'] = result.r2_test.rank(ascending=False)
    result.to_csv(out_path, index=None)
