###################################################################################
###                            search.svr.cell.py                               ###
###################################################################################

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
def SVR_model(X_train, X_val, X_test, y_train, y_val, y_test, 
              kernel, gamma, C, epsilon, save_path, return_dict, random_state=1234, verbose=1):
    '''
    Gradient Boosting Regressor model.
    '''
    if os.path.isdir(os.path.dirname(save_path)):
        model_path = '{}.kernel.{}.gamma.{}.C.{}.epsilon.{}.joblib'.format(save_path, kernel, gamma, C, epsilon)
    else:
        raise ValueError('Folder not exist!')
    # model
    print('Building model with SVR: C={}, epsilon={}'.format(C, epsilon))
    start_time = dt.datetime.now()
    svr = SVR(kernel=kernel, gamma=gamma, C=C, epsilon=epsilon, verbose=verbose, cache_size=12000)
    svr.fit(X_train, y_train)
    joblib.dump(svr, model_path)

    # predict
    pred_train = svr.predict(X_train)
    pred_val = svr.predict(X_val)
    pred_test = svr.predict(X_test)
    r2_train = r2_score(y_train, pred_train)
    r2_val = r2_score(y_val, pred_val)
    r2_test = r2_score(y_test, pred_test)

    # result
    end_time = dt.datetime.now()
    print('\n' + '*'*60 + '\nModel {} finished! \nStart time: {}, End time: {}, Time used: {}'.format(os.path.basename(model_path), start_time, end_time, end_time-start_time))
    print('Traing r2: {}, Validation r2: {}, Testing r2: {}\n'.format(r2_train, r2_val, r2_test) + '*'*60 + '\n')
    return_dict[(kernel, gamma, C, epsilon)] = kernel, gamma, C, epsilon, r2_train, r2_val, r2_test

#####################################    main    ######################################
if __name__ == '__main__':
    proj_path = '/work/bioinformatics/s418336/projects/DLMed'
    # proj_path = 'D:/projects/DLMed'
    data_path = os.path.join(proj_path, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_MutGeneExpr_cancergene_drugRes.gene_level_cell.pkl')
    model_path = os.path.join(proj_path, 'code/predict/model_repo/svr/SVR.ccle_utsw.cell.MutExpr')
    out_path = os.path.join(proj_path, 'result/genomic.drugsens/ccle_utsw.cell.MutExpr.svr.r2_corr.rbf.csv')

    with open(data_path, 'rb') as f:
        data = pkl.load(f)
    (X_train, y_train), (X_val, y_val), (X_test, y_test) = data['train'], data['val'], data['inference']

    # hyperparameters
    kernel = ['rbf']
    gamma = ['auto']
    C = [10, 100, 1000]
    epsilon = [0.1, 0.5]
    param = {'C': C, 'epsilon': epsilon}

    # multiprocess training
    manager = mp.Manager()
    return_dict = manager.dict()
    jobs = []
    for k_ in kernel:
        for g_ in gamma:
            for c_ in C:
                for e_ in epsilon:
                    p = mp.Process(target=SVR_model, args=(X_train, X_val, X_test, y_train, y_val, y_test, k_, g_, c_, e_, model_path, return_dict,))
                    jobs.append(p)
                    p.start()
    
    for proc in jobs:
        proc.join()

    result = pd.DataFrame(list(return_dict.values()), columns=['kernel', 'gamma', 'C', 'eplison', 'r2_train', 'r2_val', 'r2_inference'])
    result['rank_val'] = result.r2_val.rank(ascending=False)
    result.to_csv(out_path, index=None)
