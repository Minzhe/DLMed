##############################################################################
###                         run.predict.svr.py                             ###
##############################################################################

import os
import numpy as np
import pandas as pd
from model.fcnn import fcnn
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error


proj_path = 'D:/projects/DLCell'
data_path = os.path.join(proj_path, 'data/curated/crispr.doubleKO.geno.pheno.csv')
model_path = os.path.join(proj_path, 'code/predict/model/fcnn.random.h5')

######################    main    #########################
geno_pheno = pd.read_csv(data_path)
X, label = geno_pheno.iloc[:, 0:-1], geno_pheno.iloc[:,-1]
# suffle label
y = np.random.choice(label, size=len(label), replace=False)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=1234)

### random guess
# coefs = []
# for i in range(10000):
#     truth = np.random.choice(y, size=int(len(y)*0.25), replace=False)
#     guess = np.random.choice(truth, size=len(truth), replace=True)
#     coefs.append(np.corrcoef(truth, guess)[0,1])

# print(np.mean(coefs), np.std(coefs))

### fcnn trained on shuffled datset
fcnn_random = fcnn(input_length=X.shape[1], learning_rate=0.0001, optimizer='RMSprop')
fcnn_random.train(X_train, y_train, model_name=model_path, validation_split=0, verbose=2)
pred_train = fcnn_random.predict(X_train)
pred_test = fcnn_random.predict(X_test)
r2_train = r2_score(y_train, pred_train[:,0])
r2_test = r2_score(y_test, pred_test[:,0])
