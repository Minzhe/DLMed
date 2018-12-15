###################################################################################
###                         predict.drugsens.fcnn.py                            ###
###################################################################################

import os
import numpy as np
import pandas as pd
from model.fcnn import fcnn
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error


proj_path = 'D:/projects/DLCell'
# data_path = os.path.join(proj_path, 'data/curated/lung/lung.genomic.drug.csv')
data_path = os.path.join(proj_path, 'data/curated/lung/lung.cnv.drug.csv')
model_path = os.path.join(proj_path, 'code/predict/model/fcnn.cnv.drugsens.h5')

####################################    main   ######################################
genome_drug = pd.read_csv(data_path, index_col=0)
X, y = genome_drug.iloc[:,:-1], genome_drug.iloc[:,-1]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=1234)

### fcnn
fcnn_reg = fcnn(input_length=X.shape[1], learning_rate=0.0001, optimizer='RMSprop')
fcnn_reg.train(X_train, y_train, model_name=model_path, validation_data=(X_test, y_test), verbose=2)
pred_train = fcnn_reg.predict(X_train)
pred_test = fcnn_reg.predict(X_test)
r2_train = r2_score(y_train, pred_train[:,0])
r2_test = r2_score(y_test, pred_test[:,0])