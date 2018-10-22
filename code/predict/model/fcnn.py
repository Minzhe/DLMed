#######################################################################
###                            fcnn.py                              ###
#######################################################################

import numpy as np
from keras.models import Model, load_model
from keras.layers import Input, Dense, BatchNormalization
from keras.layers.advanced_activations import ReLU
from keras.optimizers import RMSprop, Adam
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras import backend as K


###########################  function  ################################
def batchGenerator(X, y, batch_size, input_length):
    batch_X = np.zeros((batch_size, input_length))


def r_square(y_true, y_pred):
    SS_res = K.sum(K.square(y_true - y_pred))
    SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
    return (1 - SS_res/(SS_tot + K.epsilon()))

def pearson_coef(y_true, y_pred):
    y_true_c = y_true - K.mean(y_true)
    y_pred_c = y_pred - K.mean(y_pred)
    return K.sum(y_true_c * y_pred_c) / (K.sqrt(K.sum(K.square(y_true_c)) * K.sum(K.square(y_pred_c))) + K.epsilon())


###########################  model  ################################
class fcnn(object):
    '''
    Fully connected nerual network
    '''
    def __init__(self, input_length, learning_rate=1e-4, optimizer='Adam'):
        self.input_length = input_length
        self.learning_rate = learning_rate
        if optimizer == 'Adam':
            self.optimizer = Adam(lr=self.learning_rate)
        elif optimizer == 'RMSprop':
            self.optimizer = RMSprop(lr=self.learning_rate, decay=1e-6)
        else:
            raise ValueError('Unrecognizable optimizer. Either Adam or RMSprop.')
        self.model = self._model()
    

    def _model(self):
        print('Initilizing fully connected nerual netowrk model ...', end='', flush=True)
        inputs = Input(shape=(self.input_length,), name='input')

        dense1 = Dense(256) (inputs)
        bn1 = BatchNormalization() (dense1)
        relu1 = ReLU(max_value=6.0) (bn1)

        dense2 = Dense(256) (relu1)
        bn2 = BatchNormalization() (dense2)
        relu2 = ReLU(max_value=6.0) (bn2)

        dense3 = Dense(256) (relu2)
        bn3 = BatchNormalization() (dense3)
        relu3 = ReLU(max_value=6.0) (bn3)

        output = Dense(1) (relu3)

        model = Model(inputs=inputs, outputs=output)
        model.compile(loss='mse', optimizer=self.optimizer, metrics=['mse', r_square, pearson_coef])
        print(' Done\nModel structure summary:', flush=True)
        print(model.summary())

        return model
    
    def train(self, X_train, y_train, model_name, validation_split=0.0, validation_data=None, batch_size=32, epochs=200, verbose=1):
        print('Start training neural network ... ', end='', flush=True)
        early_stopper = EarlyStopping(patience=5, verbose=1)
        check_pointer = ModelCheckpoint(model_name, verbose=1, save_best_only=True)
        result = self.model.fit(X_train, y_train, 
                                validation_split=validation_split, 
                                validation_data=validation_data,
                                batch_size=batch_size, 
                                epochs=epochs, 
                                verbose=verbose,
                                shuffle=True,
                                callbacks=[early_stopper, check_pointer])
        self.model = load_model(model_name)
        print('Done')


    def loadModel(self, path):
        print('Loading trained neural network model ... ', end='', flush=True)
        self.model = load_model(path)
        print('Done')
    
    def predict(self, X):
        print('Predicting with nerual network model ... ', flush=True)
        y = self.model.predict(X, verbose=1)
        return y