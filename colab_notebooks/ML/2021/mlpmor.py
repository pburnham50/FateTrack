# run multilayer perceptron for multi-output regression, input=features

# import various libraries
import numpy as np
import csv
import pandas as pd
import argparse
import sklearn
import keras
from sklearn.datasets import make_regression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RepeatedKFold
from keras.models import Sequential
from keras.layers import Dense, Conv2D, Flatten
from keras.losses import Huber
from pathlib import Path

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("dir", help="directory where training/test data is saved")
parser.add_argument("suff", help="data suffix (usually date)")
parser.add_argument("out", help="output suffix (usually same as suff)")
parser.add_argument("epochs", help="number of epochs for training")
args = parser.parse_args()

# get the dataset
def get_dataset(): #change this to load in data
    X_train = pd.read_csv(args.dir + "train_" + args.suff + ".csv", sep = ',', header = None)
    y_train = pd.read_csv(args.dir + "train_id_" + args.suff + ".csv", sep = ',', header = None)
    X_test = pd.read_csv(args.dir + "test_" + args.suff + ".csv", sep = ',', header = None)
    y_test = pd.read_csv(args.dir + "test_id_" + args.suff + ".csv", sep = ',', header = None)

    #y_train = np.ravel(y_train, order = "C") # have to do this to make data the right shape

    # scale data
    sclr = StandardScaler()
    X_train = sclr.fit_transform(X_train)
    X_test = sclr.transform(X_test)

    return X_train, y_train, X_test, y_test

#loss_huber = Huber(delta = 0.2)

# get the model
def get_model(n_inputs, n_outputs):
    model = Sequential()
    model.add(Dense(20, kernel_initializer='variance_scaling', activation='relu'))
    model.add(Dense(20, kernel_initializer='variance_scaling', activation='relu'))
    model.add(Dense(n_outputs, kernel_initializer='variance_scaling'))
    model.compile(loss='mse', optimizer='adam', metrics = ['mean_squared_error']) #loss = 'mse'
    return model

# load dataset
X_train, y_train, X_test, y_test = get_dataset()
n_inputs, n_outputs = X_train.shape[1], y_train.shape[1]
# get model
model = get_model(n_inputs, n_outputs)
# fit the model on all data
history = model.fit(X_train, y_train, validation_data = (X_test, y_test), verbose = 2, epochs = int(args.epochs))
# make a prediction for new data
y_pred = model.predict(X_test)
to_export = pd.DataFrame(data=y_pred)
to_export.to_csv(args.dir + "results_" + args.out + ".csv")
# write mse and validation mse out to csv
mse = pd.DataFrame(data=history.history['mean_squared_error'])
val_mse = pd.DataFrame(history.history['val_mean_squared_error'])
mse = pd.concat([mse, val_mse], axis=1)
mse.to_csv(args.dir + 'mse_' + args.out + '.csv')
  
