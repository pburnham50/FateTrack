# run convolutional neural network for multi-output regression, input=images

# import various libraries
import numpy as np
import csv
import pandas as pd
import argparse
import sklearn
import keras
from sklearn.datasets import make_regression
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RepeatedKFold
from keras.utils import to_categorical
from keras.models import Sequential
from keras.layers import Dense, Conv2D, Flatten
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
    X_train, y_train = make_regression(n_samples=1000, n_features=10, n_informative=5, n_targets=3, random_state=2)
    return X, y

# get the model
def get_model(n_inputs, n_outputs):
    model = Sequential()
    model.add(Conv2D(64, kernel_size=3, activation='relu', input_shape=(28,28,1), padding = 'same')) #change input shape based on images
    model.add(Conv2D(32, kernel_size=3, activation='relu'))
    model.add(Flatten())
    model.add(Dense(n_outputs))
    model.compile(loss = 'mse', optimizer = 'adam', metrics = ['mse'])
    return model

# load dataset
X, y = get_dataset()
n_inputs, n_outputs = X.shape[1], y.shape[1]
# get model
model = get_model(n_inputs, n_outputs)
# fit the model on all data
model.fit(X, y, verbose=0, epochs=100)
# make a prediction for new data
row = [-0.99859353,2.19284309,-0.42632569,-0.21043258,-1.13655612,-0.55671602,-0.63169045,-0.87625098,-0.99445578,-0.3677487]
newX = asarray([row])
yhat = model.predict(newX)
print('Predicted: %s' % yhat[0])
