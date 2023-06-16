# run support vector machine on HCR static images

# import various libraries
import numpy as np
import csv
import pandas as pd
import sklearn
from sklearn import svm
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV

# import data
X_train = pd.read_csv("~/Desktop/ml_data/data_0404/improps/train_0407_2.csv", sep = ',', header = None)
y_train = pd.read_csv("~/Desktop/ml_data/data_0404/improps/train_id_0407_2.csv", sep = ',', header = None)
X_test = pd.read_csv("~/Desktop/ml_data/data_0404/improps/test_0407_2.csv", sep = ',', header = None)
y_test = pd.read_csv("~/Desktop/ml_data/data_0404/improps/test_id_0407_2.csv", sep = ',', header = None)

y_train = np.ravel(y_train, order = "C") # have to do this to make data the right shape

# scaling pipeline, create and train SVM
pipe = make_pipeline(StandardScaler(), svm.SVC(C = 2, kernel='rbf', coef0 = 1, degree = 2, gamma = 0.005))
pipe.fit(X_train, y_train)  # apply scaling on training data
#pipe.score(X_test, y_test)  # apply scaling on testing data, without leaking training data.
test_results = pipe.predict(X_test)
# create and train SVM
#classifier = svm.SVC()
#classifier.fit(X_train, y_train)

#cross-validation (I performed this in a jupyter notebook)
#param_grid = {'svc__C': [1, 2, 3, 4, 5, 10, 50], 'svc__gamma': [0.0001, 0.0005, 0.001, 0.005]}
#grid = GridSearchCV(pipe, param_grid)
#grid.fit(X_train, y_train)
#print(grid.best_params_)

# test SVM and export results
#test_results = classifier.predict(X_test)
to_export = pd.DataFrame()
to_export['Pred'] = test_results
to_export['True'] = y_test
to_export.to_csv("~/Desktop/ml_data/data_0404/improps/results_0407_7.csv")
