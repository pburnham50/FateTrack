# run random forest classifier

# import various libraries
import numpy as np
import csv
import pandas as pd
import argparse
import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.inspection import permutation_importance
from sklearn import metrics
from pathlib import Path

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("dir", help="directory where training/test data is saved")
parser.add_argument("suff", help="data suffix (usually date)")
parser.add_argument("out", help = "output suffix (usually same as suff)")
args = parser.parse_args()

# import data
#data = pd.read_csv("~/Desktop/ml_data/data_0404/improps/data_0413.csv", sep = ',', header = None)
#target = pd.read_csv("~/Desktop/ml_data/data_0404/improps/target_0413.csv", sep = ',', header = None)

#X_train, X_test, y_train, y_test = train_test_split(data, target, random_state=0)
X_train = pd.read_csv(args.dir + "train_" + args.suff + ".csv", sep = ',', header = None)
y_train = pd.read_csv(args.dir + "train_id_" + args.suff + ".csv", sep = ',', header = None)
X_test = pd.read_csv(args.dir + "test_" + args.suff + ".csv", sep = ',', header = None)
y_test = pd.read_csv(args.dir + "test_id_" + args.suff + ".csv", sep = ',', header = None)

y_train = np.ravel(y_train, order = "C") # have to do this to make data the right shape

# random forest
model = RandomForestClassifier(n_estimators=500)
model.fit(X_train, y_train)
y_pred = model.predict(X_test)

#importances = model.feature_importances_
#importances = pd.DataFrame(importances)
#importances.to_csv(args.dir + "importances_" + args.out + ".csv")

#perm_imp = permutation_importance(model, X_train, y_train)
#perm_imp = pd.DataFrame(perm_imp.importances_mean)
#perm_imp.to_csv("~/Desktop/ml_data/data_0415/perm_importances_0415_75.csv")

# export results
to_export = pd.DataFrame()
to_export['Pred'] = y_pred
to_export['True'] = y_test
to_export.to_csv(args.dir + "results_" + args.out + ".csv")

y_pred = model.predict_proba(X_test)
y_pred = pd.DataFrame(y_pred)
y_pred.to_csv(args.dir + "proba_" + args.out + ".csv")

#rf_metrics = metrics.classification_report(y_pred, y_test)
#filename = Path('~/Desktop/ml_data/data_0404/improps/metrics_0413.txt')
#filename.touch(exist_ok=True)  # will create file, if it exists will do nothing
#metrics_file = open(filename, 'w')
#metrics_file = open(r'~/Desktop/ml_data/data_0404/improps/metrics_0413.txt', 'x')
#metrics_file.write(rf_metrics)
#metrics_file.close()
