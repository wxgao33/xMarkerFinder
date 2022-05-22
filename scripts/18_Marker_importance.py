#!/usr/bin/env python3

import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC,LinearSVC
from sklearn.metrics import roc_curve,auc,recall_score,precision_score,f1_score,accuracy_score,roc_auc_score
from sklearn.inspection import permutation_importance
from scipy import interp
import argparse

#import data
parser = argparse.ArgumentParser(description = "Marker importance")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--metadata','-m',help = 'input file : metadata')
parser.add_argument('--profile','-p',help = 'input file : microbial profile')
parser.add_argument('--exposure','-e',help = 'input param : the experiment group(exposure) of interest')
parser.add_argument('--group','-g',help = 'input param : the column name of experimental interest(group) in metadata')
parser.add_argument('--classifier','-c',help = 'input param : selected classifier')
parser.add_argument('--hyperparameter','-r',help = 'input param : tuned hyperparameters')
parser.add_argument('--seed','-s',help = 'input param : random seed')
parser.add_argument('--output','-o',help = 'output file prefix: feature importance result ')
args = parser.parse_args()

metadata = pd.read_table(args.Workplace+args.metadata,sep = '\t',index_col = 0)
opt_biomarker = pd.read_table(args.Workplace+args.profile,sep = '\t',index_col=0)
data_group = np.array([1 if i== str(args.exposure) else 0 for i in metadata[str(args.group)]])
RANDOM_SEED = int(args.seed)
opt_clf = args.classifier

params = {}
file = open(args.Workplace+args.hyperparameter,'r')
for line in file.readlines():
    line = line.strip()
    k = line.split(' ')[0]
    v = line.split(' ')[1]
    params[k] = v
file.close()

best_param= [{k: int(v) if v and '.' not in v else float(v) if v else None for k, v in d.items()}for d in [params]][0]

dict = {'Logistic(l1)':LogisticRegression(penalty='l1', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
'Logistic(l2)':LogisticRegression(penalty='l2', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
'DecisionTree':DecisionTreeClassifier(class_weight='balanced', random_state=RANDOM_SEED),
'RandomForest':RandomForestClassifier(oob_score=True, class_weight='balanced'),
'GradientBoost':GradientBoostingClassifier(random_state=RANDOM_SEED),
'KNeighbors':KNeighborsClassifier(n_neighbors=3),
'SVC':SVC(class_weight='balanced',random_state=RANDOM_SEED,probability = True)}

#calculate permutation importance
def feature_imps(param,data,y_data):
    clf = dict[opt_clf].set_params(**param).fit(data.values,y_data)
    result = permutation_importance(clf, data.values, y_data,n_repeats = 10, random_state = RANDOM_SEED)
    perm_sorted_idx = result.importances_mean.argsort()

    return result,perm_sorted_idx

result,_ = feature_imps(best_param,opt_biomarker,data_group)
feature_perm_imp = pd.DataFrame.from_dict(result,orient = 'index',columns = opt_biomarker.columns)
feature_perm_imp.to_csv(args.Workplace+args.output+"_"+args.classifier+"_marker_importance.txt",sep = '\t')

print("FINISH")