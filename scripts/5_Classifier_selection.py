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
import argparse


parser = argparse.ArgumentParser(description = "Classifier Selection")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--metadata','-m',help = 'input file : metadata')
parser.add_argument('--profile','-p',help = 'input file : microbial profile')
parser.add_argument('--exposure','-e',help = 'input param : the experiment group(exposure) of interest')
parser.add_argument('--group','-g',help = 'input param : the column name of experimental interest(group) in metadata')
parser.add_argument('--seed','-s',help = 'input param : random seed')
parser.add_argument('--output','-o',help = 'output file prefix: classifier selection result')
args = parser.parse_args()

#import data
metadata = pd.read_table(args.Workplace+args.metadata,sep = '\t',index_col = 0)
data = pd.read_table(args.Workplace+args.profile,sep = '\t',index_col=0)
data = data.loc[metadata.index,:]
data_group = np.array([1 if i== str(args.exposure) else 0 for i in metadata[str(args.group)]])
RANDOM_SEED = int(args.seed)

class machine_learning:
    
    def __init__(self):
        self.Method = {'Logisticl1':LogisticRegression(penalty='l1', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
                  'Logisticl2':LogisticRegression(penalty='l2', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
                  'DecisionTree':DecisionTreeClassifier(class_weight='balanced', random_state=RANDOM_SEED),
                  'RandomForest':RandomForestClassifier(oob_score=True, class_weight='balanced'),
                  'GradientBoost':GradientBoostingClassifier(random_state=RANDOM_SEED),
                  'KNeighbors':KNeighborsClassifier(n_neighbors=3),
                  'SVC':SVC(class_weight='balanced',random_state=RANDOM_SEED,probability = True)
                  }
        

    def scoring(self,clf, x, y):
        proba = clf.predict(x) 
        pred = np.array([1 if x>0.5 else 0 for x in proba])
        TP = ((pred==1) & (y==1)).sum()
        FP = ((pred==1) & (y==0)).sum()
        TN = ((pred==0) & (y==0)).sum()
        FN = ((pred==0) & (y==1)).sum()
        spe = TN / float(FP + TN)
        sen = recall_score(y, pred)
        precision = precision_score(y, pred)
        accuracy = accuracy_score(y, pred)
        auc = roc_auc_score(y, proba)
        f1 = f1_score(y, pred)
        return [sen, spe, precision,accuracy, f1, auc]
            
    def try_Classifiers(self,X, Y): 
        results=[]
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED)
        for method, clf in self.Method.items():
            self_scores = []
            cross_scores = []
            for train_index, test_index in skf.split(X, Y):
                X_train, X_test = X[train_index], X[test_index]
                Y_train, Y_test = Y[train_index], Y[test_index]
                clf.fit(X_train, Y_train)
                self_scores.append(self.scoring(clf, X_train, Y_train))
                cross_scores.append(self.scoring(clf, X_test, Y_test))
            self_scores = list(pd.DataFrame(self_scores).mean())
            cross_scores = list(pd.DataFrame(cross_scores).mean())
            row = [method]
            row.extend(self_scores)
            row.extend(cross_scores)
            results.append(row)
        select_results = pd.DataFrame(results, columns=['classifier', 'self_sensitivity', 'self_specifity',
                                                 'self_precision','self_accuracy', 'self_f1', 'self_auc', 'cross_sensitivity',
                                                 'cross_specifity', 'cross_precision', 'cross_accuracy', 'cross_f1',
                                                 'cross_auc'])
        return select_results

#classifier selection
ML = machine_learning()

result = ML.try_Classifiers(data.values,data_group)
result.to_csv(args.Workplace+args.output+"_classifier_selection.txt",sep = '\t')
print("FINISH")