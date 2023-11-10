#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.model_selection import GroupKFold, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC,LinearSVC
from sklearn.metrics import roc_curve, auc
from numpy import interp
import copy
import argparse

#import data

parser = argparse.ArgumentParser(description = "Feature effectiveness evaluation")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--metadata','-m',help = 'input file : metadata')
parser.add_argument('--profile','-p',help = 'input file : microbial profile')
parser.add_argument('--exposure','-e',help = 'input param : the experiment group(exposure) of interest')
parser.add_argument('--group','-g',help = 'input param : the column name of experimental interest(group) in metadata')
parser.add_argument('--cohort','-b',help = 'input param : batch(cohort) column name')
parser.add_argument('--classifier','-c',help = 'input param : selected classifier')
parser.add_argument('--seed','-s',help = 'input param : random seed')
parser.add_argument('--threshold','-t',help = 'input param : auc threshold')
parser.add_argument('--output','-o',help = 'output file prefix: feature effectiveness evaluation result')
args = parser.parse_args()

#import data
metadata = pd.read_table(args.Workplace+args.metadata,sep = '\t',index_col = 0)
feat_diff = pd.read_table(args.Workplace+args.profile,sep = '\t',index_col = 0)
data_group = np.array([1 if i== str(args.exposure) else 0 for i in metadata[str(args.group)]])
data_cohort = metadata[args.cohort]
RANDOM_SEED = int(args.seed)
opt_clf = args.classifier
auc_threshold = float(args.threshold)


#def function
class machine_learning:
    
    def __init__(self):
        self.Method = {'LRl1':LogisticRegression(penalty='l1', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
                  'LRl2':LogisticRegression(penalty='l2', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
                  'DT':DecisionTreeClassifier(class_weight='balanced', random_state=RANDOM_SEED),
                  'RF':RandomForestClassifier(oob_score=True, class_weight='balanced', random_state=RANDOM_SEED),
                  'GB':GradientBoostingClassifier(random_state=RANDOM_SEED),
                  'KNN':KNeighborsClassifier(n_neighbors=3),
                  'SVC':SVC(class_weight='balanced',random_state=RANDOM_SEED,probability = True)
                  }

    def crossvalidation_auc(self,data, data_group, data_cohort):
        aucs = []
        tprs = []
        mean_fpr = np.linspace(0, 1, 100)
        i = 0
        
        clf = self.Method[opt_clf]
        if len(set(data_cohort)) > 0:
            splitor = GroupKFold(n_splits=len(set(data_cohort))).split(data, data_group, data_cohort)
        else:
            splitor = StratifiedKFold(n_splits=5, shuffle=True,random_state=RANDOM_SEED).split(data, data_group)
        
        for train_index, test_index in splitor:
            y_train, y_test = data_group[train_index], data_group[test_index]
            X_train, X_test = np.array(data)[train_index], np.array(data)[test_index]
            
            probas = clf.fit(X_train, y_train).predict_proba(X_test)
            fpr, tpr, thresholds = roc_curve(y_test, probas[:, 1])
            roc_auc = auc(fpr, tpr)
            aucs.append(roc_auc)
            tprs.append(interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            i += 1
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        return mean_auc

    def get_feature_auc(self,data, data_group, data_cohort):
        feature_auc = pd.DataFrame(columns = ['auc'],index = list(data.columns))
        for i in list(data.columns):
            feature_auc.loc[i,:] = self.crossvalidation_auc(pd.DataFrame(data.loc[:,i]),data_group, data_cohort)
        feature_auc = feature_auc.sort_values('auc',ascending = False)
        return feature_auc


#Stage3A Feature effectiveness evaluation
ML = machine_learning()
feature_auc= ML.get_feature_auc(feat_diff, data_group,data_cohort)
effective_feature = feat_diff.loc[:,feature_auc[feature_auc['auc']>auc_threshold].index]

feature_auc.to_csv(args.Workplace+args.output+"_feature_auc.txt",sep = '\t')
effective_feature.to_csv(args.Workplace+args.output+'_effective_feature.txt', sep = '\t')

print("FINISH")

