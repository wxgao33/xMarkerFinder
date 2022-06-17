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
from sklearn.inspection import permutation_importance
from numpy import interp
import copy
import argparse


parser = argparse.ArgumentParser(description = "Recursive feature elimination")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--metadata','-m',help = 'input file : metadata')
parser.add_argument('--profile','-p',help = 'input file : microbial profile')
parser.add_argument('--exposure','-e',help = 'input param : the experiment group(exposure) of interest')
parser.add_argument('--group','-g',help = 'input param : the column name of experimental interest(group) in metadata')
parser.add_argument('--classifier','-c',help = 'input param : selected classifier')
parser.add_argument('--seed','-s',help = 'input param : random seed')
parser.add_argument('--output','-o',help = 'output file prefix: recursive feature elimination result')
args = parser.parse_args()


#import data
metadata = pd.read_table(args.Workplace+args.metadata,sep = '\t',index_col = 0)
data = pd.read_table(args.Workplace+args.profile,sep = '\t',index_col = 0)
data_group = np.array([1 if i== str(args.exposure) else 0 for i in metadata[str(args.group)]])
RANDOM_SEED = int(args.seed)
opt_clf = args.classifier

#def function
class machine_learning:
    
    def __init__(self):
        self.Method = {'LRl1':LogisticRegression(penalty='l1', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
                  'LRl2':LogisticRegression(penalty='l2', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
                  'DT':DecisionTreeClassifier(class_weight='balanced', random_state=RANDOM_SEED),
                  'RF':RandomForestClassifier(oob_score=True, class_weight='balanced'),
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
        if len(data_cohort) > 0:
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

    def feature_imps(self,data,data_group):
        clf = self.Method[opt_clf].fit(data.values,data_group)
        result = permutation_importance(clf,data.values,data_group,n_repeats = 10, random_state = RANDOM_SEED)
        perm_sorted_idx = result.importances_mean.argsort()

        return result,perm_sorted_idx

    def RFE(self,data, data_group):
        select = list(data.columns)
        best_auc = 0
        best_features = []
        aucs = []
        while(len(select)>1):
            aucs = []
            temp = copy.deepcopy(select)
            auc = self.crossvalidation_auc(data.loc[:,temp],data_group,[])
            aucs.append([temp,auc])
            _,perm_sorted_idx = self.feature_imps(data.loc[:,temp],data_group)
            temp.remove(temp[perm_sorted_idx.tolist()[0]])
            select = temp
            if auc >= best_auc:
                best_auc = auc
                best_feature = select
            print('Feature number: ',len(select), 'Best_AUC: ',round(best_auc, 3), round(auc, 3))
        return best_auc, best_feature

    def RFE2(self,data, data_group):
        select = list(data.columns)
        best_auc = 0
        best_features = []
        while(len(select)>1):
            aucs = []
            for i in select:
                temp = copy.deepcopy(select)
                temp.remove(i)
                auc = self.crossvalidation_auc(data.loc[:,temp],data_group,[])
                aucs.append([temp,auc])
            select, auc = sorted(aucs,key = lambda x:x[1],reverse = True)[0]
            if auc >= best_auc:
                best_auc = auc
                best_feature = select
            print('Feature number: ',len(select), 'Best_AUC: ',round(best_auc, 3), round(auc, 3))
        return best_auc, best_feature


#Stage3C Recursive feature elimination
ML = machine_learning()
_,optimal_biomarkers = ML.RFE2(data,data_group)
data.loc[:,optimal_biomarkers].to_csv(args.Workplace+args.output+"_candidate_biomarker.txt", sep = '\t')
print("FINISH")
