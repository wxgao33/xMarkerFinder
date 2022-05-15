#!/usr/bin/env python3
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC,LinearSVC
from sklearn.metrics import roc_curve,auc,recall_score,precision_score,f1_score,accuracy_score,roc_auc_score
from numpy import interp
import matplotlib.pyplot as plt
import argparse


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

#import data
parser = argparse.ArgumentParser(description = "Cross validation plot")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--metadata','-m',help = 'input file : metadata')
parser.add_argument('--profile','-p',help = 'input file : microbial profile')
parser.add_argument('--exposure','-e',help = 'input param : the experiment group(exposure) of interest')
parser.add_argument('--group','-g',help = 'input param : the column name of experimental interest(group) in metadata')
parser.add_argument('--classifier','-c',help = 'input param : selected classifier')
parser.add_argument('--hyperparameter','-r',help = 'input param : tuned hyperparameters')
parser.add_argument('--seed','-s',help = 'input param : random seed')
parser.add_argument('--output','-o',help = 'output file prefix: cross validation auc plot ')
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


#def function
class machine_learning:   
    def __init__(self):
        self.Method = {'Logistic(l1)':LogisticRegression(penalty='l1', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
                  'Logistic(l2)':LogisticRegression(penalty='l2', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
                  'DecisionTree':DecisionTreeClassifier(class_weight='balanced', random_state=RANDOM_SEED),
                  'RandomForest':RandomForestClassifier(oob_score=True, class_weight='balanced'),
                  'GradientBoost':GradientBoostingClassifier(random_state=RANDOM_SEED),
                  'KNeighbors':KNeighborsClassifier(n_neighbors=3),
                  'SVC':SVC(class_weight='balanced',random_state=RANDOM_SEED,probability = True)
                  }

    def model_construction(self,data, data_group, params, k_fold):
        aucs = []
        tprs = []
        mean_fpr = np.linspace(0, 1, 100)
        plot_data = []
        i = 0
        sens = []
        spes = []
        pres = []
        f1s = []
        accus = []
        splitor = StratifiedKFold(n_splits=k_fold, shuffle=True,random_state=RANDOM_SEED) 
        clf = self.Method[opt_clf].set_params(**params)
        
        for train_index, test_index in splitor.split(data, data_group):
            y_train, y_test = data_group[train_index], data_group[test_index]
            X_train, X_test = np.array(data)[train_index], np.array(data)[test_index]
            clf.fit(X_train,y_train)
            
            probas = clf.predict_proba(X_test)
            pred = clf.predict(X_test)
            sen = recall_score(y_test,pred)
            sens.append(sen)
            TP = ((pred==1) & (y_test==1)).sum()
            FP = ((pred==1) & (y_test==0)).sum()
            TN = ((pred==0) & (y_test==0)).sum()
            FN = ((pred==0) & (y_test==1)).sum()
            spe = TN / float(FP + TN)
            spes.append(spe)
            
            pre = precision_score(y_test, pred)
            pres.append(pre)
            f1 = f1_score(y_test, pred)
            f1s.append(f1)
            fpr, tpr, thresholds = roc_curve(y_test, probas[:, 1])
            roc_auc = auc(fpr, tpr)
            spes.append(spe)
            aucs.append(roc_auc)
            accu = accuracy_score(y_test, pred)
            accus.append(accu)
            
            
            ### plot data
            tprs.append(interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            plot_data.append([fpr, tpr, 'ROC Fold %d(AUC = %0.2f)' %(i+1, roc_auc)])
            i += 1
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        mean_spe = np.mean(spes)
        mean_sen = np.mean(sens)
        mean_pre = np.mean(pres)
        mean_f1 = np.mean(f1s)
        mean_accu = np.mean(accus)
        
        return clf, mean_auc,mean_spe,mean_sen,mean_pre,mean_f1,mean_accu,(plot_data, mean_fpr, mean_tpr, tprs,aucs, np.std(aucs))

    def plot_auc(self,plot_data):
        lines, mean_fpr, mean_tpr, tprs, aucs, std_auc = model_result[7]
        mean_auc = model_result[1]
        font1 = {'family' : 'Arial','weight' : 'normal','size': 12}
        font2 = {'family' : 'Arial','weight' : 'normal','size': 16}
        plt.figure(1, (6, 5.5), dpi=300)
        fig, ax = plt.subplots(figsize=(6,5.5))
        for fpr, tpr, label in lines:
            plt.plot(fpr, tpr, lw=1, alpha=0.3)
        plt.plot([0, 1], [0, 1], ls='--', lw=2, color='r')
        plt.plot(mean_fpr, mean_tpr, color='b', lw=2.5,label='Mean ROC(AUC=%0.2f $\pm$ %0.2f)'%(mean_auc, std_auc))

        std_tpr = np.std(tprs, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,label=r'$\pm$ 1 std. dev.')

        plt.xlabel('False Positive Rate', fontsize=13)
        plt.ylabel('True Positive Rate', fontsize=13)

        plt.legend(loc='lower right',prop=font1)
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('False Positive Rate', font2)
        plt.ylabel('True Positive Rate', font2)
        plt.yticks(fontproperties = 'Arial', size = 13)
        plt.xticks(fontproperties = 'Arial', size = 13)

        return plt

ML = machine_learning()

#prepare for plot
model_result = ML.model_construction(opt_biomarker,data_group, best_param, k_fold=5)

#plot
auc_fig = ML.plot_auc(model_result)
auc_fig.savefig(args.Workplace+args.output+"_"+args.classifier+"_cross_validation_auc.pdf",bbox_inches = 'tight')

print("FINISH")
