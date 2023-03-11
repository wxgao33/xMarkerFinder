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
from bayes_opt import BayesianOptimization, UtilityFunction 
from scipy import interp
import matplotlib.pyplot as plt
import argparse

#import data
parser = argparse.ArgumentParser(description = "Hyperparameter adjustment and model construction")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--metadata','-m',help = 'input file : metadata')
parser.add_argument('--profile','-p',help = 'input file : microbial profile')
parser.add_argument('--exposure','-e',help = 'input param : the experiment group(exposure) of interest')
parser.add_argument('--group','-g',help = 'input param : the column name of experimental interest(group) in metadata')
parser.add_argument('--classifier','-c',help = 'input param : selected classifier')
parser.add_argument('--hyperparameter','-r',help = 'input param : tuned hyperparameters')
parser.add_argument('--importance','-a',help = 'input file : feature importances')
parser.add_argument('--seed','-s',help = 'input param : random seed')
parser.add_argument('--output','-o',help = 'output file prefix: hyperparameter adjustment and model construction result')
args = parser.parse_args()


metadata = pd.read_table(args.Workplace+args.metadata,sep = '\t',index_col = 0)
opt_biomarker = pd.read_table(args.Workplace+args.profile,sep = '\t',index_col=0)
data_group = np.array([1 if i== str(args.exposure) else 0 for i in metadata[str(args.group)]])
importances = pd.read_table(args.Workplace+args.importance,sep='\t',index_col=0).T
RANDOM_SEED = int(args.seed)
opt_clf = args.classifier

opt_biomarker=opt_biomarker.fillna(0)

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
        self.Method = {'LRl1':LogisticRegression(penalty='l1', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
                  'LRl2':LogisticRegression(penalty='l2', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
                  'DT':DecisionTreeClassifier(class_weight='balanced', random_state=RANDOM_SEED),
                  'RF':RandomForestClassifier(oob_score=True, class_weight='balanced',random_state=RANDOM_SEED),
                  'GB':GradientBoostingClassifier(random_state=RANDOM_SEED),
                  'KNN':KNeighborsClassifier(n_neighbors=3),
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

ML = machine_learning()

#Minimal panel selection

feat_import = importances.index[importances['importances_mean'].argsort()]
feat_comparison = pd.DataFrame(index=list(feat_import)[::-1],columns=['AUC'])

feat_select = []
i = len(feat_import)-1
while i>-1:
    feat_select.append(feat_import[i])
    data_select = opt_biomarker.loc[:,feat_select]
    temp_result = ML.model_construction(data_select,data_group,best_param,5)
    feat_comparison.iloc[len(feat_import)-1-i,0] = temp_result[2]

    i=i-1

feat_comparison.to_csv(args.Workplace+args.output+"_minimal_panel.txt", sep = '\t')

plot = pd.read_table(args.Workplace+args.output+"_minimal_panel.txt", sep = '\t')

#plot
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.xticks(np.arange(1,len(feat_comparison.index)+1,1))

plt.xlabel('Number of biomarkers', fontsize=13)
plt.ylabel('AUC Value', fontsize=13)
plt.plot(np.arange(1,len(feat_comparison.index)+1,1),plot.AUC,"o:")

plt.savefig(args.Workplace+args.output+"_"+args.classifier+"_minimal_panel.pdf",bbox_inches = 'tight')

print("FINISH")


