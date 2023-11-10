#!/usr/bin/env python3
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC,LinearSVC
from sklearn.metrics import roc_curve,auc,recall_score,precision_score,precision_recall_curve,f1_score,accuracy_score,roc_auc_score,matthews_corrcoef
from numpy import interp
import matplotlib.pyplot as plt
import argparse


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

#import data
parser = argparse.ArgumentParser(description = "External test result & plot")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--metadata','-m',help = 'input file : training set metadata')
parser.add_argument('--profile','-p',help = 'input file : training set microbial profile')
parser.add_argument('--external_metadata','-a',help = 'input file : test set metadata')
parser.add_argument('--external_profile','-x',help = 'input file : test set microbial profile')
parser.add_argument('--exposure','-e',help = 'input param : the experiment group(exposure) of interest')
parser.add_argument('--group','-g',help = 'input param : the column name of experimental interest(group) in metadata')
parser.add_argument('--classifier','-c',help = 'input param : selected classifier')
parser.add_argument('--hyperparameter','-r',help = 'input param : tuned hyperparameters')
parser.add_argument('--seed','-s',help = 'input param : random seed')
parser.add_argument('--output','-o',help = 'output file prefix: External test result & plot ')
args = parser.parse_args()

metadata = pd.read_table(args.Workplace+args.metadata,sep = '\t',index_col = 0)
opt_biomarker = pd.read_table(args.Workplace+args.profile,sep = '\t',index_col=0)
data_group = np.array([1 if i== str(args.exposure) else 0 for i in metadata[str(args.group)]])

ex_metadata = pd.read_table(args.Workplace+args.external_metadata,sep = '\t',index_col = 0)
ex_data = pd.read_table(args.Workplace+args.external_profile,sep = '\t',index_col=0)
ex_data_group = np.array([1 if i== str(args.exposure) else 0 for i in ex_metadata[str(args.group)]])
ex_data = pd.DataFrame(ex_data,columns=opt_biomarker.columns)
ex_data = ex_data.fillna(0)

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
        self.Method = {'LRl1':LogisticRegression(penalty='l1', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
                  'LRl2':LogisticRegression(penalty='l2', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
                  'DT':DecisionTreeClassifier(class_weight='balanced', random_state=RANDOM_SEED),
                  'RF':RandomForestClassifier(oob_score=True, class_weight='balanced', random_state=RANDOM_SEED),
                  'GB':GradientBoostingClassifier(random_state=RANDOM_SEED),
                  'KNN':KNeighborsClassifier(n_neighbors=3),
                  'SVC':SVC(class_weight='balanced',random_state=RANDOM_SEED,probability = True)
                  }

    def test_model(self,data, data_group, ex_data,ex_data_group,params):
        clf = self.Method[opt_clf].set_params(**params).fit(opt_biomarker,data_group)
        probas = clf.predict_proba(ex_data)
        pred = clf.predict(ex_data)
        sen = recall_score(ex_data_group,pred)

        precision, recall, _ = precision_recall_curve(ex_data_group,pred)
        aupr = auc(recall,precision)
        mcc = matthews_corrcoef(ex_data_group,pred)

        TP = ((pred==1) & (ex_data_group==1)).sum()
        FP = ((pred==1) & (ex_data_group==0)).sum()
        TN = ((pred==0) & (ex_data_group==0)).sum()
        FN = ((pred==0) & (ex_data_group==1)).sum()
        spe = TN / float(FP + TN)
        
        pre = precision_score(ex_data_group, pred)
        f1 = f1_score(ex_data_group, pred)
        fpr, tpr, thresholds = roc_curve(ex_data_group, probas[:, 1])
        roc_auc = auc(fpr, tpr)
        accu = accuracy_score(ex_data_group, pred)
        
        return clf, roc_auc,aupr, mcc,spe,sen,pre,f1,accu,fpr, tpr

    def plot_auc(self,fpr,tpr,auc):
        font1 = {'weight' : 'normal','size': 12}
        font2 = {'weight' : 'normal','size': 16}

        plt.figure(1, (6, 5.5), dpi=300)
        fig, ax = plt.subplots(figsize=(6,5.5))

        plt.plot([0, 1], [0, 1], ls='--', lw=2, color='r')
        plt.plot(fpr, tpr, color='b', lw=2.5,label='AUC=%0.2f'%(roc_auc))


        plt.xlabel('False Positive Rate', fontsize=13)
        plt.ylabel('True Positive Rate', fontsize=13)

        plt.legend(loc='lower right',prop=font1)
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('False Positive Rate', font2)
        plt.ylabel('True Positive Rate', font2)
        plt.yticks(size = 13)
        plt.xticks(size = 13)

        return plt

ML = machine_learning()

#Test
_, roc_auc,aupr,mcc,spe,sen,pre,f1,accu,fpr, tpr= ML.test_model(opt_biomarker,data_group,ex_data,ex_data_group, best_param)
test_result = pd.DataFrame([roc_auc,aupr,mcc,spe,sen,pre,f1,accu],columns = ['Value'],index = ['AUC','AUPR','MCC','Specificity','Sensitivity','Precision','F1 Score','Accuracy'])
test_result.to_csv(args.Workplace+args.output+"_"+args.classifier+"_test_result.txt", sep = '\t')


#Plot
auc_fig = ML.plot_auc(fpr,tpr,roc_auc)
auc_fig.title("ROC curve of external test", fontsize=20, fontweight='bold', pad=20)
auc_fig.savefig(args.Workplace+args.output+"_"+args.classifier+"_test_auc.pdf",bbox_inches = 'tight')
auc_fig.savefig(args.Workplace+args.output+"_"+args.classifier+"_test_auc.svg",bbox_inches = 'tight',format = 'svg')

print("FINISH")
