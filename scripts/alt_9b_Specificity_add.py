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
import seaborn as sns
from mpl_toolkits.axes_grid1 import ImageGrid


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

#import data
parser = argparse.ArgumentParser(description = "Alternative specificity assessment & plot")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--metadata','-m',help = 'input file : training set metadata')
parser.add_argument('--profile','-p',help = 'input file : optimal biomarkers')
parser.add_argument('--external_metadata','-q',help = 'input file : test set metadata')
parser.add_argument('--external_profile','-l',help = 'input file : test set microbial profile')
parser.add_argument('--other_metadata','-a',help = 'input file : test set metadata')
parser.add_argument('--other_profile','-x',help = 'input file : test set microbial profile')
parser.add_argument('--exposure','-e',help = 'input param : the control group name')
parser.add_argument('--group','-g',help = 'input param : the column name of experimental interest(group) in metadata')
parser.add_argument('--batch','-b',help = 'input param : the column name of cohort(dataset)')
parser.add_argument('--classifier','-c',help = 'input param : selected classifier')
parser.add_argument('--hyperparameter','-r',help = 'input param : tuned hyperparameters')
parser.add_argument('--number','-n',help = 'input param : number of samples for validation each time')
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
number = int(args.number)

other_metadata = pd.read_table(args.Workplace+args.other_metadata,sep = '\t',index_col = 0)
other_data = pd.read_table(args.Workplace+args.other_profile,sep = '\t',index_col=0)
other_data = pd.DataFrame(other_data,columns=opt_biomarker.columns)
other_data = other_data.fillna(0)

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
                  'RF':RandomForestClassifier(oob_score=True, class_weight='balanced'),
                  'GB':GradientBoostingClassifier(random_state=RANDOM_SEED),
                  'KNN':KNeighborsClassifier(n_neighbors=3),
                  'SVC':SVC(class_weight='balanced',random_state=RANDOM_SEED,probability = True)
                  }

    def test_model(self,data, data_group, ex_data,ex_data_group,params):
        clf = self.Method[opt_clf].set_params(**params).fit(opt_biomarker,data_group)
        probas = clf.predict_proba(ex_data)
        pred = clf.predict(ex_data)
        sen = recall_score(ex_data_group,pred)

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
        
        return clf, roc_auc

ML = machine_learning()

#Specificity assessment
cases = set(other_metadata[args.group])
cases.remove(str(args.exposure))
column_name = []
for i in cases:
    column_name.append(i+"_control")
    column_name.append(i+"_case")
auc_comparison = pd.DataFrame(columns = column_name,index = range(1,11,1))
ls=[]
ls.extend("0"*number)

ex_data_group_add = np.append(ex_data_group,ls).astype(int)


for i in cases:
    print(i+" testing")
    temp_set = list(set(other_metadata[other_metadata[args.group]==i][args.batch]))
    temp_meta = other_metadata[other_metadata[args.batch].isin(temp_set)]
    temp_control = other_data[other_data.index.isin(temp_meta[temp_meta[args.group]==args.exposure].index)]
    temp_case = other_data[other_data.index.isin(temp_meta[temp_meta[args.group]!=args.exposure].index)]
    
    for j in range (1,11,1):
        add_control = temp_control.sample(n=number,axis=0,random_state=j)
        ex_add_control = pd.concat((ex_data,temp_control.sample(n=number,axis=0,random_state=j)),axis=0)
        ex_add_case = pd.concat((ex_data,temp_case.sample(n=number,axis=0,random_state=j)),axis=0)
        _, roc_auc_control= ML.test_model(opt_biomarker,data_group,ex_add_control,ex_data_group_add, best_param)
        _, roc_auc_case= ML.test_model(opt_biomarker,data_group,ex_add_case,ex_data_group_add, best_param)
        auc_comparison.at[j,i+"_control"]=roc_auc_control
        auc_comparison.at[j,i+"_case"]=roc_auc_case

auc_comparison.to_csv(args.Workplace+args.output+"_specificity_add_result.txt", sep = '\t')

#fig = plt.figure(figsize=(8,6))
fig = sns.set_theme(style="white")
fig = sns.boxplot(data = auc_comparison)
fig = sns.swarmplot(data = auc_comparison)
fig.set_ylabel('AUC')
fig.set_xticklabels(auc_comparison.columns,rotation=30)
#plt.tight_layout()
plt.savefig(args.Workplace+args.output+'_specificity_add_auc.pdf',bbox_inches = 'tight')

print("FINISH")