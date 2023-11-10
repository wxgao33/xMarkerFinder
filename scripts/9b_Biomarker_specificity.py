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
from scipy.stats import wilcoxon
from mpl_toolkits.axes_grid1 import ImageGrid


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

#import data
parser = argparse.ArgumentParser(description = "Biomarker specificity assessment & plot")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
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
parser.add_argument('--seed','-s',help = 'input param : random seed')
parser.add_argument('--output','-o',help = 'output file prefix: External test result & plot ')
args = parser.parse_args()


opt_biomarker = pd.read_table(args.Workplace+args.profile,sep = '\t',index_col=0)

test_metadata = pd.read_table(args.Workplace+args.external_metadata,sep = '\t',index_col = 0)
test_data = pd.read_table(args.Workplace+args.external_profile,sep = '\t',index_col=0)
test_data_group = np.array([0 if i== str(args.exposure) else 1 for i in test_metadata[str(args.group)]])
test_data = test_data.loc[:,test_data.columns.isin(opt_biomarker.columns)]
test_data = test_data.fillna(0)


ex_metadata = pd.read_table(args.Workplace+args.other_metadata,sep = '\t',index_col = 0)
ex_data = pd.read_table(args.Workplace+args.other_profile,sep = '\t',index_col=0)
ex_data = ex_data.loc[:,ex_data.columns.isin(opt_biomarker.columns)]
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

    def model_construction(self,data, data_group,params,SEED,k_fold):
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
        splitor = StratifiedKFold(n_splits=k_fold, shuffle=True,random_state=SEED) 
        clf = self.Method[opt_clf].set_params(**params)

        for train_index, test_index in splitor.split(data, data_group):
            y_train, y_test = data_group[train_index], data_group[test_index]
            X_train, X_test = np.array(data)[train_index], np.array(data)[test_index]
            clf.fit(X_train,y_train)
            
            probas = clf.predict_proba(X_test)
            pred = clf.predict(X_test)
            
        
            fpr, tpr, thresholds = roc_curve(y_test, probas[:, 1])
            roc_auc = auc(fpr, tpr)
            
            aucs.append(roc_auc)
         
            ### plot data
            tprs.append(interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            plot_data.append([fpr, tpr, 'ROC Fold %d(AUC = %0.2f)' %(i+1, roc_auc)])
            i += 1
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        
        return clf, aucs,mean_auc


ML = machine_learning()

cases = set(ex_metadata[args.group])
cases.remove(str(args.exposure))
auc_comparison = pd.DataFrame(columns = cases,index = range(1,11,1))

for i in cases:
    print(i+" testing")
    temp_set = list(set(ex_metadata[ex_metadata[args.group]==i][args.batch]))
    temp_meta = ex_metadata[ex_metadata[args.batch].isin(temp_set)]
    temp_group = np.array([0 if i== str(args.exposure) else 1 for i in temp_meta[args.group]])
    temp = ex_data[ex_data.index.isin(temp_meta.index)]
    for j in range(1,11,1):
        temp_result = ML.model_construction(temp,temp_group,best_param,j,5)
        auc_comparison.at[j,i]=temp_result[2]

test_groups = list(set(test_metadata[args.group]))
test_case_group = test_groups[0] if test_groups[1] == str(args.exposure) else test_groups[1]

auc_comparison.insert(0,test_case_group,np.nan)
for j in range(1,11,1):
    temp_result = ML.model_construction(test_data,test_data_group,best_param,j,5)
    auc_comparison.at[j,test_case_group]=temp_result[2]

auc_comparison.to_csv(args.Workplace+args.output+"_specificity_result.txt", sep = '\t')

p_values = []
first_column_data = auc_comparison.iloc[:, 0]

for i in range(1, len(auc_comparison.columns)):
    _, p_val = wilcoxon(first_column_data, auc_comparison.iloc[:, i], alternative='two-sided')
    p_values.append(p_val)

significance_markers = []
for p in p_values:
    if p < 0.001:
        significance_markers.append('***')
    elif p < 0.01:
        significance_markers.append('**')
    elif p < 0.05:
        significance_markers.append('*')
    else:
        significance_markers.append('ns')  


sns.set_theme(style="white")
fig, ax = plt.subplots(figsize=(10, 6))
sns.boxplot(data=auc_comparison, ax=ax)
sns.swarmplot(data=auc_comparison, ax=ax)

ax.set_ylabel('AUC')
ax.set_xticklabels(auc_comparison.columns)
ax.set_title("Biomarker specificity assessment", fontsize=20, fontweight='bold', pad=20)


y_max = auc_comparison.max().max()  

for i, marker in enumerate(significance_markers):
    y_text = y_max + 0.01 + i*0.02  
    #x_offset = 0.1 + i*0.05  
    x1, x2 = 0 , i+1  
    ax.plot([x1, x1, x2, x2], [y_text-0.01, y_text, y_text, y_text-0.01], lw=1.5, color='black')
    ax.text((x1+x2)*.5, y_text, marker, ha='center', va='center', color='black', fontsize=15)

plt.savefig(args.Workplace+args.output+'_specificity_auc.pdf',bbox_inches = 'tight')
plt.savefig(args.Workplace+args.output+'_specificity_auc.svg',bbox_inches = 'tight',format = 'svg')

print("FINISH")
