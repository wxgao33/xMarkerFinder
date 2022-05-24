#!/usr/bin/env python3
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import argparse

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
#plt.rcParams['font.size'] = 8

#import data
parser = argparse.ArgumentParser(description = "Marker importance plot")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--plot','-p',help = 'input file : feature importance table')
parser.add_argument('--output','-o',help = 'output file prefix: feature importance plot ')
args = parser.parse_args()

plot= pd.read_table(args.Workplace+args.plot,sep = '\t',index_col = 0).T

#def function
def plot_feature_imps(plot):
    perm_sorted_idx = plot['importances_mean'].argsort()

    result = pd.DataFrame(index = plot.index,columns = range(1,11))
    for i in range(len(plot.index)):
        for j in range(10):
            result.iloc[i,j] = re.split(r"[ ]+",plot['importances'][i][1:-1])[j]
    result = result.replace(r'^\s*$', np.nan, regex=True).astype('float')

    fig,ax = plt.subplots()
    ax.boxplot(result.iloc[perm_sorted_idx].T,vert=False,)
    ax.set_title('Permutation feature importance')
    ax.set_yticklabels(result.index[perm_sorted_idx])
    return plt

#plot
feature_imp = plot_feature_imps(plot)
feature_imp.savefig(args.Workplace+args.output+'_marker_importance.pdf',bbox_inches = 'tight')

print("FINISH")