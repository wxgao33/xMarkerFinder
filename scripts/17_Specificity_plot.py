#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import argparse

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.sans-serif'] = 'Arial'

#import data
parser = argparse.ArgumentParser(description = "Specificity plot")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--plot','-p',help = 'input file : grid plot table')
parser.add_argument('--output','-o',help = 'output file prefix: internal validation plot ')
args = parser.parse_args()

auc_comparison= pd.read_table(args.Workplace+args.plot,sep = '\t',index_col = 0)

fig = plt.figure(figsize=(8,6))
fig = sns.set_theme(style="white")
fig = sns.boxplot(data = auc_comparison)
fig = sns.swarmplot(data = auc_comparison)
fig.set_ylabel('AUC')
fig.set_xticklabels(auc_comparison.columns)

plt.savefig(args.Workplace+args.output+'_specificity_auc.pdf',bbox_inches = 'tight')

print("FINISH")