#!/usr/bin/env python3

import numpy as np
import pandas as pd
import math
from bioinfokit import analys, visuz
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

parser = argparse.ArgumentParser(description = "Differential signature volcano plot")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--input','-i',help = 'input file : differential significance table')
parser.add_argument('--threshold','-t',help = 'input param : pvalue threshold')
parser.add_argument('--output','-o',help = 'output file prefix: differential volvano plot')
args = parser.parse_args()

#import data
pvalues = pd.read_table(args.Workplace+args.input,sep = '\t')
threshold = float(args.threshold)

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

#prepare plot data
plot = pd.DataFrame(pvalues.loc[:,['pval','coef']])
plot['group'] = '#b8b8b8'
plot['log10pval'] = [-math.log(x,10) for x in list(plot['pval'])]
pval_threshold = -math.log(threshold,10)
plot.loc[(plot.coef>0)&(plot.log10pval>pval_threshold),'group'] = 'tab:red'
plot.loc[(plot.coef<0)&(plot.log10pval>pval_threshold),'group'] = 'tab:blue'
plot.loc[plot.log10pval<pval_threshold,'group'] = '#b8b8b8'


xmax=max(abs(plot.coef))
xmin=0-xmax
ymin=0
ymax=max(plot.log10pval)+2

def volcano_plot(plot_data):
    fig = plt.figure(figsize=(7,8))
    ax = fig.add_subplot()
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='')
    ax.scatter(plot['coef'], plot['log10pval'], s=10, c=plot['group'],)
    ax.set_ylabel('-Log10(P value)')
    ax.set_xlabel('Coefficient')
    ax.spines['right'].set_visible(False) 
    ax.spines['top'].set_visible(False) 
    ax.vlines(0, 0, 5, color='dimgrey',linestyle='dashed', linewidth=1) 
    ax.hlines(pval_threshold, -0.0025, 0.0025, color='dimgrey',linestyle='dashed', linewidth=1)

    ax.scatter(plot['coef'], plot['log10pval'], s=10, c=plot['group'],marker = 'o')
    return plt

#plot
fig = volcano_plot(plot)
fig.savefig(args.Workplace+args.output+"_differential_volcano.pdf",bbox_inches = 'tight')
print("FINISH")