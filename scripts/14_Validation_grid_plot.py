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
parser = argparse.ArgumentParser(description = "Internal validation grid plot")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--plot','-p',help = 'input file : grid plot table')
parser.add_argument('--output','-o',help = 'output file prefix: internal validation plot ')
args = parser.parse_args()

plot_heatmap = pd.read_table(args.Workplace+args.plot,sep = '\t',index_col = 0)

max_plot = np.array(plot_heatmap.max()).max()
min_plot = np.array(plot_heatmap.min()).min()
cohorts = plot_heatmap.columns[:-1]
lc = len(cohorts)

def validation_plot(plot_heatmap):
    sns.set(font_scale=1.5)
    grid_kws = {"height_ratios":(lc,1,1),"width_ratios":(lc,1)}
    f, axs= plt.subplots(3,2, gridspec_kw=grid_kws)

    sns.heatmap(plot_heatmap.iloc[0:lc,0:lc],cmap="YlGnBu", annot = True, ax = axs[0,0],vmin=min_plot,vmax=max_plot,cbar = False,fmt='.2')
    axs[0,0].xaxis.set_ticks_position('top')
    sns.heatmap(pd.DataFrame(plot_heatmap.loc['LOCO',cohorts]).T,cmap="YlGnBu", annot = True, ax = axs[2,0],vmin=min_plot,vmax=max_plot, cbar = False,xticklabels=False,fmt='.2')
    axs[2,0].set_yticklabels(["LOCO"],rotation='horizontal')
    sns.heatmap(pd.DataFrame(plot_heatmap.iloc[lc,0:lc]).T,cmap="YlGnBu", annot = True, ax = axs[1,0],vmin=min_plot,vmax=max_plot,cbar = False,xticklabels=False,fmt='.2')
    axs[1,0].set_yticklabels(["Average"],rotation='horizontal')
    sns.heatmap(pd.DataFrame(plot_heatmap.iloc[0:lc,lc]),cmap="YlGnBu", annot = True, ax = axs[0,1],vmin=min_plot,vmax=max_plot,cbar = False,yticklabels=False,fmt='.2')
    axs[0,1].xaxis.set_ticks_position('top')
    sns.heatmap(pd.DataFrame(plot_heatmap.loc['Average','Average'],index = ["Average"],columns = ["Average"]),cmap="YlGnBu", annot = True, ax = axs[1,1],vmin=min_plot,vmax=max_plot,cbar = False,xticklabels=False,yticklabels=False,fmt='.2')
    sns.heatmap(pd.DataFrame(plot_heatmap.loc['LOCO','Average'],index = ["LODO"],columns = ["Average"]),cmap="YlGnBu", annot = True, ax = axs[2,1],vmin=min_plot,vmax=max_plot,cbar = False,xticklabels = False,yticklabels=False,fmt='.2')
    plt.subplots_adjust(wspace =.05, hspace =.1)
    return plt

fig = validation_plot(plot_heatmap)
fig.savefig(args.Workplace+args.output+'_validation.pdf',bbox_inches = 'tight')

print("FINISH")