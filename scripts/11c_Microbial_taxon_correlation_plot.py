#!/usr/bin/env python3

import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
import argparse

#import data
parser = argparse.ArgumentParser(description = "Microbial taxon correlation plot")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--correlation','-c',help = 'input file : correlation profile')
parser.add_argument('--pvalues','-p',help = 'input file : pvalue profile')
parser.add_argument('--threshold','-t',help = 'input param : correlation threshold')
parser.add_argument('--output','-o',help = 'output file prefix: feature importance result ')
args = parser.parse_args()

corr = pd.read_csv(args.Workplace+args.correlation,sep = '\t',index_col=0)
pval = pd.read_csv(args.Workplace+args.pvalues,sep = '\t',index_col=0)
corr_threshold = float(args.threshold)

pval[pval>0.05]=-1
pval[pval>0]=1
pval[pval<0]=0

corr_filter = corr * pval
corr_filter[abs(corr_filter)==0]=0
corr_filter[abs(corr_filter)<corr_threshold]=0
corr_filter=corr_filter[corr_filter.sum(axis=1) != 0]
corr_filter=corr_filter.loc[:,corr_filter.index]


corr_filter.to_csv(args.Workplace+args.output+"_correlation.csv")

print("FINISH")