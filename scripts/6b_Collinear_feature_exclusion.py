#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse

#imprt data

parser = argparse.ArgumentParser(description = "Collinear feature exclusion")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--profile','-p',help = 'input file : microbial profile')
parser.add_argument('--threshold','-t',help = 'input param : corr threshold')
parser.add_argument('--output','-o',help = 'output file prefix: collinear feature exclusion result')
args = parser.parse_args()

#import data
data = pd.read_table(args.Workplace+args.profile,sep = '\t',index_col=0)
corr_threshold = float(args.threshold)


#def function
def feature_correlation(data,corr_threshold):
    corrdf = data.corr(method = 'spearman')
    corrdf[corrdf==1]= 0
    i=0
    while i <len(corrdf.index):
        corrdf = corrdf[abs(corrdf.iloc[i,])<corr_threshold]
        corrdf = corrdf.loc[:,(corrdf.index)]
        i +=1
    self_effective_feature = data.loc[:,corrdf.index]
    return corrdf, self_effective_feature

#Stage3B Collinear feature exclusion
corrdf, self_effective_feature = feature_correlation(data,corr_threshold)

corrdf.to_csv(args.Workplace+args.output+"_feature_correlation.txt", sep = '\t')
self_effective_feature.to_csv(args.Workplace+args.output+"_uncorrelated_effective_feature.txt", sep = '\t')

print("FINISH")
