import pandas as pd
import argparse

#import data
parser = argparse.ArgumentParser(description = "Convert files for Step 20")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--profile','-p',help = 'input file : microbial profile')
parser.add_argument('--select','-s',help = 'input file : selected signatures for Step 20')
parser.add_argument('--output','-o',help = 'output file prefix: convert result')
args = parser.parse_args()

data = pd.read_table(args.Workplace+args.profile,sep = '\t',index_col = 0)
select = pd.read_table(args.Workplace+args.select,sep = '\t',index_col = 0)
data = data.loc[:,select.columns]
data = data.T 
data.index.name="#OTU ID"
data = data.fillna(0)
data.to_csv(args.Workplace+args.output+"_convert.tsv",sep = '\t')
print("FINISH")