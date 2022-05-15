#!/usr/bin/env Rscript

if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}

if (TRUE){
  option_list = list(
    make_option(c("-W", "--workplace"), type="character", default="/Users/evelyn/Downloads/research/01Project/protocols/data/",
                help="Input workplace [default %default]"),
    make_option(c("-m", "--metadata"), type="character", default="metadata.txt",
                help="metadata file [default %default]"),
    make_option(c("-p", "--feature_profile"), type="character", default="feature_abundance.txt",
                help="feature abundance profile[default %default]"),
    make_option(c("-s", "--seed"), type="character", default="0",
                help="random_seed[default %default]"),
    make_option(c("-g", "--exposure"), type="character", default="Group",
                help="the column name of experimental interest(group) [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output file prefix [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))

  # 显示输入输出确认是否正确
  print(paste("The workplace is ", opts$workplace,  sep = ""))
  print(paste("The metadata file is ", opts$metadata,  sep = ""))
  print(paste("The profile of features is ", opts$feature_profile,  sep = ""))
  print(paste("The exposure column is ", opts$exposure,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}

package_list = c("dplyr","Boruta")
###CHECK PACKAGES
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

#import data
metadata <- read.table(file = paste(opts$workplace,opts$metadata,sep=''),sep = '\t',header =  TRUE,row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
feat_abd <- read.csv(file = paste(opts$workplace,opts$feature_profile,sep=''),sep = '\t',header =  TRUE,row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
#feat_abd = t(feat_abd)
set.seed(opts$seed)
exposure <- opts$exposure
Group <- factor(metadata[,exposure])

boruta <- Boruta(x=feat_abd, y=Group, pValue=0.05, mcAdj=T,maxRuns=1000)
print(table(boruta$finalDecision))

#extract feature
boruta.finalVarsWithTentative <- data.frame(Item=getSelectedAttributes(boruta, withTentative = F), Type="Boruta_with_tentative")
feature <- boruta.finalVarsWithTentative$Item
boruta.variable.imp.use <- boruta$ImpHistory[,feature]
feature_importance <- apply(boruta.variable.imp.use,2,mean)
feature_importance <- data.frame(sort(feature_importance,decreasing = TRUE))
selected_features <- feat_abd[,rownames(feature_importance)]
write.table(feature_importance,file = paste(opts$workplace,opts$output,"_boruta_feature_imps.txt",sep=''),sep = '\t')
write.table(selected_features,file = paste(opts$workplace,opts$output,"_boruta_selected_features.txt",sep=''),sep = '\t')

print("FINISH")

