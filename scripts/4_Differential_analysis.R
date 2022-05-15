#!/usr/bin/env Rscript

if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}

if (TRUE){
  option_list = list(
    make_option(c("-W", "--workplace"), type="character", default="/Workplace/",
                help="Input workplace [default %default]"),
    make_option(c("-m", "--metadata"), type="character", default="metadata.txt",
                help="metadata file [default %default]"),
    make_option(c("-f", "--feature_profile"), type="character", default="feature_abundance.txt",
                help="feature abundance profile[default %default]"),
    make_option(c("-g", "--group"), type="character", default="Group",
                help="the column name of experimental interest(group) in metadata [default %default]"),
    make_option(c("-b", "--batch"), type="character", default="Cohort",
                help="major confounder [default %default]"),
    make_option(c("-c", "--covariates"), type="character", default="covariates.txt",
                help="minor confounders [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output file prefix [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))

  # 显示输入输出确认是否正确
  print(paste("The workplace is ", opts$workplace,  sep = ""))
  print(paste("The metadata file is ", opts$metadata,  sep = ""))
  print(paste("The profile of features is ", opts$feature_profile,  sep = ""))
  print(paste("The exposure is ", opts$exposure,  sep = ""))
  print(paste("The batch is ", opts$batch,  sep = ""))
  print(paste("The covariates file is ", opts$covariates,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}


package_list = c("MMUPHin","dplyr")
###CHECK PACKAGES
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

#import data
metadata <- read.csv(file = paste(opts$workplace,opts$metadata,sep=''),sep = '\t',header =  TRUE,row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
feat_abd <- read.csv(file = paste(opts$workplace,opts$feature_profile,sep=''),sep = '\t',header =  TRUE,row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
feat_abd = t(feat_abd)
feat_abd[is.na(feat_abd)]<-0
batch <- opts$batch
exposure <- opts$group
covariates <- row.names(read.table(file = paste(opts$workplace,opts$covariates,sep=''),row.names = 1))


#Differential analysis
fit_meta <- lm_meta(feature_abd = feat_abd,
                    exposure = exposure,
                    batch = batch,
                    covariates = covariates,
                    data = metadata)

meta_results <- fit_meta$meta_fits
maaslin_results <- fit_meta$maaslin_fits

write.table(meta_results,file = paste(opts$workplace,opts$output,"_feature_significance.txt",sep=''),sep = '\t')
write.table(maaslin_results,file = paste(opts$workplace,opts$output,"_feature_significance_single_cohort.txt",sep=''),sep = '\t')

feat_diff <- t(feat_abd[(meta_results %>% filter(pval<0.05))[,'feature'],])
write.table(feat_diff,file = paste(opts$workplace,opts$output,"_differential_signature.txt",sep=''),sep = '\t')

print("FINISH")

