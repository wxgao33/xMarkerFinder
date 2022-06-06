#!/usr/bin/env Rscript

if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}

if (TRUE){
  option_list = list(
    make_option(c("-W", "--workplace"), type="character", default="/Workplace/",
                help="Input workplace [default %default]"),
    make_option(c("-p", "--profile"), type="character", default="abundance.txt",
                help="feature abundance profile[default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output file prefix [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))

  # 显示输入输出确认是否正确
  print(paste("The workplace is ", opts$workplace,  sep = ""))
  print(paste("The microbial profile is ", opts$profile,  sep = ""))
}


package_list = c("dplyr")
###CHECK PACKAGES
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

#import data and preprocess
feat_abd <- read.csv(file = paste(opts$workplace,opts$profile,sep=''),sep = '\t',header =  TRUE,row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
feat_abd[is.na(feat_abd)]<-0 
rela <- t(apply(feat_abd,1, function(x){x/sum(x)}))
write.table(rela,file = paste(opts$workplace,opts$output,"_relative_abundance.txt",sep=''),sep = '\t')

print("FINISH")

