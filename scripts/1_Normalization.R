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
    make_option(c("-m", "--method"), type="character", default="REL",
                help="normalization method[default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output file prefix [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))

  # 显示输入输出确认是否正确
  print(paste("The workplace is ", opts$workplace,  sep = ""))
  print(paste("The microbial profile is ", opts$profile,  sep = ""))
  print(paste("The normalization method is ", opts$method, sep = ""))
}


package_list = c("dplyr","compositions","edgeR")
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
norm_method <- opts$method

#normalization
if(norm_method == "REL"){
  data_norm <- t(apply(feat_abd,1, function(x){x/sum(x)}))
} else if (norm_method =="AST"){
  data_norm <- t(apply(feat_abd,1, function(x){x/sum(x)}))
  data_norm <- t(apply(data_norm, 1, function(x) asin(sqrt(x))))
} else if (norm_method =="CLR"){
  data_norm <- clr(feat_abd)
} else if (norm_method =="TMM"){
  data_norm = t(feat_abd)
  dge <- DGEList(counts = data_norm)
  dge <- calcNormFactors(dge,lib.size=NULL,method = "TMM")
  data_norm <- t(data_norm)/(dge$samples$norm.factors)
}


write.table(data_norm,file = paste(opts$workplace,opts$output,"_normalized_abundance.txt",sep=''),sep = '\t',col.names = NA)

print("FINISH")

