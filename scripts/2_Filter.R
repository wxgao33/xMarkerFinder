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
    make_option(c("-b", "--batch"), type="character", default="Cohort",
                help="the column name of batch(cohort) in metadata [default %default]"),
    make_option(c("-t", "--threshold"), type="character", default="2",
                help="minimal cohort number for feature filtering [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output file prefix [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))

  # 显示输入输出确认是否正确
  print(paste("The workplace is ", opts$workplace,  sep = ""))
  print(paste("The metadata file is ", opts$metadata,  sep = ""))
  print(paste("The profile of features is ", opts$feature_profile,  sep = ""))
  print(paste("The batch is ", opts$batch,  sep = ""))
  print(paste("The threshold is ", opts$threshold,  sep = ""))
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
metadata <- read.table(file = paste(opts$workplace,opts$metadata,sep=''),sep = '\t',header =  TRUE,row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
feat_abd <- read.csv(file = paste(opts$workplace,opts$feature_profile,sep=''),sep = '\t',header =  TRUE,row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
batch <- opts$batch
feat_abd[is.na(feat_abd)]<-0
threshold <-opts$threshold  
#rela <- t(apply(feat_abd,1, function(x){x/sum(x)}))
#write.table(rela,file = paste(opts$workplace,opts$output,"_feature_relative_abundance.txt",sep=''),sep = '\t')
feat_abd = t(feat_abd)

filter.f <- function(dat, Num){
  SD <- apply(dat,1,sd)
  num_0 <- apply(dat, 1, function(x) length(which(x == 0)))
  ave_abun <- apply(dat,1,mean)
  tmp <- cbind(dat,data.frame(SD),data.frame(num_0),data.frame(ave_abun))
  colnames(tmp)[(ncol(dat)+1):(ncol(dat)+3)] <- c("sd","count0",'avebun')
  #dat_filter <- tmp %>% filter(count0 <= as.numeric(Num*0.9) & sd >0) 
  dat_filter <- tmp[(tmp$count0 <= as.numeric(Num*0.8)) & (tmp$sd >0) & (tmp$avebun > 0.0001),]
  #dat_filter <- tmp[(tmp$count0 <= as.numeric(Num*0.8)) & (tmp$sd >0),]
  return(dat_filter[,1:ncol(dat)])
}

temp.mean = function(feat_abd,metadata){
    t(sapply(row.names(feat_abd),FUN=function(feature){
        sapply(unique(metadata[batch]),FUN=function(feature,cohort){
            mean.ab=mean(feat_abd[feature,which(metadata[batch]==cohort)])
            },feature=feature)
        }))
}

temp.mean.data <- temp.mean(feat_abd,metadata)

f.idx.data = rowSums(temp.mean.data >= 0.00001) >= threshold & row.names(feat_abd) != '-1'

feat_filter = feat_abd[f.idx.data,]

feat.filter <- t(filter.f(feat_filter,ncol(feat_filter)))

cat('Retaining', ncol(feat.filter), 'features after quality control...\n')
write.table(feat.filter,file = paste(opts$workplace,opts$output,"_filtered_feature_abundance.txt",sep=''),sep = '\t')

print("FINISH")

