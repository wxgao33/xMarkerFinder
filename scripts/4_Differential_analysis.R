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
    make_option(c("-p", "--profile"), type="character", default="filtered_abundance.txt",
                help="filtered reletive abundance profile[default %default]"),
    make_option(c("-g", "--group"), type="character", default="Group",
                help="the column name of experimental interest(group) in metadata [default %default]"),
    make_option(c("-e", "--exposure"), type="character", default="CRC",
                help="the experiment group (exposure) of interest in metadata [default %default]"),
    make_option(c("-b", "--batch"), type="character", default="Cohort",
                help="major confounder [default %default]"),
    make_option(c("-c", "--covariates"), type="character", default="covariates.txt",
                help="minor confounders [default %default]"),
    make_option(c("-a", "--adjust"), type="character", default="F",
                help="the option to use adjusted pvalues instead of pvalues[default %default]"),
    make_option(c("-t", "--threshold"), type="character", default="0.05",
                help="pvalue threshold [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output file prefix [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))

  # 显示输入输出确认是否正确
  print(paste("The workplace is ", opts$workplace,  sep = ""))
  print(paste("The metadata file is ", opts$metadata,  sep = ""))
  print(paste("The profile of features is ", opts$profile,  sep = ""))
  print(paste("The group is ", opts$group,  sep = ""))
  print(paste("The group of interest is ", opts$exposure,  sep = ""))
  print(paste("The batch is ", opts$batch,  sep = ""))
  print(paste("The covariates file is ", opts$covariates,  sep = ""))
  print(paste("The adjust option for pvalues is ", opts$adjust,  sep = ""))
  print(paste("The pvalue threshold is ", opts$threshold, sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}


package_list = c("MMUPHin","dplyr","ggplot2")
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
feat_abd <- read.csv(file = paste(opts$workplace,opts$profile,sep=''),sep = '\t',header =  TRUE,row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
feat_abd = t(feat_abd)
feat_abd[is.na(feat_abd)]<-0
batch <- opts$batch
exposure <- opts$group
case <- opts$exposure
adjust <- opts$adjust
p.threshold <- as.numeric(opts$threshold)
options(warn=-1)

metadata <- metadata %>% 
    mutate(Group = case_when(
        Group == case ~ paste0("z",case),
        TRUE ~ paste0("a",Group)
    ))

metadata <- metadata[colnames(feat_abd),]

scale_0_to_1 <- function(data) {
  t(apply(data, 1, function(row) {
    min_val <- min(row, na.rm = TRUE)
    max_val <- max(row, na.rm = TRUE)
    return((row - min_val) / (max_val - min_val))
  }))
}

feat_abd_scale <- scale_0_to_1(feat_abd)

#Differential analysis
file_path <- paste(opts$workplace,opts$covariates,sep='')
file_details <- file.info(file_path)
lines <- readLines(file_path, warn = FALSE)

if (file_details$size > 0 && length(lines) > 0) {
    covariates <- row.names(read.table(file_path,row.names = 1))
    fit_meta <- lm_meta(feature_abd = feat_abd_scale,
                        exposure = exposure,
                        batch = batch,
                        covariates = covariates,
                        data = metadata)
} else {
    fit_meta <- lm_meta(feature_abd = feat_abd_scale,
                        exposure = exposure,
                        batch = batch,
                        #covariates = covariates,
                        data = metadata)

}


meta_results <- fit_meta$meta_fits
maaslin_results <- fit_meta$maaslin_fits

write.table(meta_results,file = paste(opts$workplace,opts$output,"_differential_significance.txt",sep=''),sep = '\t',col.names = NA)
write.table(maaslin_results,file = paste(opts$workplace,opts$output,"_differential_significance_single_cohort.txt",sep=''),sep = '\t',col.names = NA)

if (adjust == "F"){
    feat_diff <- t(feat_abd[(meta_results %>% filter(pval<p.threshold))[,'feature'],])
    diff.plot <- meta_results[,c('feature','coef','pval')]
    diff.plot$label = ifelse(diff.plot$pval<p.threshold, ifelse(diff.plot$coef>0,paste0("Upregulated in ", case),paste0("Downregulated in ", case)),"Not significant")
    p <- ggplot(
      diff.plot, aes(x = coef, y = -log10(pval), colour=label)) +
      geom_point(alpha=0.4, size=3.5) +
      scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
      geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
      geom_hline(yintercept = -log10(p.threshold),lty=4,col="black",lwd=0.8) +
      labs(x="Coef",
           y="-log10 (p-value)")+
      theme_bw()+
      theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  
            legend.position="right", 
            legend.title = element_blank())+
      ggtitle("Volcano plot of all features")

} else if (adjust == "bonf"){
    feat_diff <- t(feat_abd[(meta_results %>% filter(pval.bonf<p.threshold))[,'feature'],])
    diff.plot <- meta_results[,c('feature','coef','pval.bonf')]
    diff.plot$label = ifelse(diff.plot$pval.bonf<p.threshold, ifelse(diff.plot$coef>0,paste0("Upregulated in ", case),paste0("Downregulated in ", case)),"Not significant")
    p <- ggplot(
      diff.plot, aes(x = coef, y = -log10(pval.bonf), colour=label)) +
      geom_point(alpha=0.4, size=3.5) +
      scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
      geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
      geom_hline(yintercept = -log10(p.threshold),lty=4,col="black",lwd=0.8) +
      labs(x="Coef",
           y="-log10 (pval.bonf)")+
      theme_bw()+
      theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
            legend.position="right", 
            legend.title = element_blank())+
      ggtitle("Volcano plot of all features")

} else if (adjust == "fdr"){
    feat_diff <- t(feat_abd[(meta_results %>% filter(qval.fdr<p.threshold))[,'feature'],])
    diff.plot <- meta_results[,c('feature','coef','qval.fdr')]
    diff.plot$label = ifelse(diff.plot$qval.fdr<p.threshold, ifelse(diff.plot$coef>0,paste0("Upregulated in ", case),paste0("Downregulated in ", case)),"Not significant")
    p <- ggplot(
      diff.plot, aes(x = coef, y = -log10(qval.fdr), colour=label)) +
      geom_point(alpha=0.4, size=3.5) +
      scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
      geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
      geom_hline(yintercept = -log10(p.threshold),lty=4,col="black",lwd=0.8) +
      labs(x="Coef",
           y="-log10 (qval.fdr)")+
      theme_bw()+
      theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
            legend.position="right", 
            legend.title = element_blank())+
      ggtitle("Volcano plot of all features")
}

write.table(feat_diff,file = paste(opts$workplace,opts$output,"_differential_signature.txt",sep=''),sep = '\t',col.names = NA)

#plot
pdf(file = paste(opts$workplace,opts$output,"_differential_plot.pdf",sep=''))
p
while (!is.null(dev.list()))  dev.off()

svg(file = paste(opts$workplace,opts$output,"_differential_plot.svg",sep=''))
p
dev.off()

print("FINISH")

