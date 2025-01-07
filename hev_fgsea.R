if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("fgsea")

setwd("~/nanjing/HEV_scRNA-seq_nankeda/deg")
library(tibble)
library(fgsea)			# GSEA分析主程序
library(data.table)	# 数据处理
library(ggplot2)		# 画图处理
library(dplyr)			# 数据处理
#library(msigdb)			# 包含基因集合，通常和GSEA分析共同使用
#library(GSEABase)		# 可以提供GSEA基础结构和函数,也会被其他包调用

#hcq_vs_con
res<-read.table("tmp05")#norm_foldchange
ranks <- deframe(res)
pathways<-gmtPathways("c2.all.v2024.1.Hs.symbols.gmt")
fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 100000)
fwrite(fgseaRes, file="nk_alf_ahe_gsea.txt", sep="\t", sep2=c("", " ", ""))

plotEnrichment(pathways.kegg[["KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS"]],
               ranks) + labs(title="SYSTEMIC_LUPUS_ERYTHEMATOSUS")


#pcq_vs_con
res<-read.table("pcq_vs_con_fgsea.txt")
ranks <- deframe(res)
fgseaRes_kegg<- fgsea(pathways = pathways.kegg, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 100000)
fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 100000)
fwrite(fgseaRes_kegg, file="pcq_vs_con_kegg01.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="pcq_vs_con_all01.txt", sep="\t", sep2=c("", " ", ""))


plotEnrichment(pathways.kegg[["KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS"]],
               ranks) + labs(title="SYSTEMIC_LUPUS_ERYTHEMATOSUS")

plotEnrichment(pathways.kegg[["KEGG_LEISHMANIA_INFECTION"]],
               ranks) + labs(title="LEISHMANIA_INFECTION")


#phpma_vs_con
res<-read.table("phpma_vs_con_fgsea.txt")
ranks <- deframe(res)
fgseaRes_kegg<- fgsea(pathways = pathways.kegg, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 100000)
fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 100000)
fwrite(fgseaRes_kegg, file="phpma_vs_con_kegg01.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="phpma_vs_con_all01.txt", sep="\t", sep2=c("", " ", ""))

