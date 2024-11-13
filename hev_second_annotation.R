#https://www.jianshu.com/p/550f62642443如何修改R包源代码及使用修改后的R包-以DoubletFinder为例
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
#install.packages('devtools')
#devtools::install_github('immunogenomics/presto')

library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(dplyr)
library(Matrix)
library(clustree)
library(devtools)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(reshape2)
library(tidyverse)
library(garnett)
library(data.table)
library(harmony)
library(patchwork)
library(presto)


#setwd("~/nanjing/HEV_scRNA-seq_nankeda/data/169940_9_cells")
setwd("/home/r2409301158/data")
pancreas.integrated<-readRDS(file="samples_id_test169940.rds")
DimPlot(pancreas.integrated, reduction = "umap",label = TRUE,raster=FALSE)

#细胞及细胞中基因与RNA数量
hms_cluster01<-pancreas.integrated
slotNames(hms_cluster01)
#assay
hms_cluster01@assays
dim(hms_cluster@meta.data)
View(hms_cluster@meta.data)

cd8<-subset(hms_cluster, idents=c('cd8'))
DimPlot(cd8, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(cd8, file = "cd8.rds")

hms_neighbor<- FindNeighbors(cd8, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,10,by=1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 10)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE, raster=FALSE)
saveRDS(hms_cluster, file = "hms_cluster_test_10_cd8.rds")

hms_cluster<-readRDS(file="hms_cluster_test_10_cd8.rds")
scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.0001, 
                                logfc.threshold = 0.00001
)
write.table(scRNA.markers,file="cellMarkers10_cd8.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) 
write.csv(file="top20_cell_markers10_cd8.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top100_marker_genes_10_cd8.csv",top20_table,row.names=F)

hms_cluster<-readRDS("hms_cluster_test_10_cd8.rds")

new.cluster.ids <- c("cytotoxic","macrophage","dc","naive","naive","cytotoxic","mast","cytotoxic","mast","naive","naive",
                     "megakaryocyte","naive","chronic_activation","naive","effector","cytotoxic","chronic_activation","naive","effector","cytotoxic",
                     "hepatic","chronic_activation","mast","cytotoxic","chronic_activation","chronic_activation","cytotoxic","cytotoxic","naive","effector",
                     "cytotoxic","cytotoxic","naive","cytotoxic","chronic_activation","mesothellial","cytotoxic","effector","effector","effector",
                     "chronic_activation","effector","naive","b","naive","cytotoxic","effector","effector","naive","effector",
                     "plasma","ec","chronic_activation","proliferation","chronic_activation","plasma","effector","effector","cytotoxic","naive",
                     "b","nk","b","b","effector","macrophage","ec","chronic_activation","ec","cd4",
                     "chronic_activation","proliferation","monocyte","proliferation","cytotoxic","megakaryocyte","red_blood","effector","mast","ec",
                     "dc","ec","b","neutrophil","pre_exhausted"

                      )
names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
hms_cluster_id[['cell_type']]<-hms_cluster_id@active.ident
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.01) 
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.01)

cd8<-subset(hms_cluster_id, idents=c("naive","effector","cytotoxic","proliferation","chronic_activation","pre_exhausted"))
DimPlot(cd8, reduction = "umap", label = TRUE, pt.size = 0.01)
DimPlot(cd8, reduction = "umap", label = FALSE, pt.size = 0.01) 

saveRDS(cd8, file = "cd8_id.rds")
a<-cd8@meta.data
write.csv(a,"a")

features = c("IFIT1","IFIT3","IFI44L","TNFRSF4", "SELL","CCR7","CENPF","CXCL9","MKI67","TOX2","TOP2A","CXCL13","MED12L", "GZMK", "CCL3L3", "TRAPPC3L","PDCD1",
 "LAG3" ,"CTLA4" ,"GZMB","PRF1",  "IL7R","TCF7", "LEF1","MAL","GNLY","TYROBP","KLRF1","FCGR3A"  )
DotPlot(cd8,features=features)+RotatedAxis()

#cd4 cells
cd4<-subset(hms_cluster, idents=c('cd4'))
DimPlot(cd4, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(cd4, file = "cd4.rds")

hms_neighbor<- FindNeighbors(cd4, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,10,by=1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 10)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(hms_cluster, file = "hms_cluster_test_cd4.rds")

hms_cluster<-readRDS(file="hms_cluster_test_cd4.rds")
scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.0001, 
                                logfc.threshold = 0.00001
)
write.table(scRNA.markers,file="cellMarkers10.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) 
write.csv(file="top20_cell_markers10.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top100_marker_genes_10.csv",top20_table,row.names=F)

hms_cluster<-readRDS("hms_cluster_test_cd4.rds")
DimPlot(hms_cluster, reduction = "umap",label = TRUE,raster=FALSE)
new.cluster.ids <- c("naive","naive","naive","cd8","naive","naive","cd8","cd8","cd8","ec","effector",
                     "naive","naive","effector","chronic_activation","naive","naive","naive","naive","effector","naive",
                     "naive","naive","effector","naive","cd8","naive","effector","naive","naive","naive",
                     "naive","naive","macrophage","naive","naive","naive","effector","chronic_activation","naive","dc",
                     "effector","effector","naive","naive","effector","effector","treg","naive","naive","naive",
                     "chronic_activation","naive","treg","effector","naive","naive","treg","dc","chronic_activation","effector",
                     "effector","naive","treg","chronic_activation","treg","chronic_activation","naive","effector","macrophage","effector",
                     "monocyte","naive","treg","treg","effector","naive","effector","cd8","naive","cd8",
                     "effector","treg"
)
names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
hms_cluster_id[['cell_type']]<-hms_cluster_id@active.ident
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.01) 
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.01)

cd4<-subset(hms_cluster_id, idents=c("naive","effector","chronic_activation","treg"))
DimPlot(cd4, reduction = "umap", label = TRUE, pt.size = 0.01)
DimPlot(cd4, reduction = "umap", label = FALSE, pt.size = 0.01) 

features = c("FOXP3","GZMK","ACTG1", "TIGIT", "IL32", "CXCL13","PDCD1", "TP53INP1","ITM2A",  "CORO1B", "MT-CO1","MT-CO2","SYNE2", "SYNE1","ZEB2",
         "NEAT1",  "GNLY","MAL" ,"NPM1","TRABD2A","CD7", "IL7R", "LEF1")
DotPlot(cd4,features=features)+RotatedAxis()

saveRDS(cd4, file = "cd4_id.rds")
a<-cd4@meta.data
write.csv(a,"a")

#b cells
b<-subset(hms_cluster, idents=c('b'))
DimPlot(b, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(b, file = "b.rds")

hms_neighbor<- FindNeighbors(b, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,10,by=1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 10)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(hms_cluster, file = "hms_cluster_test_b.rds")

hms_cluster<-readRDS(file="hms_cluster_test_b.rds")
scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.0001, 
                                logfc.threshold = 0.00001
)
write.table(scRNA.markers,file="cellMarkers10.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) 
write.csv(file="top20_cell_markers10.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top100_marker_genes_b.csv",top20_table,row.names=F)

hms_cluster<-readRDS("b_test.rds")
DimPlot(hms_cluster, reduction = "umap",label = TRUE,raster=FALSE)
new.cluster.ids <- c("naive","active","active","cd8","gcb","active","naive","proliferation","naive","naive","naive","gcb","gcb","naive","naive","t","gcb","proliferation","active","naive","naive","naive","active","active","naive","active","proliferation","naive","gcb","memory","active","memory","naive","memory","gcb","naive","naive","naive","memory","plasma","naive","active","cd8","naive","active","active","plasma","cd8","proliferation","plasma","active","active","cd8","naive","proliferation","plasma","plasma","proliferation","cd8","active","plasma","t","plasma","active","gcb","plasma","plasma","naive","active","cd8","proliferation","naive","active","plasma","plasma","plasma","t","memory","naive","proliferation","proliferation"
)

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
hms_cluster_id[['cell_type']]<-hms_cluster_id@active.ident
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.01) 
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.01)

cd4<-subset(hms_cluster_id, idents=c("naive","memory","gcb","active","plasma","proliferation"))
DimPlot(cd4, reduction = "umap", label = TRUE, pt.size = 0.01)
DimPlot(cd4, reduction = "umap", label = FALSE, pt.size = 0.01) 

features = c("IGKC","IGHM",  "SSR4", "IGLC2" , "IGLC3","CD27","IGHG1","MS4A1", "SOX5","RGS2","COTL1", "PLD4",
             "TUBA1B","MKI67","TUBB1","TYMS","STMN1",
                "MT-CO1","MT-CO2",  "NEAT1","TCF7", "AIM2", "MARCKS","JUNB",
              "VPREB3","AREG","FCER2",  "SELL","TCL1A","IGHD"
      )
DotPlot(cd4,features=features)+RotatedAxis()

saveRDS(cd4, file = "b_id.rds")
a<-cd4@meta.data
write.csv(a,"a")

#nk cells
nk<-subset(hms_cluster, idents=c('nk'))
DimPlot(nk, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(nk, file = "nk.rds")

hms_neighbor<- FindNeighbors(nk, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,14,by=1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 14)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(hms_cluster, file = "hms_cluster_test_nk.rds")

hms_cluster<-readRDS(file="hms_cluster_test_nk.rds")
scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.0001, 
                                logfc.threshold = 0.00001
)
write.table(scRNA.markers,file="cellMarkers10.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) 
write.csv(file="top20_cell_markers10.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top100_marker_genes_nk.csv",top20_table,row.names=F)


hms_cluster<-readRDS("nk_test.rds")
DimPlot(hms_cluster, reduction = "umap",label = TRUE,raster=FALSE)
new.cluster.ids <- c( 
  "nkt_mature","nk_cytokine","nkt_mature","nkt_exhausted","nkt_exhausted","nk_cytokine","nk_mature","nk_cytokine","nk_mature","nkt_mature","nkt_mature","nk_actived","nkt_mature","nk_cytokine","nk_mature","nk_mature","nk_cytokine","nkt_mature","nk_mature","nk_cytokine","nkt_naive","nkt_mature","nkt_mature","nk_mature","nk_cytokine","nk_cytokine","nk_mature","nk_mature","nk_actived","nk_cytokine","nkt_mature","nk_actived","nkt_mature","nk_cytokine","nkt_mature","nk_actived","nk_mature","nk_cytokine","nk_cytokine","nkt_mature","nk_actived","nk_mature","nk_cytokine","nk_actived","nk_actived","nk_cytokine","nk_cytokine","nk_mature","nk_mature","nk_cytokine","nk_cytokine","nk_actived","nk_actived","nk_mature","nk_mature","nk_actived","nk_mature","nk_mature","nk_actived","nk_mature","nk_mature","nk_mature","nkt_mature","nk_mature","nk_actived","nk_cytokine","nkt_mature","nk_mature","nk_actived","nk_mature","nk_mature","nk_mature","nk_mature","nkt_mature","nkt_mature","nk_mature","nk_mature","nk_actived","nk_actived","nk_mature","nk_actived","nk_actived","nk_cytokine","nk_actived","nk_cytokine","nkt_exhausted","nk_actived","nk_mature","nk_actived","nkt_mature","nk_mature","nk_mature","nk_cytokine","nkt_mature","nk_cytokine","nk_mature","nk_mature","nkt_mature","nk_actived","nkt_mature","nk_actived","nkt_mature","nkt_naive","nk_cytokine","nk_cytokine","nk_mature","nk_mature","nk_actived","nk_actived","nk_mature","nk_mature","nk_mature","nk_mature","nk_mature","nk_mature","nk_actived","nk_actived","nkt_exhausted","nkt_mature"
  
  )

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
hms_cluster_id[['cell_type']]<-hms_cluster_id@active.ident
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.01) 
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.01)


#DimPlot(cd4, reduction = "umap", label = TRUE, pt.size = 0.01)
#DimPlot(cd4, reduction = "umap", label = FALSE, pt.size = 0.01) 

cd4<-hms_cluster_id

features = c("PSME2","CALR","ACTB","PFN1", "CD69",
             "IL2RB", "GZMK","CLDND1","NFKBIA","DUSP1","IGFBP4",
             "IER2",    "PTGDS","LAIR2","MYOM2", "MT-CO1","MT-CO2","SPON2",
           "KLRG1","IKZF2","CCL4L2","HLA-DMA", "PADI4","GNLY" ,"LGALS3","CCL5", "KLRC2", 
          "ZNF683","CD4",   "TRBV3-1","SUB1","LCNL1","NPDC1","TRBC2","ITM2A","ANXA5",  "CD82"
)

DotPlot(cd4,features=features)+RotatedAxis()

saveRDS(cd4, file = "nk_id.rds")
a<-cd4@meta.data
write.csv(a,"a")

#monocyte cells
monocyte<-subset(hms_cluster, idents=c('monocyte'))
DimPlot(monocyte, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(monocyte, file = "monocyte.rds")

hms_neighbor<- FindNeighbors(monocyte, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(1,8,by=1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 8)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(hms_cluster, file = "hms_cluster_test_monocyte.rds")

hms_cluster<-readRDS(file="hms_cluster_test_monocyte.rds")
scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.000001, 
                                logfc.threshold = 0.0000001
)
write.table(scRNA.markers,file="cellMarkers10.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) 
write.csv(file="top20_cell_markers10.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top50_marker_genes_monocyte.csv",top20_table,row.names=F)

hms_cluster<-readRDS("monocyte_test.rds")
DimPlot(hms_cluster, reduction = "umap",label = TRUE,raster=FALSE)
new.cluster.ids <- c( 
  "cd16","t","t","cd14","cd14_cd16","cd14","cd14","t","cd14","cd16","cd14","cd16","b","cd14_cd16","cd14_cd16","cd16","cd14","b","cd14_cd16","cd16","t","cd14","cd14","cd16","astrocyte","cd14_cd16","ec","cd14","cd14","cd16","cd14","cd14_cd16","t","t","cd14","cd14","ec","b","cd14","cd14","cd14_cd16","cd14","cd14_cd16","b","cd14_cd16","cd14_cd16","cd14","cd14","cd14","cd14","cd14","t","megakaryocyte","cd14","cd14","cd14","cd14","cd14","t","cd16","cd14","cd16","cd14","cd14","cd14_cd16","cd16","t","cd14_cd16","cd14","cd16","cd14_cd16","cd16","cd14","cd14_cd16","cd4","cd14","t","cd16","cd14","b","cd16","t","mast"
  
)

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
hms_cluster_id[['cell_type']]<-hms_cluster_id@active.ident
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.01) 
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.01)

cd4<-subset(hms_cluster_id, idents=c("cd16","cd14_cd16","cd14"))
DimPlot(cd4, reduction = "umap", label = TRUE, pt.size = 0.01)
DimPlot(cd4, reduction = "umap", label = FALSE, pt.size = 0.01) 

features = c("FCGR3A","MS4A7","NAMPT","IFITM3","PRELID1","S100A11","CDKN1C","FDPS","PTMA","ATP5J2","TCF7L2","AIF1","FTH1","HLA-DRA",
             "CD14","LYZ","TYROBP",
            "CEBPD","FCN1","CD36","SELL","S100A8",   "S100A9","S100A10","S100A12","MS4A6A","FTL",
            "FCER1G","S100A4","S100A6","AGTRAP","RBP7","VNN2",
          "CFD","MNDA","SRGN","CDC42EP3","CD300E",   "LGALS2","CPVL","TMEM176A","SERPINA1",
         "DUSP6","MARCKS","NEAT1","FGL2","TNFAIP2","FOS","PSAP","MCL1","LCP1"
)

DotPlot(cd4,features=features)+RotatedAxis()

saveRDS(cd4, file = "monocyte_id.rds")
a<-cd4@meta.data
write.csv(a,"a")



#macrophage cells
macrophage<-subset(hms_cluster, idents=c('macrophage'))
DimPlot(macrophage, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(macrophage, file = "macrophage.rds")
macrophage<-readRDS(file="macrophage.rds")
hms_neighbor<- FindNeighbors(macrophage, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(1,5,by=1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 5)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(hms_cluster, file = "macrophage_test.rds")

hms_cluster<-readRDS(file="macrophage_test.rds")
scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.000001, 
                                logfc.threshold = 0.0000001
)
write.table(scRNA.markers,file="cellMarkers10.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) 
write.csv(file="top20_cell_markers10.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top50_marker_genes_macrophage.csv",top20_table,row.names=F)

#megakaryocyte cells
megakaryocyte<-subset(hms_cluster, idents=c('megakaryocyte'))
DimPlot(megakaryocyte, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(megakaryocyte, file = "megakaryocyte.rds")
megakaryocyte<-readRDS(file="megakaryocyte.rds")
hms_neighbor<- FindNeighbors(megakaryocyte, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(1,5,by=1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 5)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(hms_cluster, file = "megakaryocyte_test.rds")

hms_cluster<-readRDS(file="megakaryocyte_test.rds")
scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.000001, 
                                logfc.threshold = 0.0000001
)
write.table(scRNA.markers,file="cellMarkers10.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) 
write.csv(file="top20_cell_markers10.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top50_marker_genes_megakaryocyte.csv",top20_table,row.names=F)

#dc cells
dc<-subset(hms_cluster01, idents=c('dc'))
DimPlot(dc, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(dc, file = "dc.rds")
dc<-readRDS(file="dc.rds")
hms_neighbor<- FindNeighbors(dc, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(1,5,by=1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 5)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(hms_cluster, file = "dc_test.rds")

hms_cluster<-readRDS(file="dc_test.rds")
scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.000001, 
                                logfc.threshold = 0.0000001
)
write.table(scRNA.markers,file="cellMarkers10.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) 
write.csv(file="top20_cell_markers10.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top50_marker_genes_dc.csv",top20_table,row.names=F)

#ec cells
ec<-subset(hms_cluster01, idents=c('ec'))
DimPlot(ec, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(ec, file = "ec.rds")
ec<-readRDS(file="ec.rds")
hms_neighbor<- FindNeighbors(ec, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(1,5,by=1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 5)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(hms_cluster, file = "ec_test.rds")

hms_cluster<-readRDS(file="ec_test.rds")
scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.000001, 
                                logfc.threshold = 0.0000001
)
write.table(scRNA.markers,file="cellMarkers10.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) 
write.csv(file="top20_cell_markers10.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top50_marker_genes_ec.csv",top20_table,row.names=F)

#red_blood
red_blood<-subset(hms_cluster01, idents=c('red_blood'))
DimPlot(red_blood, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(red_blood, file = "red_blood.rds")
red_blood<-readRDS(file="red_blood.rds")
hms_neighbor<- FindNeighbors(red_blood, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(1,5,by=1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 5)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(hms_cluster, file = "red_blood_test.rds")

hms_cluster<-readRDS(file="red_blood_test.rds")
scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.000001, 
                                logfc.threshold = 0.0000001
)
write.table(scRNA.markers,file="cellMarkers10.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) 
write.csv(file="top20_cell_markers10.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top50_marker_genes_red_blood.csv",top20_table,row.names=F)

#mast
mast<-subset(hms_cluster01, idents=c('mast'))
DimPlot(mast, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(mast, file = "mast.rds")
mast<-readRDS(file="mast.rds")
hms_neighbor<- FindNeighbors(mast, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(1,2,by=1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 2)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE,raster=FALSE)
saveRDS(hms_cluster, file = "mast_test.rds")

hms_cluster<-readRDS(file="mast_test.rds")
scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.000001, 
                                logfc.threshold = 0.0000001
)
write.table(scRNA.markers,file="cellMarkers10.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 12, wt = avg_log2FC) 
write.csv(file="top20_cell_markers10.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top50_marker_genes_mast.csv",top20_table,row.names=F)

