library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(dplyr)
library(Matrix)
library(clustree)

setwd("~/HEV_scRNA-seq_nankeda")

#input each cluster
AHE01<-readRDS("AHE01.rds")
AHE02<-readRDS("AHE02.rds")
AHE03<-readRDS("AHE03.rds")
AHE04<-readRDS("AHE04.rds")
AHE05<-readRDS("AHE05.rds")
AHE06<-readRDS("AHE06.rds")
AHE07<-readRDS("AHE07.rds")
AHE08<-readRDS("AHE08.rds")
ALF01<-readRDS("ALF01.rds")
ALF02<-readRDS("ALF02.rds")
ALF03<-readRDS("ALF03.rds")
ALF04<-readRDS("ALF04.rds")
ALF05<-readRDS("ALF05.rds")
ALF06<-readRDS("ALF06.rds")
ALF07<-readRDS("ALF07.rds")
ALF08<-readRDS("ALF08.rds")
HC01<-readRDS("HC01.rds")
HC02<-readRDS("HC02.rds")
HC03<-readRDS("HC03.rds")
HC04<-readRDS("HC04.rds")

AHE01=UpdateSeuratObject(object=AHE01)
AHE02=UpdateSeuratObject(object=AHE02)
AHE03=UpdateSeuratObject(object=AHE03)
AHE04=UpdateSeuratObject(object=AHE04)
AHE05=UpdateSeuratObject(object=AHE05)
AHE06=UpdateSeuratObject(object=AHE06)
AHE07=UpdateSeuratObject(object=AHE07)
AHE08=UpdateSeuratObject(object=AHE08)
ALF01=UpdateSeuratObject(object=ALF01)
ALF02=UpdateSeuratObject(object=ALF02)
ALF03=UpdateSeuratObject(object=ALF03)
ALF04=UpdateSeuratObject(object=ALF04)
ALF05=UpdateSeuratObject(object=ALF05)
ALF06=UpdateSeuratObject(object=ALF06)
ALF07=UpdateSeuratObject(object=ALF07)
ALF08=UpdateSeuratObject(object=ALF08)
HC01=UpdateSeuratObject(object=HC01)
HC02=UpdateSeuratObject(object=HC02)
HC03=UpdateSeuratObject(object=HC03)
HC04=UpdateSeuratObject(object=HC04)

AHE01<-RenameCells(AHE01,add.cell.id="AHE01",for.merge=T)
AHE01@meta.data$tech<-"AHE"
AHE01@meta.data$celltype<-"AHE01"

AHE02<-RenameCells(AHE02,add.cell.id="AHE02",for.merge=T)
AHE02@meta.data$tech<-"AHE"
AHE02@meta.data$celltype<-"AHE02"

AHE03<-RenameCells(AHE03,add.cell.id="AHE03",for.merge=T)
AHE03@meta.data$tech<-"AHE"
AHE03@meta.data$celltype<-"AHE03"

AHE04<-RenameCells(AHE04,add.cell.id="AHE04",for.merge=T)
AHE04@meta.data$tech<-"AHE"
AHE04@meta.data$celltype<-"AHE04"

AHE05<-RenameCells(AHE05,add.cell.id="AHE05",for.merge=T)
AHE05@meta.data$tech<-"AHE"
AHE05@meta.data$celltype<-"AHE05"

AHE06<-RenameCells(AHE06,add.cell.id="AHE06",for.merge=T)
AHE06@meta.data$tech<-"AHE"
AHE06@meta.data$celltype<-"AHE06"

AHE07<-RenameCells(AHE07,add.cell.id="AHE07",for.merge=T)
AHE07@meta.data$tech<-"AHE"
AHE07@meta.data$celltype<-"AHE07"

AHE08<-RenameCells(AHE08,add.cell.id="AHE08",for.merge=T)
AHE08@meta.data$tech<-"AHE"
AHE08@meta.data$celltype<-"AHE08"

HC01<-RenameCells(HC01,add.cell.id="HC01",for.merge=T)
HC01@meta.data$tech<-"HC"
HC01@meta.data$celltype<-"HC01"

HC02<-RenameCells(HC01,add.cell.id="HC02",for.merge=T)
HC02@meta.data$tech<-"HC"
HC02@meta.data$celltype<-"HC02"

HC03<-RenameCells(HC01,add.cell.id="HC03",for.merge=T)
HC03@meta.data$tech<-"HC"
HC03@meta.data$celltype<-"HC03"

HC04<-RenameCells(HC01,add.cell.id="HC04",for.merge=T)
HC04@meta.data$tech<-"HC"
HC04@meta.data$celltype<-"HC04"

ALF01<-RenameCells(ALF01,add.cell.id="ALF01",for.merge=T)
ALF01@meta.data$tech<-"ALF"
ALF01@meta.data$celltype<-"ALF01"

ALF02<-RenameCells(ALF02,add.cell.id="ALF02",for.merge=T)
ALF02@meta.data$tech<-"ALF"
ALF02@meta.data$celltype<-"ALF02"

ALF03<-RenameCells(ALF03,add.cell.id="ALF03",for.merge=T)
ALF03@meta.data$tech<-"ALF"
ALF03@meta.data$celltype<-"ALF03"

ALF04<-RenameCells(ALF04,add.cell.id="ALF04",for.merge=T)
ALF04@meta.data$tech<-"ALF"
ALF04@meta.data$celltype<-"ALF04"

ALF05<-RenameCells(ALF05,add.cell.id="ALF05",for.merge=T)
ALF05@meta.data$tech<-"ALF"
ALF05@meta.data$celltype<-"ALF05"

ALF06<-RenameCells(ALF06,add.cell.id="ALF06",for.merge=T)
ALF06@meta.data$tech<-"ALF"
ALF06@meta.data$celltype<-"ALF06"

ALF07<-RenameCells(ALF07,add.cell.id="ALF07",for.merge=T)
ALF07@meta.data$tech<-"ALF"
ALF07@meta.data$celltype<-"ALF07"

ALF08<-RenameCells(ALF08,add.cell.id="ALF08",for.merge=T)
ALF08@meta.data$tech<-"ALF"
ALF08@meta.data$celltype<-"ALF08"

ahe_hc<-merge(AHE01, y = c(AHE02, AHE03, AHE04, AHE05, AHE06, AHE07, AHE08, HC01, HC02, HC03, HC04), 
              add.cell.ids = c("AHE01","AHE02", "AHE03", "AHE04", "AHE05", "AHE06", "AHE07", "AHE08", "HC01", "HC02",
                               "HC03", "HC04"), project = "ahe_hc")

#hev<- merge(AHE01, y = c(AHE02, AHE03, AHE04, AHE05, AHE06, AHE07, AHE08, ALF01, ALF02, ALF03, ALF04, ALF05,
#                          ALF06, ALF07, ALF08, HC01, HC02, HC03, HC04), 
#              add.cell.ids = c("AHE01","AHE02", "AHE03", "AHE04", "AHE05", "AHE06", "AHE07", "AHE08", "ALF01",
#                               "ALF02", "ALF03", "ALF04", "ALF05", "ALF06", "ALF07", "ALF08", "HC01", "HC02",
#                               "HC03", "HC04"), project = "hev")

saveRDS(ahe_hc, file="ahe_hc.rds")
hms<-ahe_hc

hms<-readRDS(file="ahe_hc.rds")
#before integrate
hms[["percent.mt"]] <- PercentageFeatureSet(hms, pattern = "^Mt-")
VlnPlot(hms, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = hms, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
  NoLegend()
plot_grid(p1,p2)

saveRDS(pancreas, file="pancreas_ahe_hc.rds")
pancreas<-readRDS(file="pancreas_ahe_hc.rds")
DimPlot(pancreas, reduction = "umap",label = TRUE, repel = TRUE)
#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}
reference.list <- pancreas.list[c("AHE01","AHE02", "AHE03", "AHE04", "AHE05", "AHE06", "AHE07", "AHE08", "HC01", "HC02",
                                  "HC03", "HC04")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
saveRDS(pancreas.anchors, file="pancreas.anchors_ahe_hc.rds")
pancreas.anchors<-readRDS(file="pancreas.anchors_ahe_hc.rds")

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
saveRDS(pancreas.integrated, file="pancreas.integrated_ahe_hc.rds")

pancreas.integrated<-readRDS(file="pancreas.integrated_ahe_hc.rds")
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_integrated.rds")

#pbmc <- JackStraw(hms_individual_integrated, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:30)
#JackStrawPlot(pbmc, dims = 1:30)
#ElbowPlot(pbmc,ndims = 30)
#find how many 15cluster
#ElbowPlot(hms_individual_integrated)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,4,by=0.1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 4)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "hms_cluster_test_4.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.05, 
                                logfc.threshold = 0.05
)
write.table(scRNA.markers,file="cellMarkers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="top20_cell_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top20_marker_genes_4.csv",top20_table,row.names=F)


alf_hc<-merge(ALF01, y = c(ALF02, ALF03, ALF04, ALF05, ALF06, ALF07, ALF08, HC01, HC02, HC03, HC04), 
              add.cell.ids = c("ALF01","ALF02", "ALF03", "ALF04", "ALF05", "ALF06", "ALF07", "ALF08", "HC01", "HC02",
                               "HC03", "HC04"), project = "alf_hc")
saveRDS(alf_hc, file="alf_hc.rds")
hms<-alf_hc

hms<-readRDS(file="alf_hc.rds")
#before integrate
hms[["percent.mt"]] <- PercentageFeatureSet(hms, pattern = "^Mt-")
VlnPlot(hms, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = hms, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
  NoLegend()
plot_grid(p1,p2)

saveRDS(pancreas, file="pancreas_alf_hc.rds")
pancreas<-readRDS(file="pancreas_alf_hc.rds")
DimPlot(pancreas, reduction = "umap",label = TRUE, repel = TRUE)
#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}
reference.list <- pancreas.list[c("ALF01","ALF02", "ALF03", "ALF04", "ALF05", "ALF06", "ALF07", "ALF08", "HC01", "HC02",
                                  "HC03", "HC04")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
saveRDS(pancreas.anchors, file="pancreas.anchors_alf_hc.rds")

pancreas.anchors<-readRDS(file="pancreas.anchors_alf_hc.rds")

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
saveRDS(pancreas.integrated, file="pancreas.integrated_alf_hc.rds")

pancreas.integrated<-readRDS(file="pancreas.integrated_alf_hc.rds")

DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)w
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_integrated01.rds")

hms_individual_integrated<-readRDS(file="hms_after_integrated01.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
p1

#pbmc <- JackStraw(hms_individual_integrated, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:30)
#JackStrawPlot(pbmc, dims = 1:30)
#ElbowPlot(pbmc,ndims = 30)
#find how many 15cluster
#ElbowPlot(hms_individual_integrated)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,4,by=0.1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 4)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "hms_cluster_test_401.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.05, 
                                logfc.threshold = 0.05
)
write.table(scRNA.markers,file="cellMarkers01.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="top20_cell_markers01.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top20_marker_genes_401.csv",top20_table,row.names=F)
      


hc<-merge(HC01, y = c(HC02, HC03, HC04), 
              add.cell.ids = c( "HC01", "HC02","HC03", "HC04"), project = "hc")
saveRDS(hc, file="hc.rds")
hms<-hc

hms<-readRDS(file="hc.rds")
#before integrate
hms[["percent.mt"]] <- PercentageFeatureSet(hms, pattern = "^Mt-")
VlnPlot(hms, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = hms, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
  NoLegend()
plot_grid(p1,p2)

saveRDS(pancreas, file="pancreas_hc.rds")
pancreas<-readRDS(file="pancreas_hc.rds")
DimPlot(pancreas, reduction = "umap",label = TRUE, repel = TRUE)
#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}
reference.list <- pancreas.list[c( "HC01", "HC02","HC03", "HC04")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
saveRDS(pancreas.anchors, file="pancreas.anchors_hc.rds")

pancreas.anchors<-readRDS(file="pancreas.anchors_hc.rds")

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
saveRDS(pancreas.integrated, file="pancreas.integrated_hc.rds")

pancreas.integrated<-readRDS(file="pancreas.integrated_hc.rds")















pancreas.integrated<-readRDS(file="pancreas.integrated_alf_hc.rds")

DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)w
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_integrated01.rds")

hms_individual_integrated<-readRDS(file="hms_after_integrated01.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
p1

#pbmc <- JackStraw(hms_individual_integrated, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:30)
#JackStrawPlot(pbmc, dims = 1:30)
#ElbowPlot(pbmc,ndims = 30)
#find how many 15cluster
#ElbowPlot(hms_individual_integrated)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,4,by=0.1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 4)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "hms_cluster_test_401.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.05, 
                                logfc.threshold = 0.05
)
write.table(scRNA.markers,file="cellMarkers01.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="top20_cell_markers01.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top20_marker_genes_401.csv",top20_table,row.names=F)






















#细胞及细胞中基因与RNA数量
slotNames(hms_cluster)
#assay
hms_cluster@assays
dim(hms_cluster@meta.data)
View(hms_cluster@meta.data)

hms_cluster<-readRDS(file="hms_cluster_test_4.rds")
DimPlot(hms_cluster, reduction = "umap",label=TRUE)

new.cluster.ids <- c("cd8_naive","nk","cd8_naive","cd8_mait","nk","b","b","cd8_exhausted","cd8_exhausted","cd8_exhausted","nk",
"nk","nk","monocyte","treg","macrophage","treg","monocyte","nk","cd8_em","b",
"ec","dc","cd8_em","macrophage","monocyte","b","monocyte","monocyte","cd4_naive ","dc",
"monocyte","platelet","treg","ec","mast","monocyte","cd8_em","cd8_mait","b","cd8_em",
"ec","monocyte","cd8_em","cd8_pre-exhausted","monocyte","ec","platelet","dc","ec","γδ_T",
"cd4_em","γδ_T","b","monocyte","ec","macrophage","dc","monocyte","b","nk",
"mast","dc","cd8_pre-exhausted","cd8_pre-exhausted","ec","macrophage","b") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "hms_cluster_id_test_.rds")

hms_cluster_id<-readRDS("hms_cluster_id_test.rds")

T<-subset(hms_cluster_id, idents=c('T'))
DimPlot(T, reduction = "umap")
saveRDS(T, file="T.rds")

Epithelial<-subset(hms_cluster_id, idents=c('Epithelial'))
DimPlot(Epithelial, reduction = "umap")
saveRDS(Epithelial, file="Epithelial.rds")

B<-subset(hms_cluster_id, idents=c('B'))
DimPlot(B, reduction = "umap")
saveRDS(B, file="B.rds")

Secretory<-subset(hms_cluster_id, idents=c('Secretory'))
DimPlot(Secretory, reduction = "umap")
saveRDS(Secretory, file="Secretory.rds")

Macrophage<-subset(hms_cluster_id, idents=c('Macrophage'))
DimPlot(Macrophage, reduction = "umap")
saveRDS(Macrophage, file="Macrophage.rds")

Monocyte<-subset(hms_cluster_id, idents=c('Monocyte'))
DimPlot(Monocyte, reduction = "umap")
saveRDS(Monocyte, file="Monocyte.rds")

Natural_Killer<-subset(hms_cluster_id, idents=c('Natural_Killer'))
DimPlot(Natural_Killer, reduction = "umap")
saveRDS(Natural_Killer, file="Natural_Killer.rds")

Progenitor<-subset(hms_cluster_id, idents=c('Progenitor'))
DimPlot(Progenitor, reduction = "umap")
saveRDS(Progenitor, file="Progenitor.rds")

#input each cluster
T<-readRDS("T.rds")
Epithelial<-readRDS("Epithelial.rds")
B_cell<-readRDS("B.rds")
Secretory<-readRDS("Secretory.rds")
Macrophage<-readRDS("Macrophage.rds")
Monocyte<-readRDS("Monocyte.rds")
Natural_Killer<-readRDS("Natural_Killer.rds")
Progenitor<-readRDS("Progenitor.rds")

#DimPlot
DimPlot(T, reduction = "umap", split.by = "tech")
DimPlot(Epithelial, reduction = "umap", split.by = "tech")
DimPlot(B, reduction = "umap", split.by = "tech")
DimPlot(Secretory, reduction = "umap", split.by = "tech")
DimPlot(Macrophage, reduction = "umap", split.by = "tech")
DimPlot(Monocyte, reduction = "umap", split.by = "tech")
DimPlot(Natural_Killer, reduction = "umap", split.by = "tech")
DimPlot(Progenitor, reduction = "umap", split.by = "tech")

#regroup T cell
T<-readRDS("T.rds")
#pbmc <- JackStraw(T, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:30)
#JackStrawPlot(pbmc, dims = 1:30)
#ElbowPlot(pbmc,ndims = 30)
#find how many 15cluster
hms_neighbor<- FindNeighbors(pbmc, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 1.2)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "T_test_1.2.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.05, 
                                logfc.threshold = 0.05)
write.table(scRNA.markers,file="TcellMarkers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="top20_cell_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="Ttop20_marker_genes_1.2.csv",top20_table,row.names=F)

hms_cluster<-readRDS("T_test_1.2.rds")
DimPlot(hms_cluster, reduction = "umap")

#细胞及细胞中基因与RNA数量
slotNames(hms_cluster)
#assay
hms_cluster@assays
dim(hms_cluster@meta.data)
View(hms_cluster@meta.data)

new.cluster.ids <- c("Naive_CD4_T", "Natural_Killer_T","Natural_Killer_T","T_Helper",
                     "Effector_Memory_CD4_T","Resident_Memory_CD8_T","Effector_Memory_CD4_T", 
                     "Terminally_Exhausted_CD8_T","Regulatory_T","Effector_Memory_CD4_T",
                     "Regulatory_T","Terminally_Exhausted_CD8_T","B","Macrophage","Regulatory_T",
                     "Cytotoxic_CD8_T","Effector_Memory_CD8_T", "Macrophage",
                     "Cycling_Natural_Killer", "Terminally_Exhausted_CD8_T",
                     "Pre_Exhausted_CD8_T","Epithelial")
names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "T_cluster_id_test.rds")

T_cluster_id<-readRDS("T_cluster_id_test.rds")

#only T cells
Only_T<-subset(hms_cluster_id, idents=c("Naive_CD4_T", "Natural_Killer_T","T_Helper",
                                        "Effector_Memory_CD4_T","Resident_Memory_CD8_T",
                                        "Terminally_Exhausted_CD8_T","Regulatory_T",
                                        "Cytotoxic_CD8_T","Effector_Memory_CD8_T",
                                        "Pre_Exhausted_CD8_T"))
DimPlot(Only_T, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(Only_T, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(Only_T, file = "Only_T_cluster_id_test.rds")



















#Plot
RidgePlot(CD8_T, features = c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"),
          cols = c("green3","orangered"), group.by="tech", ncol = 4) + theme(axis.title.y = element_blank())

markers.to.plot<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2")
DoHeatmap(subset(hms_cluster_id,downsample=50000),features = markers.to.plot,size=5)

genes11_heatmap<-DotPlot(hms_cluster_id,features = c("CD8A","CD8B","CD38","CD69","ENTPD1",
                                                     "GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"))+RotatedAxis()
genes11_heatmap<-genes11_heatmap$data
write.csv(genes11_heatmap,"genes11_heatmap")

#expression level in each patients
CD8_T<-readRDS("CD8_T.rds")
a<-DoHeatmap(subset(CD8_T,downsample=50000),features = markers.to.plot,size=5,group.by="celltype")
a1<-a$data
write.table(a1,"a1")

Only_T <- readRDS(file="Only_T_cluster_id_test.rds") 

DimPlot(Only_T, reduction = "umap", label = TRUE, pt.size = 0.5) 
#DimPlot(Only_T, reduction = "umap", label = FALSE, pt.size = 0.5) 

#doheatmap
markers.to.plot<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2")
p1<-DoHeatmap(subset(Only_T,downsample=50000),features = markers.to.plot,size=5)
p1
a<-p1$data
write.table(a,"a")

markers.to.plot<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2")
p1<-DoHeatmap(subset(Only_T,downsample=50000),features = markers.to.plot,size=5,group.by="tech")
p1
a<-p1$data
write.table(a,"a")


genes11_heatmap<-DotPlot(Only_T,features = c("CD8A","CD8B","CD38","CD69","ENTPD1",
                                             "GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"))+RotatedAxis()

#Plot
RidgePlot(Only_T, features = c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"),
          cols = c("green3","orangered"), group.by="tech", ncol = 4) + theme(axis.title.y = element_blank())

markers.to.plot<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2")
DoHeatmap(subset(Only_T,downsample=50000),features = markers.to.plot,size=5)

genes11_heatmap<-DotPlot(hms_cluster_id,features = c("CD8A","CD8B","CD38","CD69","ENTPD1",
                                                     "GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"))+RotatedAxis()
genes11_heatmap<-genes11_heatmap$data
write.csv(genes11_heatmap,"genes11_heatmap")

#expression level in each patients
CD8_T<-readRDS("CD8_T.rds")

a1<-a$data
write.table(a1,"a1")

#Dotplot
setwd("~/gse162498/data_finally")
hms_cluster<-readRDS("hms_cluster_id_test.rds")
DimPlot(hms_cluster, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster, reduction = "umap", label = FALSE, pt.size = 0.5) 
features = c("NNMT","ANXA1","IL7R","KLRB1","NKG7","AIF1","S100A4","STMN1",
             "FTL", "APOE","APOC1","JCHAIN", "CD79A", 
             "SCGB1A1","SFTPC", "BATF","CXCL13",
             "FOXP3", "CCR7","IL32","GZMB")
DotPlot(hms_cluster,features=features)+RotatedAxis()

hms_cluster<-readRDS("Only_T_cluster_id_test.rds")
DimPlot(hms_cluster, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster, reduction = "umap", label = FALSE, pt.size = 0.5) 
features = c("CCR7","SELL","NFKBIZ","SLC2A3","CXCL13","KLRB1","IL7R","ANXA1",
             "GZMK", "CCL4", "GZMB","CD8A","LAG3","BATF",
             "FOXP3","NKG7","PRF1","GNLY","CD8B","CD79A",
             "IFI6","ISG15","IFITM1")
DotPlot(hms_cluster,features=features)+RotatedAxis()





P57<-merge(P57_B,P57_T)
P58<-merge(P58_B,P58_T)
P60<-merge(P60_B,P60_T)
P61<-merge(P61_B,P61_T)
P5758<-merge(P57,P58)
P6061<-merge(P60,P61)
P6061_J <- merge(P60_J, P61_J)
P57586061<-merge(P5758,P6061)
P57586061J<-merge(P57586061,P6061_J)

saveRDS(P57586061J, file="P57586061J_before_integrate.rds")
saveRDS(P57586061, file="P57586061_before_integrate.rds")
hms01<-P57586061J
hms<-P57586061

#before integrate
hms[["percent.mt"]] <- PercentageFeatureSet(hms, pattern = "^Mt-")
VlnPlot(hms, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = hms, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
  NoLegend()
plot_grid(p1,p2)
#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}
#reference.list <- pancreas.list[c("Blood_P57","Tumor_P57","Blood_P58","Tumor_P58",
   #                               "Blood_P60","Tumor_P60", "Blood_P61", "Tumor_P61","Juxta_P60","Juxta_P61")]
reference.list <- pancreas.list[c("Blood_P57","Tumor_P57","Blood_P58","Tumor_P58",
                                  "Blood_P60","Tumor_P60", "Blood_P61", "Tumor_P61")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_integrated.rds")
saveRDS(pancreas.integrated, file = "P57586061_after_integrated.rds")

hms_individual_integrated<-readRDS(file="P57586061_after_integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
p1
#find how many 15cluster
ElbowPlot(hms_individual_integrated)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:10)
hms_cluster <- FindClusters( hms_neighbor, resolution = 0.5)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:10)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "hms_cluster_test.rds")

hms_cluster<-readRDS(file="hms_cluster_test.rds")
cluster0.markers <- FindMarkers(hms_cluster, ident.1=0, min.pcr=0.25)
head(cluster0.markers, n=10)
cluster1.markers <- FindMarkers(hms_cluster, ident.1=1, min.pcr=0.25)
head(cluster1.markers, n=10)
cluster2.markers <- FindMarkers(hms_cluster, ident.1=2, min.pcr=0.25)
head(cluster2.markers, n=10)
cluster3.markers <- FindMarkers(hms_cluster, ident.1=3, min.pcr=0.25)
head(cluster3.markers, n=10)
cluster4.markers <- FindMarkers(hms_cluster, ident.1=4, min.pcr=0.25)
head(cluster4.markers, n=10)
cluster5.markers <- FindMarkers(hms_cluster, ident.1=5, min.pcr=0.25)
head(cluster5.markers, n=10)
cluster6.markers <- FindMarkers(hms_cluster, ident.1=6, min.pcr=0.25)
head(cluster6.markers, n=10)
cluster7.markers <- FindMarkers(hms_cluster, ident.1=7, min.pcr=0.25)
head(cluster7.markers, n=10)
cluster8.markers <- FindMarkers(hms_cluster, ident.1=8, min.pcr=0.25)
head(cluster8.markers, n=10)
cluster9.markers <- FindMarkers(hms_cluster, ident.1=9, min.pcr=0.25)
head(cluster9.markers, n=10)
cluster10.markers <- FindMarkers(hms_cluster, ident.1=10, min.pcr=0.25)
head(cluster10.markers, n=10)
cluster11.markers <- FindMarkers(hms_cluster, ident.1=11, min.pcr=0.25)
head(cluster11.markers, n=10)
cluster12.markers <- FindMarkers(hms_cluster, ident.1=12, min.pcr=0.25)
head(cluster12.markers, n=10)
cluster13.markers <- FindMarkers(hms_cluster, ident.1=13, min.pcr=0.25)
head(cluster13.markers, n=10)

new.cluster.ids <- c("Memory CD4+ T", "Memory CD4+ T", "Naive CD4+ T","Memory CD4+ T",
                     "Natural Killer T", "CD8+ T","B","Dendritic Cell","Astrocyte",
                     "Natural Killer T", "Dendritic Cell","B") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(hms_cluster_id, file = "hms_cluster_id_test.rds")

hms_cluster_id<-readRDS("hms_cluster_id_test.rds")
  
#hms_cluster_id<-readRDS(file="hms_cluster_id_test.rds")
Natural_Killer_T <-subset(hms_cluster_id, idents=c('Natural Killer T'))
DimPlot(Natural_Killer_T, reduction = "umap")
saveRDS(Natural_Killer_T, file="Natural_Killer_T.rds")

Naive_CD4_T<-subset(hms_cluster_id, idents=c('Naive CD4+ T'))
DimPlot(Naive_CD4_T, reduction = "umap")
saveRDS(Naive_CD4_T, file="Naive_CD4_T.rds")

B<-subset(hms_cluster_id, idents=c('B'))
DimPlot(B, reduction = "umap")
saveRDS(B, file="B.rds")

Astrocyte<-subset(hms_cluster_id, idents=c('Astrocyte'))
DimPlot(Astrocyte, reduction = "umap")
saveRDS(Astrocyte, file="Astrocyte.rds")

Dendritic_cell<-subset(hms_cluster_id, idents=c('Dendritic Cell'))
DimPlot(Dendritic_cell, reduction = "umap")
saveRDS(Dendritic_cell, file="Dendritic_cell.rds")

CD8_T<-subset(hms_cluster_id, idents=c('CD8+ T'))
DimPlot(CD8_T, reduction = "umap")
saveRDS(CD8_T, file="CD8_T.rds")

Memory_CD4_T<-subset(hms_cluster_id, idents=c('Memory CD4+ T'))
DimPlot(Memory_CD4_T, reduction = "umap")
saveRDS(Memory_CD4_T, file="Memory_CD4_T.rds")

#input each cluster
Naive_CD4_T<-readRDS("Naive_CD4_T.rds")
Natural_killer_T<-readRDS("Natural_Killer_T.rds")
Astrocyte<-readRDS("Astrocyte.rds")
B_cell<-readRDS("B.rds")
Memory_CD4_T<-readRDS("Memory_CD4_T.rds")
Dendritic_cell<-readRDS("Dendritic_cell.rds")
CD8_T<-readRDS("CD8_T.rds")
hms_cluster_id<-readRDS("hms_cluster_id_test.rds")

#DimPlot
DimPlot(Naive_CD4_T, reduction = "umap", split.by = "tech")
DimPlot(Natural_killer_T, reduction = "umap", split.by = "tech")
DimPlot(Astrocyte, reduction = "umap", split.by = "tech")
DimPlot(B_cell, reduction = "umap", split.by = "tech")
DimPlot(Memory_CD4_T, reduction = "umap", split.by = "tech")
DimPlot(Dendritic_cell, reduction = "umap", split.by = "tech")
DimPlot(CD8_T, reduction = "umap", split.by = "tech")

RidgePlot(CD8_T, features = c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"),cols = c("green3","orangered"), group.by="tech", ncol = 4) + theme(axis.title.y = element_blank())

markers.to.plot<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2")
DoHeatmap(subset(hms_cluster_id,downsample=50000),features = markers.to.plot,size=5,label = FALSE)

genes11_heatmap<-DotPlot(hms_cluster_id,features = c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"))+RotatedAxis()
genes11_heatmap<-genes11_heatmap$data
write.csv(genes11_heatmap,"genes11_heatmap")

#expression level in each patients
CD8_T<-readRDS("CD8_T.rds")
a<-DoHeatmap(subset(CD8_T,downsample=50000),features = markers.to.plot,size=5,group.by="celltype",label=FALSE)
a1<-a$data
write.table(a1,"a1")





