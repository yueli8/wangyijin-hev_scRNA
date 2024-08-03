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

setwd("~/nanjing/HEV_scRNA-seq_nankeda")

#input each cluster
AHE01<-readRDS("AHE01.rds")
AHE02<-readRDS("AHE02.rds")
AHE03<-readRDS("AHE03.rds")
AHE04<-readRDS("AHE04.rds")
AHE05<-readRDS("AHE05.rds")
AHE06<-readRDS("AHE06.rds")
AHE07<-readRDS("AHE07.rds")
AHE08<-readRDS("AHE08.rds")
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

ahe_hc<-merge(AHE01, y = c(AHE02, AHE03, AHE04, AHE05, AHE06, AHE07, AHE08, HC01, HC02, HC03, HC04), 
              add.cell.ids = c("AHE01","AHE02", "AHE03", "AHE04", "AHE05", "AHE06", "AHE07", "AHE08", "HC01", "HC02",
                               "HC03", "HC04"), project = "ahe_hc")

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

saveRDS(pancreas, file="pancreas_ahe.rds")
pancreas<-readRDS(file="pancreas_ahe.rds")
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
saveRDS(pancreas.anchors, file="pancreas.anchors_ahe.rds")
pancreas.anchors<-readRDS(file="pancreas.anchors_ahe.rds")

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
saveRDS(pancreas.integrated, file="pancreas.integrated_ahe.rds")

pancreas.integrated<-readRDS(file="pancreas.integrated_ahe.rds")
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_integrated_ahe.rds")

hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,4,by=0.1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 4)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "hms_cluster_test_ahe_4.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.05, 
                                logfc.threshold = 0.05
)
write.table(scRNA.markers,file="cellMarkers_ahe.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="top20_cell_markers_ahe.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top20_marker_genes_ahe_4.csv",top20_table,row.names=F)

#细胞及细胞中基因与RNA数量
slotNames(hms_cluster)
#assay
hms_cluster@assays
dim(hms_cluster@meta.data)
View(hms_cluster@meta.data)

hms_cluster<-readRDS(file="hms_cluster_test_ahe_4.rds")
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
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.1) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.1) 
saveRDS(hms_cluster_id, file = "hms_cluster_id_test_ahe_4.rds")



#alf04.rds,alf06.rds used company provided
setwd("~/nanjing/HEV_scRNA-seq_nankeda")
a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/alf01")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset =  percent.mt <4)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "alf01.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/alf02")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset =  percent.mt <4)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "alf02.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/alf03")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = percent.mt <4)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "alf03.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/alf05")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = percent.mt <4)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "alf05.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/alf07")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset =  percent.mt <4)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "alf07.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/alf08")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset =percent.mt <4)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "alf08.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/hc01")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = percent.mt <4)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "hc01.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/hc02")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = percent.mt <4)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "hc02.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/hc03")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc,subset = percent.mt <4)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "hc03.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/hc04")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = percent.mt <4)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "hc04.rds")

#input each cluster
ALF01<-readRDS("alf01.rds")
ALF02<-readRDS("alf02.rds")
ALF03<-readRDS("alf03.rds")
ALF04<-readRDS("alf04.rds")
ALF05<-readRDS("alf05.rds")
ALF06<-readRDS("alf06.rds")
ALF07<-readRDS("alf07.rds")
ALF08<-readRDS("alf08.rds")
HC01<-readRDS("hc01.rds")
HC02<-readRDS("hc02.rds")
HC03<-readRDS("hc03.rds")
HC04<-readRDS("hc04.rds")

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

HC01<-RenameCells(HC01,add.cell.id="HC01",for.merge=T)
HC01@meta.data$tech<-"HC"
HC01@meta.data$celltype<-"HC01"

HC02<-RenameCells(HC02,add.cell.id="HC02",for.merge=T)
HC02@meta.data$tech<-"HC"
HC02@meta.data$celltype<-"HC02"

HC03<-RenameCells(HC03,add.cell.id="HC03",for.merge=T)
HC03@meta.data$tech<-"HC"
HC03@meta.data$celltype<-"HC03"

HC04<-RenameCells(HC04,add.cell.id="HC04",for.merge=T)
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

saveRDS(pancreas, file="pancreas_alf.rds")
pancreas<-readRDS(file="pancreas_alf.rds")
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
saveRDS(pancreas.anchors, file="pancreas.anchors_alf.rds")
pancreas.anchors<-readRDS(file="pancreas.anchors_alf.rds")

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
saveRDS(pancreas.integrated, file="pancreas.integrated_alf.rds")

pancreas.integrated<-readRDS(file="pancreas.integrated_alf.rds")
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_integrated.rds")

hms_individual_integrated<-readRDS(file="hms_after_integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
p1

#pbmc <- JackStraw(hms_individual_integrated, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:30)
#JackStrawPlot(pbmc, dims = 1:30)
#ElbowPlot(pbmc,ndims = 30)
#find how many 15cluster
#ElbowPlot(hms_individual_integrated)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,5,by=0.1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 5)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "hms_cluster_test_5.rds")

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
write.csv(file="top20_marker_genes_5.csv",top20_table,row.names=F)


