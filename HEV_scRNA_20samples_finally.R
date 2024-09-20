#https://www.jianshu.com/p/550f62642443如何修改R包源代码及使用修改后的R包-以DoubletFinder为例
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

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
library(DoubletFinder)
library(harmony)
library(patchwork)

setwd("~/nanjing/HEV_scRNA-seq_nankeda")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/alf01")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc,  subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
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
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc,  subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
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
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "alf03.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/alf04")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "alf04.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/alf05")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc,  subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "alf05.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/alf06")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "alf06.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/alf07")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
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
pbmc <- subset(pbmc, subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
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
pbmc <- subset(pbmc,  subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
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
pbmc <- subset(pbmc,  subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
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
pbmc <- subset(pbmc,  subset = nFeature_RNA >1276 & nFeature_RNA < 2500 & percent.mt <3)
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
pbmc <- subset(pbmc,  subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "hc04.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/ahe01")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc,  subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "ahe01.rds")

#not use ahe02,use AHE02.rda
#a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/ahe02")
#pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
#pbmc
#pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 1300 & nFeature_RNA < 2500 & percent.mt <3)
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#all.genes <- rownames(pbmc)
#pbmc <- ScaleData(pbmc, features = all.genes)
#pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=48)
#pbmc <- RunUMAP(pbmc, dims = 1:20)
#saveRDS(pbmc, file = "ahe02.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/ahe03")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "ahe03.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/ahe04")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "ahe04.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/ahe05")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc,  subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "ahe05.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/ahe06")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc,  subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "ahe06.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/ahe07")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "ahe07.rds")

a1 <- Read10X(data.dir = "~/nanjing/HEV_scRNA-seq_nankeda/ahe08")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 1276 & nFeature_RNA < 2500 & percent.mt <3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "ahe08.rds")

AHE01<-readRDS("ahe01.rds")
#AHE02<-readRDS("ahe02.rds")
AHE03<-readRDS("ahe03.rds")
AHE04<-readRDS("ahe04.rds")
AHE05<-readRDS("ahe05.rds")
AHE06<-readRDS("ahe06.rds")
AHE07<-readRDS("ahe07.rds")
AHE08<-readRDS("ahe08.rds")
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
AHE02<-readRDS("AHE02.rds")

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

HC02<-RenameCells(HC02,add.cell.id="HC02",for.merge=T)
HC02@meta.data$tech<-"HC"
HC02@meta.data$celltype<-"HC02"

HC03<-RenameCells(HC03,add.cell.id="HC03",for.merge=T)
HC03@meta.data$tech<-"HC"
HC03@meta.data$celltype<-"HC03"

HC04<-RenameCells(HC04,add.cell.id="HC04",for.merge=T)
HC04@meta.data$tech<-"HC"
HC04@meta.data$celltype<-"HC04"

hev<-merge(ALF01, y = c(ALF02, ALF03, ALF04, ALF05, ALF06, ALF07, ALF08,
                        AHE01, AHE02, AHE03, AHE04, AHE05, AHE06, AHE07, AHE08,
                        HC01, HC02, HC03, HC04), 
              add.cell.ids = c("ALF01", "ALF02", "ALF03", "ALF04", "ALF05", "ALF06", "ALF07", "ALF08",
                               "AHE01","AHE02", "AHE03", "AHE04", "AHE05", "AHE06", "AHE07", "AHE08",
  "HC01", "HC02", "HC03", "HC04"), project = "hev")

saveRDS(hev, file="hev.rds")
hms<-hev

hms<-readRDS(file="hev.rds")
#before integrate
hms[["percent.mt"]] <- PercentageFeatureSet(hms, pattern = "^Mt-")
VlnPlot(hms, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = hms, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) 
 
plot_grid(p1,p2)

saveRDS(pancreas, file="pancreas_1276hev.rds")
pancreas<-readRDS(file="pancreas_1276hev.rds")
DimPlot(pancreas, reduction = "umap",label = TRUE, repel = TRUE)
#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}
reference.list <- pancreas.list[c("ALF01", "ALF02", "ALF03", "ALF04", "ALF05", "ALF06", "ALF07", "ALF08",
                                  "AHE01","AHE02", "AHE03", "AHE04", "AHE05", "AHE06", "AHE07", "AHE08",
                                  "HC01", "HC02", "HC03", "HC04")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
saveRDS(pancreas.anchors, file="pancreas.anchors_1276hev.rds")

pancreas.anchors<-readRDS(file="pancreas.anchors_1276hev.rds")
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
saveRDS(pancreas.integrated, file="pancreas.integrated_1276hev.rds")





pancreas.integrated<-readRDS(file="pancreas.integrated_1276hev.rds")
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_1276integrated.rds")

hms_individual_integrated<-readRDS(file="hms_after_1276integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype",label=TRUE,pt.size=0.001)
p1
p2 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype",label=FALSE,pt.size=0.001)
p2
#pbmc <- JackStraw(hms_individual_integrated, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:30)
#JackStrawPlot(pbmc, dims = 1:30)
#ElbowPlot(pbmc,ndims = 30)
#find how many 15cluster
#ElbowPlot(hms_individual_integrated)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,15,by=1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 15)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE)
saveRDS(hms_cluster, file = "hms_cluster_test_15.rds")

hms_cluster<-readRDS(file="hms_cluster_test_15.rds")
scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.0001, 
                                logfc.threshold = 0.00001
)
write.table(scRNA.markers,file="cellMarkers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) 
write.csv(file="top20_cell_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top100_marker_genes_15.csv",top20_table,row.names=F)

hms_cluster<-readRDS("hms_cluster_test_15.rds")
DimPlot(hms_cluster, reduction = "umap",label = TRUE)

#细胞及细胞中基因与RNA数量
hms_cluster<-pancreas.integrated
slotNames(hms_cluster)
#assay
hms_cluster@assays
dim(hms_cluster@meta.data)
View(hms_cluster@meta.data)

new.cluster.ids <- c("cd4_naive","cd4_naive","cd4_naive","cd4_naive","cd4_tm","cd4_tm","cd4_tm","cd4_naive","cd8_naive","cd8_effector","cd4_naive","cd4_naive","b","cd4_naive","nk","cd4_naive","cd8_effector","cd8_effector","cd8_effector","cd4_naive","cd4_tm","cd8_naive","cd4_naive","nk","b","cd4_naive","cd8_effector","cd8_effector","cd4_naive","cd4_tm","cd8_effector","nk","b","nk","cd4_naive","cd4_naive","cd8_naive","cd8_effector","monocyte","monocyte","nk","cd8_effector","monocyte","nk","cd8_effector","cd8_cytotoxic","monocyte","b","cd8_effector","dc","cd4_naive","cd4_tm","cd4_naive","b","monocyte","cd8_naive","monocyte","cd8_naive","cd4_naive","dc","b","ilc","monocyte","nk","cd4_naive","nk","treg","cd4_naive","nk","nk","cd8_cytotoxic","red_blood","cd4_naive","nk","cd8_effector","nk","monocyte","nk","cd8_cytotoxic","treg","b","cd8_cytotoxic","cd4_naive","monocyte","treg","megakaryocyte","cd4_effector","b","megakaryocyte","cd8_effector","cd4_naive","b","cd4_naive","dc","nk","cd4_helper","cd4_naive","cd8_effector","megakaryocyte","cd8_effector","cd8_effector","cd8_effector","dc","plasma","plasma","megakaryocyte","cd4_naive","cd4_naive","red_blood","megakaryocyte","b"
)
names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
hms_cluster_id[['cell_type']]<-hms_cluster_id@active.ident
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.01) 
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.01) 
saveRDS(hms_cluster_id, file = "20samples_id_test.rds")

a<-hms_cluster_id@meta.data
write.csv(a,"a")


features = c("MZB1","LILRA4","PLD4", "SERPINF1","CD74","TRAC","AQP3","GZMK", "TRGC1","LYAR","PPBP","TUBB1","PF4", 
            "HBA1","HBA2","FOXP3",
              "FCGR3A","IFITM3","CST3", "LGALS1","GZMH","CCL5","CST7","CMC1","LYZ",
              "S100A8","S100A9","CD7","CTSW","FGFBP2","PRF1","GNLY","NKG7","CD79A","BANK1",          
             "CD8A","CD8B","CCR7","TCF7","EEF1A1","PASK","TSHZ2","IL7R"
         )
DotPlot(hms_cluster_id,features=features)+RotatedAxis()
