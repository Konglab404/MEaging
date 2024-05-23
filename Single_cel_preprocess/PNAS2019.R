setwd("../PNAS19/")
library(Seurat)
library(Matrix)
library(Matrix.utils)
library(RColorBrewer)
library(dplyr)
library(data.table)

counts<-fread("01.UMI.txt.gz")
counts<-as.data.frame(counts)
rownames(counts)<-counts$V1
counts<-counts[,2:ncol(counts)]
metainfo<-read.table("03.Cell.Barcodes.txt")
sampleinfo<-read.table("Sampleinfo.txt",header = T)
pnas<-CreateSeuratObject(counts = counts,project = "Japan2019",min.cells = 1, min.features = 0)

metainfo<-merge(metainfo,sampleinfo,by.x="V2",by.y="ID")
rownames(metainfo)<-metainfo$V1
metainfo<-metainfo[rownames(pnas2019@meta.data),]
all(rownames(pnas2019@meta.data) == metainfo$V1)

pnas2019<-subset(pnas,cells=metainfo$V1)
pnas2019<-AddMetaData(pnas2019,metainfo,col.name = c("ID","barcode","AgeGroup","Age","Gender","Group"))

genes <- read.table("./genes.tsv")
mt.ensemble <- genes[grep("^MT-", genes$V2),"V1"]
pnas.percent.mito <- colSums(expm1(pnas2019@assays$RNA@data[mt.ensemble, ])) / colSums(expm1(pnas2019@assays$RNA@data))
pnas2019 <- AddMetaData(pnas2019, pnas.percent.mito, "percent.mito")

##
pnas2019.seurat<-NormalizeData(pnas2019,
                               normalization.method = "LogNormalize", scale.factor = 10000)

##
## Do a preliminary clustring for QC-plots
##
## Calculata variable genes
pnas2019.seurat <- FindVariableFeatures(pnas2019.seurat, x.low.cutoff = 0, y.cutoff = 1, do.plot = F,selection.method = "mvp")
## Regress out % mitochondrial genes and nFeature_RNA
pnas2019.seurat <- ScaleData(pnas2019.seurat, features  = VariableFeatures(pnas2019.seurat),
                             vars.to.regress = c("percent.mito","nFeature_RNA"))

saveRDS(pnas2019.seurat,file="PNAS2019_processed_seurat.rds")
