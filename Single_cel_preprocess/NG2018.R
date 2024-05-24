###########################################################################################################################
setwd("/data1/wanght/Single_cell/NG2018")
library(Seurat)
library(Matrix)
library(Matrix.utils)
library(RColorBrewer)
library(dplyr)

###########################################################################################################################
#
# Functions
#
###########################################################################################################################

extract_field <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}


# Name: get.violin.data
# Function: returns the expression data for the createn of violin plots
# Input:
#   Name      Type          Description
#   seurat    Seurat        Seurat object containing normalized expression data
#   genes     vector        Vector containing the gene names for plotting
# Output:
#   Melted dataframe containing al the info for plotting
get.violin.data <- function(seurat, genes) {
  output <- data.frame(gene = character(0), value= numeric(0), ident = character(0))
  for (gene in genes) {
    data.use = data.frame(FetchData(seurat,gene))
    data.use = t(data.use)
    data.melt=data.frame(rep(gene, length(seurat@ident)))
    colnames(data.melt)[1]="gene"
    data.melt$value=as.numeric(data.use[1,1:length(seurat@ident)])
    data.melt$id=names(data.use)[1:length(seurat@ident)]
    data.melt$ident=seurat@ident
    noise <- rnorm(length(data.melt$value))/100000
    data.melt$value=as.numeric(as.character(data.melt$value))+noise
    output <- rbind(output, data.melt)
  }
  return(output)
}

# Name: add.data
# Function: Merges cellranger count matrices into one sparse matrix
# Input:
#   Name        Type            Description
#   orig.data   Sparse matrix   Sparse matrix where the data needs to be added to.
#   lane        Integer         Number of the lane which needs to be added
# Output:
#   Merged data
add.data <- function(orig.data = NULL, lane) {
  data <- readMM(paste0("./count_matrices_per_lane/lane_", lane, "/matrix.mtx"))
  rownames(data) <- make.names(sapply(readLines(paste0("./count_matrices_per_lane/lane_", lane, "/genes.tsv")), Extr, 1, delim = "\\t"), unique=TRUE)
  colnames(data) <- paste0(sapply(readLines(paste0("./count_matrices_per_lane/lane_", lane, "/barcodes.tsv")), extract_field, 1, delim = "-"), "_lane", lane)
  
  if (is.null(orig.data)) return(data)
  else return(merge.Matrix(orig.data, data, by.x=rownames(orig.data), by.y=rownames(data), all.x=T, all.y=T))
}

# Name: DoClusterAnalysis
# Function: Run all the Seurat function for a clustering
# Input:
#   Name        Type            Description
#   seurat      Seurat object   Seurat object to do the analysis on.
#   pc.use      Integer         Number of principal componentes to include in the clustering
#   y.cutoff    Numeric         Cut-off for the variable genes function
#   ssn.res     Numeric         Resolution used for the SNN clustering
# Output:
#   The Seurat object
DoClusterAnalysis <- function(seurat, pc.use=16, y.cutoff = 1, ssn.res = 1.2) {
  seurat <- MeanVarPlot(seurat, x.low.cutoff = 0, y.cutoff = y.cutoff, do.plot = F)
  seurat <- RegressOut(seurat, genes.regress = seurat@var.genes, latent.vars = c("percent.mito", "nFeature_RNA"))
  seurat <- PCAFast(seurat, pc.genes = seurat@var.genes, pcs.compute = pc.use, do.print = F)
  PCElbowPlot(seurat, num.pc = pc.use)
  seurat <- RunTSNE(seurat, dims.use = 1:pc.use, do.fast = T)
  seurat <- FindClusters(seurat, pc.use = 1:pc.use, resolution = ssn.res, save.SNN = T, do.sparse = T)
  seurat <- BuildClusterTree(seurat, do.reorder = T, reorder.numeric = T)
  FeaturePlot(seurat, "llkDoublet.llkSinglet")
  return(seurat)
}

###########################################################################################################################
#
# Main
#
###########################################################################################################################
pilot3.merged <- add.data(lane = 1)
pilot3.merged <- add.data(pilot3.merged, lane = 2)
pilot3.merged <- add.data(pilot3.merged, lane = 3)
pilot3.merged <- add.data(pilot3.merged, lane = 4)
pilot3.merged <- add.data(pilot3.merged, lane = 5)
pilot3.merged <- add.data(pilot3.merged, lane = 6)
pilot3.merged <- add.data(pilot3.merged, lane = 7)
pilot3.merged <- add.data(pilot3.merged, lane = 8)

dim(pilot3.merged)

##
## Combined clustering
##
combined.seurat <- CreateSeuratObject(pilot3.merged, 
                                      min.cells = 1, min.genes = 0, 
                                      project = "pilot3", do.scale = F, 
                                      do.center = F, names.field = 1, names.delim = "\\-")
dim(combined.seurat)

# ## Load deAnonymized data
# pilot3.de.anonymize <- data.frame(read.table("clustering/deanonymize.txt"))
# row.names(pilot3.de.anonymize) <- pilot3.de.anonymize[,1]
# pilot3.de.anonymize[,1] <- NULL
# colnames(pilot3.de.anonymize) <- c("llkDoublet.llkSinglet", "nSNPs_tested", "sample", "lane")
# combined.seurat <- AddMetaData(combined.seurat, pilot3.de.anonymize)

## Mitochondrial genes
genes <- read.table("./count_matrices_per_lane/lane_1/genes.tsv")
mt.ensemble <- genes[grep("^MT-", genes$V2),"V1"]
combined.percent.mito <- colSums(expm1(combined.seurat@assays$RNA@data[mt.ensemble, ])) / colSums(expm1(combined.seurat@assays$RNA@data))
combined.seurat <- AddMetaData(combined.seurat, combined.percent.mito, "percent.mito")

##
combined.seurat<-NormalizeData(combined.seurat, 
                               normalization.method = "LogNormalize", scale.factor = 10000)

##
## Do a preliminary clustring for QC-plots
##
## Calculata variable genes
combined.seurat <- FindVariableFeatures(combined.seurat, x.low.cutoff = 0, y.cutoff = 1, do.plot = F,selection.method = "mvp")
## Regress out % mitochondrial genes and nFeature_RNA
combined.seurat <- ScaleData(combined.seurat, features  = VariableFeatures(combined.seurat), 
                             vars.to.regress = c("percent.mito","nFeature_RNA"))
#combined.seurat <- ScaleData(combined.seurat, genes.use = combined.seurat@var.genes) # scale all data if wanted
combined.seurat <- RunPCA(combined.seurat, features = VariableFeatures(combined.seurat), npcs = 40)
## Use the first 16 PCs for t-SNE
combined.seurat <- RunTSNE(combined.seurat, dims = 1:16, do.fast = T)

combined.meta.data <- FetchData(combined.seurat, c("nFeature_RNA", "percent.mito"))
#true.doublet <- combined.meta.data$llkDoublet.llkSinglet > 0 & !combined.meta.data$lane %in% c(2,3)
is.mito.cut <- combined.meta.data$percent.mito > 0.05
gene.high <- combined.meta.data$nFeature_RNA > 3500
gene.low <- combined.meta.data$nFeature_RNA < 500

qc.colors <- rep("Singlet", 28855)
qc.colors[gene.high] <- ">3500 genes"
qc.colors[is.mito.cut] <- ">5% mitochondrial"
#qc.colors[true.doublet] <- "Doublet"
qc.colors <- factor(qc.colors)
qc.colors <- factor(qc.colors, levels(qc.colors)[c(3, 2, 1)])

colors <- c("#153057","#009ddb","#e64b50")

library(ggplot2)
tsne.plot <- ggplot(as.matrix(combined.seurat@reductions$tsne@cell.embeddings), 
                    aes(x=tSNE_1,y=tSNE_2, colour=qc.colors)) +
  geom_point(size=0.8) +
  scale_color_manual(values=colors) +
  theme_minimal(base_size = 14) +
  ylab("t-SNE 1") + xlab("t-SNE 2") +
  theme(legend.title=element_blank(),
        legend.text=element_text(size=14),
        axis.line = element_line(size = 0.5)) +
  guides(colour = guide_legend(override.aes = list(size=8, alpha=1)))
tsne.plot

combined.seurat.subset <- SubsetData(combined.seurat, subset.name = "llkDoublet.llkSinglet", accept.high = -2)
dim(combined.seurat.subset@data)

# combined.meta.data <- FetchData(combined.seurat.subset, c("nUMI", "nGene", "percent.mito", "llkDoublet.llkSinglet", "lane"))
# ggplot(combined.meta.data, aes(nUMI, percent.mito)) + geom_hex(bins=100) + 
#   scale_fill_distiller(palette = "Spectral", name="Cell frequencies") + 
#   ylab("Fraction mtDNA-encoded genes") + xlab("Number of reads") + 
#   theme(axis.text=element_text(size=12), axis.title=element_text(size=18)) + 
#   geom_hline(yintercept = 0.06, colour="red")

combined.seurat.subset <- subset(combined.seurat, subset = percent.mito<0.05)
dim(combined.seurat.subset@data)
# combined.meta.data <- FetchData(combined.seurat.subset, c("nUMI", "nGene", "percent.mito", "llkDoublet.llkSinglet", "lane"))
# # Gene / UMI hexagon plot
# ggplot(combined.meta.data, aes(nUMI, nGene)) + 
#   geom_hex(bins=100) + 
#   scale_fill_distiller(palette = "Spectral", name="Cell frequencies") + 
#   ylab("Number of genes") + xlab("Number of reads") + 
#   theme(axis.text=element_text(size=12), axis.title=element_text(size=18)) + 
#   geom_hline(yintercept = 3500, colour="red")

combined.seurat.subset <- subset(combined.seurat.subset, subset= nFeature_RNA < 3500)
dim(combined.seurat.subset)

# Cluster_info
clusterinfo<-read.table("barcodes_to_cell_types.tsv",sep="\t",header = T)
combines.final<-subset(combined.seurat,cells = clusterinfo$barcode)
dim(combines.final)
all(dimnames(combines.final)[[2]] == clusterinfo$barcode)

combines.final<-AddMetaData(combines.final,clusterinfo$cell_type,"cell.type")

# Sampleinfo
sampleinfo<-read.table("pilot3_persons.tsv",header = F)
combines.final<-AddMetaData(combines.final,sampleinfo$V2,"Donor_ID")
sampledetail<-read.table("sample_detail",header = T)
rownames(sampledetail)<-sampledetail$Sample
combines.final<-AddMetaData(combines.final,sampledetail[combines.final$Donor_ID,]$Age,"Age")
combines.final<-AddMetaData(combines.final,sampledetail[combines.final$Donor_ID,]$Sex,"Sex")

saveRDS(combines.final,file="NG2018_Seurat_processed.RDS")
