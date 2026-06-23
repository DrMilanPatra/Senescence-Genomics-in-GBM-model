#!/usr/bin/env Rscript
library(Seurat)
library(data.table)
library(Matrix)
# Set path to folder 1
data1_dir <- "/home/user/Documents/Milan/Epigenomics/data/scRNA_seq_GBM4set_mice2set/GSE182109"



# Load data
data1 <- ReadMtx(mtx = file.path(data1_dir, "Raw_matrix.mtx.gz"), features = file.path(data1_dir, "Raw_genes.tsv.gz"),   cells = file.path(data1_dir, "Raw_barcodes.tsv.gz"))


# Create Seurat object
gbm<- CreateSeuratObject(counts = data1, project = "GBM")

metadata <- read.table(paste0(data1_dir, "/Meta_GBM.txt"), header = TRUE, sep = ",", row.names = 1)
length(intersect(rownames(metadata), colnames(gbm)))

common_cells <- intersect(rownames(metadata), colnames(gbm))

gbm <- subset(gbm, cells = common_cells)

metadata_sub <- metadata[colnames(gbm), , drop = FALSE]

gbm <- AddMetaData(gbm, metadata = metadata_sub)

#head(gbm@meta.data)

#saveRDS(gbm, file = "gbm_raw_GSE182109.rds")


# Split Seurat object by GSMID
gbm_list <- SplitObject(gbm, split.by = "GSMID")

# Export each sample for Scrublet
for (gsm in names(gbm_list)) {
  message("Exporting sample: ", gsm)
  dir.create(paste0("scrublet_input/", gsm), recursive = TRUE, showWarnings = FALSE)
  
  # Extract raw counts
  counts <- GetAssayData(gbm_list[[gsm]], slot = "counts")
  
  # Write matrix
  writeMM(counts, file = paste0("scrublet_input/", gsm, "/matrix.mtx"))
  
  # Write gene list
  write.table(rownames(counts), paste0("scrublet_input/", gsm, "/genes.tsv"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Write barcodes
  write.table(colnames(counts), paste0("scrublet_input/", gsm, "/barcodes.tsv"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}



########################################################################
##filter for singlet
library(Seurat)


gbm <- readRDS(file = "/home/user/Documents/Milan/Epigenomics/data/scRNA_seq_GBM4set_mice2set/GSE182109/gbm_raw_GSE182109.rds")



# Load scrublet results
scrublet_results <- read.csv("/home/user/scrublet_combined_results.csv", stringsAsFactors = FALSE)


rownames(scrublet_results) <- scrublet_results$barcode

gbm_cells <- colnames(gbm)

matched_results <- scrublet_results[gbm_cells, c("doublet_score", "predicted_doublet")]

all(rownames(matched_results) == colnames(gbm))  # should return TRUE

intersect(colnames(gbm@meta.data), colnames(matched_results))

gbm <- AddMetaData(gbm, metadata = matched_results)

##check for doublets
table(gbm@meta.data$predicted_doublet)

# Keep only singlets
gbm$predicted_doublet <- gbm$predicted_doublet == "True"
gbm_singlets <- subset(gbm, subset = predicted_doublet == FALSE)


# Save filtered object
saveRDS(gbm_singlets, file = "gbm_singlets_scrublet.rds")


##############################################
#additional fitering
library(Seurat)
library(data.table)
library(Matrix)
library(tidyverse)

se1 <- readRDS(file = "/home/user/Documents/Milan/Epigenomics/data/scRNA_seq_GBM4set_mice2set/GSE182109/gbm_singlets_scrublet.rds")

DefaultAssay(object = se1) <- "RNA"

se1$nFeature_RNA <- as.numeric(se1$nFeature_RNA)
se1$nCount_RNA <- as.numeric(se1$nCount_RNA)

# Add number of genes per UMI for each cell to metadata
se1$log10GenesPerUMI <- log10(se1$nFeature_RNA) / log10(se1$nCount_RNA)
# Compute percent mito ratio
se1$percent.mt <- PercentageFeatureSet(object = se1, pattern = "^MT-")

se1$mitoRatio <- se1@meta.data$percent.mt / 100

#Compute percent robo ration:
se1$percent.rp <- PercentageFeatureSet(object = se1, pattern = "^RP[SL]")
se1$riboRatio <-se1@meta.data$percent.rp / 100


metadata <- se1@meta.data

metadata$cells <- rownames(metadata)

metadata <- metadata %>%
        dplyr::rename(nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
					  

se1@meta.data <- metadata
                           


# Visualize the number of cell counts per sample
metadata %>% 
  	ggplot(aes(x=orig.ident,  fill=orig.ident)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells")
	

metadata %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 200) +
  	facet_wrap(~orig.ident)
	metadata %>% 
  	ggplot(aes(color=orig.ident, x=nGene, fill= orig.ident)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300)
	

metadata %>% 
  	ggplot(aes(color=orig.ident, x=mitoRatio, fill=orig.ident)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2)
	
## Visualize the distribution of ribosomal gene expression detected per cell
#metadata %>% 
  	ggplot(aes(color=orig.ident, x=riboRatio, fill=orig.ident)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2)
	
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  	ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.8)

#Cell-level filtering	
# Filter out low quality reads using selected thresholds - these will change with experiment
 se1 <- subset(x = se1, 
               subset= (nUMI >= 200) & 
               (nGene > 500) & 
               (log10GenesPerUMI > 0.8) &
	       (riboRatio < 0.5) & (mitoRatio < 0.20))






## nomalization

se1 <- NormalizeData(se1, normalization.method = "RC", scale.factor = 1000000)
cpm<-as.data.frame(t(se1[["RNA"]]@data))

cpm1 <- log2((cpm/10)+1)
se1[["RNA"]]@data <- as.sparse(t(as.matrix(cpm1)))

#Remove genes that are less expressed in any cell 
keep_feature <- rowSums(se1[["RNA"]]@data > 1) > ncol(se1[["RNA"]]@data)*0.01
se1[["RNA"]]@data <- se1[["RNA"]]@data[keep_feature, ]

y<- rownames(se1[["RNA"]]@data) 
se <- subset(x = se1, features = y)
#saveRDS(se, file = "gbm_GSE182109_processed2.rds")


##############clustering 

library(Seurat)
library(Matrix)

# Load data
seurat_merged <- readRDS("/home/user/Documents/Milan/Epigenomics/data/scRNA_seq_GBM4set_mice2set/GSE182109/processed/gbm_GSE182109_processed2.rds")

seurat_merged <- FindVariableFeatures(seurat_merged, nfeatures = 2000)

seurat_merged <- ScaleData(seurat_merged)
seurat_merged <- RunPCA(seurat_merged, npcs = 50)

library(harmony)
seurat_merged <- RunHarmony(seurat_merged, group.by.vars = "GSMID")
ElbowPlot(seurat_merged)

seurat_merged <- FindNeighbors(seurat_merged, reduction = "harmony", dims = 1:7, nn.method = "annoy", n.trees = 50, annoy.metric = "euclidean")
seurat_merged <- FindClusters(seurat_merged, algorithm = 1, resolution = 0.25)
seurat_merged <- RunUMAP(seurat_merged, reduction = "harmony", dims = 1:7)

saveRDS(seurat_merged, file = "gbm_GSE182109_processed_3.rds")



##sample exclusion for p16 genomic alteration

library(tidyverse)
library(ggrepel)
library(ggpubr)
library("ggsci")
library("ggplot2")
library(DESeq2)
library(MAST)
library(Seurat)
library("gridExtra")
library(gplots)
library(RColorBrewer)
library(harmony)

df = read.table("hgnc_symbol_ensembl99.txt",sep="\t",header = F, col.names = c("gene"))
se <- readRDS("C:/Users/mlnpa/Documents/glioblastoma_scrnaseq/gbm_GSE182109_processed_3.rds")

se <- subset(se, subset = Assignment == "Glioma")
se<- subset(x = se, Patient != "ndGBM-06")
se<- subset(x =se, features=df$gene)
saveRDS(se, file="C:/Users/mlnpa/Documents/glioblastoma_scrnaseq/gbm_GSE182109_filtered.rds")

#######signature scoring
rm(list = ls())
library(Seurat)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(gplots)
library("ggsci")
library("ggplot2")
library("gridExtra")
library(DESeq2)
library(MAST)
library(SingleCellExperiment)


se <- readRDS(file="C:/Users/mlnpa/Documents/glioblastoma_scrnaseq/gbm_GSE182109_filtered.rds")

beta <- read.table(file="genesets_milan/genesets.txt", sep="\t",  header = T)
beta$V2 = 1

y<- rownames(se[["RNA"]]@data)
se <- ScaleData(se, features = y, assay = "RNA", do.center=TRUE, do.scale=TRUE)



common_genes_ccle<- rownames(se[["RNA"]]@scale.data)
beta_gene <- vector("list", length = length(colnames(beta)))
for (m in colnames(beta)) {
  print(beta[[m]]);
  beta_gene[[m]] <- beta[[m]][is.element(beta[[m]], common_genes_ccle)]
  
}
filtered_list <- Filter(function(x) !is.null(x) && !is.na(x) && x != "", beta_gene)



for (name in names(filtered_list)) {
  
  print(name)
  score_name <- paste0(name, ".score")
  se <- AddModuleScore(se, features = list(filtered_list[[name]]), slot = 'scale.data', name = score_name, seed = 1)

}


#save(se, file="C:/Users/mlnpa/Documents/glioblastoma_scrnaseq/gbm_GSE182109_filtered_scores.RData")
#######################################################################
