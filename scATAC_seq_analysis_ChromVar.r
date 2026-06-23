##scATAC-data processing
library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(harmony)
setwd("~/Documents/Milan/Epigenomics/data/scatac")

fragment_files <- list.files(
  path = "Cancer_scATACseq_data/flat_fragments",
  pattern = "\\.fragments\\.tsv\\.gz$",
  full.names = TRUE
)

print(length(fragment_files)) 
print(fragment_files)          




peaks_df <- read.csv("Cancer_scATACseq_data/flat_fragments/GBM_peakset.csv")
peaks <- makeGRangesFromDataFrame(peaks_df, keep.extra.columns = TRUE,
                                  seqnames.field = "seqnames",
                                  start.field = "start",
                                  end.field = "end")

for (i in seq_along(fragment_files)) {
  fragment_path <- fragment_files[i]
  
  sample_name <- gsub(".*scATAC_(GBMx_.*)\\.fragments.*", "\\1", fragment_path)
  message("Processing sample: ", sample_name)
  
  
  frag_obj <- CreateFragmentObject(path = fragment_path)
  
  
  chrom_assay <- CreateChromatinAssay(
    counts = FeatureMatrix(
      fragments = frag_obj,
      features = peaks,
      cells = NULL
    ),
    sep = c(":", "-"),
    fragments = frag_obj,
    min.cells = 0,
    min.features = 0
  )
  
  
  seurat_obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    project = sample_name
  )
  

  saveRDS(seurat_obj, file = paste0("Cancer_scATACseq_data/flat_fragments/", sample_name, "_seurat.rds"))
  
  
  seurat_list[[sample_name]] <- seurat_obj
}


########################################################################################################

rm(list = ls())
library(Signac)
library(Seurat)
library(GenomicRanges)
#library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(harmony)
library(AnnotationHub)
# Set working directory
setwd("~/Documents/Milan/Epigenomics/data/scatac")




rds_files <- list.files(
  path = "Cancer_scATACseq_data/flat_fragments/",
  pattern = "_seurat\\.rds$",
  full.names = TRUE
)


seurat_list <- lapply(rds_files, readRDS)

names(seurat_list) <- gsub(".*/(GBMx_.*)_seurat\\.rds$", "\\1", rds_files)


ah <- AnnotationHub()
ensdb_v98 <- ah[["AH75011"]]

annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

for (i in seq_along(seurat_list)) {
  Annotation(seurat_list[[i]]) <- annotations
}


#QC (TSSEnrichment, NucleosomeSignal)
for (i in seq_along(seurat_list)) {
  message("Running QC for: ", names(seurat_list)[i])
  seurat_list[[i]] <- TSSEnrichment(seurat_list[[i]])
  seurat_list[[i]] <- NucleosomeSignal(seurat_list[[i]])
}


output_dir <- "processed_seurat_objects/"
dir.create(output_dir, showWarnings = FALSE)

for (i in seq_along(seurat_list)) {
  sample_id <- names(seurat_list)[i]
  saveRDS(
    seurat_list[[i]],
    file = file.path(output_dir, paste0(sample_id, "_no_renamed.rds"))
  )
}


filtered_list <- list()

# Loop through each sample
for (i in seq_along(seurat_list)) {
  sample_id <- names(seurat_list)[i]
  obj <- seurat_list[[i]]
  
  cat("\n=== Sample:", sample_id, "===\n")
  
  # Total cells before any filtering
  total_cells <- ncol(obj)
  cat("Total cells:", total_cells, "\n")
  
 
  obj <- subset(obj, subset = !is.na(nucleosome_signal))
  cat("After removing NA nucleosome_signal:", ncol(obj), "\n")
  
 
  count_min <- sum(obj$nCount_peaks > 5)
  count_max <- sum(obj$nCount_peaks < 150000)
  count_nuc <- sum(obj$nucleosome_signal < 4)
  count_tss <- sum(obj$TSS.enrichment > 0)

  cat("Cells with nCount_peaks > 5:", count_min, "\n")
  cat("Cells with nCount_peaks < 150000:", count_max, "\n")
  cat("Cells with nucleosome_signal < 4:", count_nuc, "\n")
  cat("Cells with TSS.enrichment > 0:", count_tss, "\n")

 
  filtered_obj <- subset(
    obj,
    subset = nCount_peaks > 5 & nCount_peaks < 150000 &
             nucleosome_signal < 4
  )
  

  
  final_cells <- ncol(filtered_obj)
  percent_kept <- round(100 * final_cells / total_cells, 1)

  cat("Final filtered cells:", final_cells, "\n")
  cat("Percentage kept:", percent_kept, "%\n")

 
  filtered_list[[sample_id]] <- filtered_obj
  
}



output_dir <- "processed_seurat_objects/"
dir.create(output_dir, showWarnings = FALSE)

for (i in seq_along(filtered_list)) {
  sample_id <- names(filtered_list)[i]
  saveRDS(
    filtered_list[[i]],
    file = file.path(output_dir, paste0(sample_id, "_filtered1.rds"))
  )
}

################################################################################
##Main label transfer and down stream analysis

rm(list = ls())

library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(harmony)
library(AnnotationHub)
library(future)
library(cicero)
library(dplyr)

# ==============================
# Gene Activity (Sequential)
# ==============================


setwd("~/Documents/Milan/Epigenomics/data/scatac")


filtered_files <- list.files(
  path = "processed_seurat_objects/",
  pattern = "_filtered1\\.rds$",
  full.names = TRUE
)
filtered_list <- lapply(filtered_files, readRDS)
names(filtered_list) <- gsub(".*/(GBMx_.*)_filtered1\\.rds$", "\\1", filtered_files)


ah <- AnnotationHub()
ensdb_v98 <- ah[["AH75011"]]
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
genome(annotations) <- "hg38"


output_dir <- "processed_seurat_objects/gene_activity/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Compute gene activity sequentially
for (i in seq_along(filtered_list)) {
  sample_id <- names(filtered_list)[i]
  cat("\n=== Processing sample:", sample_id, "===\n")
  
  atac <- filtered_list[[i]]
  
  if (is.null(Annotation(atac))) Annotation(atac) <- annotations
  
  if (!"ACTIVITY" %in% names(atac@assays)) {
    cat("Computing gene activity for sample:", sample_id, "...\n")
    gene.activities <- GeneActivity(atac)  # uses TSS ±2kb + gene body
    atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
    atac <- NormalizeData(atac, assay = "ACTIVITY")
    atac <- ScaleData(atac, assay = "ACTIVITY")
    
    rm(gene.activities)
    gc()
  } else {
    cat("Gene activity already exists for sample:", sample_id, "\n")
  }
  
  saveRDS(atac, file = file.path(output_dir, paste0(sample_id, "_gene_activity.rds")))
  rm(atac)
  gc()
  
  cat("Sample", sample_id, "gene activity done.\n")
}

cat("\nAll samples: Gene activity completed.\n")


##Reference building
##########################################
#!/usr/bin/env Rscript

rm(list = ls())

library(Seurat)
library(dplyr)

setwd("~/Documents/Milan/Epigenomics/data/scatac")

# ==============================
# Load reference
# ==============================
load("/home/user/Documents/Milan/Epigenomics/data/scRNA_seq_GBM4set_mice2set/GSE182109/gbm_GSE182109_cellstate_subtype_cellphonedb.RData")

DefaultAssay(sce) <- "RNA"

cat("Original reference:", ncol(sce), "cells\n")

# ==============================
# State-aware downsampling
# ==============================
set.seed(123)

sce_list <- SplitObject(sce, split.by = "celltype_state")

sce_down_list <- lapply(names(sce_list), function(state) {
  
  x <- sce_list[[state]]
  n_cells <- ncol(x)
  
  is_glioma <- grepl("Glioma", state)
  
  if (is_glioma) {
    n_keep <- min(1500, n_cells)
  } else {
    n_keep <- ifelse(n_cells > 1000, 750, n_cells)
  }
  
  cat(state, ":", n_cells, "→", n_keep, "\n")
  
  subset(x, cells = sample(colnames(x), n_keep))
})

sce <- merge(sce_down_list[[1]], y = sce_down_list[-1])

cat("Final reference:", ncol(sce), "cells\n")


# ==============================
# Normalize + PCA
# ==============================
cat("Running normalization + PCA...\n")

sce <- NormalizeData(sce, verbose = FALSE)
sce <- FindVariableFeatures(sce, nfeatures = 1500, verbose = FALSE)
sce <- ScaleData(sce, verbose = FALSE)
sce <- RunPCA(sce, npcs = 25, verbose = FALSE)

# ==============================
# Precompute neighbors
# ==============================
cat("Precomputing PCA neighbors...\n")


sce <- FindNeighbors(
  sce,
  reduction = "pca",
  dims = 1:7,
  return.neighbor = TRUE,
  compute.SNN = FALSE,
  nn.method = "annoy",
  k.param=50,
  verbose = FALSE
)

names(sce@graphs)[1] <- "RNA.nn"
# ==============================
# sanity check
# ==============================
cat("Available neighbor:\n")
print(names(sce@neighbors))

# ==============================
# Save optimized reference
# ==============================
saveRDS(sce, "optimized_reference.rds")

cat("\nReference ready: optimized_reference.rds\n")


##Rscript to run 
#######################coarse.R

#!/usr/bin/env Rscript

rm(list = ls())


library(Seurat)
library(Signac)
library(dplyr)
library(future)


Sys.setenv(OMP_NUM_THREADS = 2)
Sys.setenv(MKL_NUM_THREADS = 2)

plan("multicore", workers = 2)
options(future.globals.maxSize = 8 * 1024^3)


setwd("~/Documents/Milan/Epigenomics/data/scatac")


sce <- readRDS("optimized_reference.rds")
DefaultAssay(sce) <- "RNA"

cat("Loaded reference:", ncol(sce), "cells\n")
cat("Available neighbors:", names(sce@neighbors), "\n")

if (!"RNA.nn" %in% names(sce@neighbors)) {
  stop("ERROR: RNA.nn not found in reference")
}


dims_coarse <- 1:12

activity_files <- list.files(
  "processed_seurat_objects/gene_activity/",
  pattern = "_gene_activity\\.rds$",
  full.names = TRUE
)

out_dir <- "processed_seurat_objects/label_transfer/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================
# SAFE WRAPPER
# ==============================
safe_run <- function(expr, sample_id) {
  tryCatch(
    expr,
    error = function(e) {
      cat("\n❌ ERROR in sample:", sample_id, "\n")
      cat("Message:", e$message, "\n")
      return(NULL)
    }
  )
}

# ==============================
# MAIN LOOP (STREAMING)
# ==============================
for (file in activity_files) {

  sample_id <- gsub(".*/(GBMx_.*)_gene_activity\\.rds$", "\\1", file)
  out_file <- file.path(out_dir, paste0(sample_id, "_coarse.rds"))

  
  if (file.exists(out_file)) {
    cat("\nSkipping (cached):", sample_id, "\n")
    next
  }

  cat("\n==============================\n")
  cat("COARSE:", sample_id, "\n")
  cat("==============================\n")

  # ------------------------------
  # Load ATAC (ONE AT A TIME)
  # ------------------------------
  atac <- readRDS(file)
  DefaultAssay(atac) <- "ACTIVITY"

  # ------------------------------
  # FEATURE SELECTION
  # ------------------------------
  shared_features <- intersect(
    VariableFeatures(sce),
    rownames(atac)
  )

  n_use <- min(1000, length(shared_features))
  features <- shared_features[1:n_use]

  cat("Features used:", length(features), "\n")

  if (length(features) < 400) {
    cat("Too few features, skipping\n")
    rm(atac); gc()
    next
  }


# ------------------------------
# ADAPTIVE ANCHOR FINDING
# ------------------------------
k_filter_use <- min(50, max(20, floor(length(features) / 50)))

cat("Using k.filter =", k_filter_use, "\n")

anchors <- safe_run(
  FindTransferAnchors(
    reference = sce,
    query = atac,
    reference.assay = "RNA",
    query.assay = "ACTIVITY",
    features = features,
    reduction = "pcaproject",
    reference.reduction = "pca",
    dims = dims_coarse,
    k.filter = k_filter_use,
    reference.neighbors = "RNA.nn"
  ),
  sample_id
)

if (is.null(anchors)) {
  rm(atac); gc()
  next
}


# ------------------------------
# ROBUST TRANSFER
# ------------------------------
n_anchors <- nrow(anchors@anchors)

cat("Total anchors:", n_anchors, "\n")

if (n_anchors == 0) {
  cat(" No anchors found — forcing minimal transfer\n")
  rm(atac, anchors); gc()
  next
}

# Always valid k.weight
k_use <- max(1, min(10, floor(n_anchors / 2)))

cat("Using k.weight =", k_use, "\n")

pred <- safe_run(
  TransferData(
    anchorset = anchors,
    refdata = sce$Assignment,
    dims = dims_coarse,
    k.weight = k_use
  ),
  sample_id
)

if (is.null(pred)) {
  rm(atac, anchors); gc()
  next
}


  # ------------------------------
  # ADD METADATA
  # ------------------------------
  atac <- AddMetaData(atac, pred)

  # ------------------------------
  # SAFE SAVE (atomic)
  # ------------------------------
  tmp_file <- paste0(out_file, ".tmp")
  saveRDS(atac, tmp_file)
  file.rename(tmp_file, out_file)

  cat("Finished:", sample_id, "\n")

  rm(atac, anchors, pred)
  gc()
}

cat("\n COARSE DONE (SAFE & STABLE)\n")

######################################################fine.R#####################
#!/usr/bin/env Rscript

rm(list = ls())


library(Seurat)
library(Signac)
library(dplyr)
library(future)


Sys.setenv(OMP_NUM_THREADS = 2)
Sys.setenv(MKL_NUM_THREADS = 2)

plan("multicore", workers = 2)
options(future.globals.maxSize = 8 * 1024^3)


setwd("~/Documents/Milan/Epigenomics/data/scatac")

coarse_dir <- "processed_seurat_objects/label_transfer/"
out_dir <- "processed_seurat_objects/label_transfer/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# ==============================
# Load reference 
# ==============================
sce <- readRDS("optimized_reference.rds")
DefaultAssay(sce) <- "RNA"

cat("Loaded reference:", ncol(sce), "cells\n")
cat("Neighbors available:", names(sce@neighbors), "\n")

if (is.null(sce@neighbors$RNA.nn)) {
  stop("ERROR: RNA.nn not found in reference object")
}


dims_use <- 1:7   

safe_run <- function(expr, sample_id) {
  tryCatch(expr,
           error = function(e) {
             cat("\n❌ ERROR:", sample_id, "\n", e$message, "\n")
             return(NULL)
           })
}


coarse_files <- list.files(
  coarse_dir,
  pattern = "_coarse.rds$",
  full.names = TRUE
)


for (file in coarse_files) {

  sample_id <- gsub(".*/(GBMx_.*)_coarse.rds$", "\\1", file)
  out_file <- file.path(out_dir, paste0(sample_id, "_labeled_state.rds"))

  if (file.exists(out_file)) {
    cat("\n Skipping:", sample_id, "\n")
    next
  }

  cat("\n==============================\n")
  cat("FINE:", sample_id, "\n")
  cat("==============================\n")

  atac <- readRDS(file)
  DefaultAssay(atac) <- "ACTIVITY"

  atac$celltype_state_final <- NA
  atac$state_score <- NA


  for (ct in unique(atac$predicted.id)) {

    cat("  Lineage:", ct, "\n")

 
    cells_query <- colnames(atac)[
      atac$predicted.id == ct &
      atac$prediction.score.max > 0.5
    ]

    if (length(cells_query) < 30) {
      cat("Too few query cells\n")
      next
    }


    ref_cells <- colnames(sce)[sce$Assignment == ct]

    if (length(ref_cells) < 50) {
      cat("Too few reference cells\n")
      next
    }

    # ------------------------------
    # Feature selection
    # ------------------------------
    features_all <- VariableFeatures(sce)
    shared_features <- features_all[features_all %in% rownames(atac)]

    if (length(shared_features) < 200) {
      cat("Too few shared features\n")
      next
    }

    n_use <- min(600, length(shared_features))
    features_use <- shared_features[1:n_use]

    # ------------------------------
    # STABLE REDUCTION 
    # ------------------------------
    reduction_use <- "pcaproject"
    k_filter_use <- min(30, max(10, floor(n_use / 50)))

    cat("   Features:", n_use,
        "| k.filter:", k_filter_use, "\n")


    # ------------------------------
    # Anchor finding
    # ------------------------------
    anchors <- safe_run(
      FindTransferAnchors(
        reference = sce,
        query = atac,
        reference.assay = "RNA",
        query.assay = "ACTIVITY",
        features = features_use,
        reduction = reduction_use,
        reference.reduction = "pca",
        dims = dims_use,
        k.filter = k_filter_use,
        reference.neighbors = "RNA.nn"
      ),
      sample_id
    )

    if (is.null(anchors)) next

    n_anchors <- nrow(anchors@anchors)

    if (n_anchors < 10) {
      cat("Too few anchors\n")
      next
    }

    k_use <- max(1, min(6, floor(n_anchors / 2)))

    cat("   Anchors:", n_anchors,
        "| k.weight:", k_use, "\n")

    # ------------------------------
    # Transfer 
    # ------------------------------
    pred <- safe_run(
      TransferData(
        anchorset = anchors,
        refdata = sce$celltype_state,
        dims = dims_use,
        k.weight = k_use
      ),
      sample_id
    )

    if (is.null(pred)) next


    atac$celltype_state_final[cells_query] <- pred$predicted.id
    atac$state_score[cells_query] <- pred$prediction.score.max

    rm(anchors, pred)
    gc()
  }

  # ==============================
  # Confidence filtering
  # ==============================
  atac$state_filtered <- ifelse(
    atac$state_score < 0.6,
    "low_confidence",
    atac$celltype_state_final
  )


  tmp_file <- paste0(out_file, ".tmp")
  saveRDS(atac, tmp_file)
  file.rename(tmp_file, out_file)

  rm(atac)
  gc()

  cat("Finished:", sample_id, "\n")
}

cat("\nFINE LABEL TRANSFER COMPLETED (STABLE & PRODUCTION READY)\n")




##Motif aalysis
###################################################################################scmotif.R
#!/usr/bin/env Rscript

library(Signac)
library(Seurat)
library(GenomicRanges)
library(AnnotationHub)
library(Matrix)
library(future)
library(chromVAR)
library(motifmatchr)
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(BiocParallel)


setwd("~/Documents/Milan/Epigenomics/data/scatac")

input_dir  <- "processed_seurat_objects/"
output_dir <- "processed_seurat_objects/chromvar/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

plan("sequential")

cat("\n chromVAR PIPELINE STARTED\n")

# ==============================
# MOTIFS 
# ==============================
cat("\n Loading JASPAR motifs...\n")

pfm <- getMatrixSet(
  JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

cat(" Motifs loaded:", length(pfm), "\n")

# ==============================
# FILES
# ==============================
files <- list.files(
  input_dir,
  pattern = "_filtered1\\.rds$",
  full.names = TRUE
)

cat("\nTotal samples:", length(files), "\n")

# ==============================
# LOOP
# ==============================
for (file in files) {

  sample_id <- gsub(".*/(GBMx_.*)_filtered1\\.rds$", "\\1", file)
  out_file <- file.path(output_dir, paste0(sample_id, "_chromvar.rds"))

  if (file.exists(out_file)) {
    cat("\n Skipping:", sample_id, "\n")
    next
  }

  cat("\n==============================\n")
  cat("Sample:", sample_id, "\n")
  cat("==============================\n")


  atac <- readRDS(file)
  DefaultAssay(atac) <- "peaks"

  cat("\n Loaded object\n")
  print(atac)

  peak_mat <- GetAssayData(atac, assay = "peaks", layer = "counts")

  cat("\n Raw matrix:", dim(peak_mat), "\n")

  # ------------------------------
  # REMOVE EMPTY PEAKS
  # ------------------------------
  peak_sums <- Matrix::rowSums(peak_mat)

cutoff <- quantile(peak_sums, 0.25)  # Q1

peak_mat <- peak_mat[peak_sums > cutoff, , drop = FALSE]


  #peak_mat <- peak_mat[Matrix::rowSums(peak_mat) > 0, , drop = FALSE]

  cat(" Non-zero peaks:", nrow(peak_mat), "\n")

  if (nrow(peak_mat) < 5000) {
    cat("Too few peaks, skipping\n")
    next
  }

  peak_names <- rownames(peak_mat)

  coords <- do.call(rbind, strsplit(peak_names, "-"))

  peak_ranges <- GRanges(
    seqnames = coords[,1],
    ranges = IRanges(
      start = as.numeric(coords[,2]),
      end   = as.numeric(coords[,3])
    )
  )

  genome(peak_ranges) <- "hg38"

  cat("\n GRanges created\n")
  print(head(peak_ranges))

  # ------------------------------
  # BUILD SE OBJECT
  # ------------------------------
  se <- SummarizedExperiment(
    assays = list(counts = peak_mat),
    rowRanges = peak_ranges
  )

 
  colData(se)$depth <- Matrix::colSums(assay(se, "counts"))

  cat("\n SE created:", dim(assay(se)), "\n")

  # ------------------------------
  # FILTER 
  # ------------------------------
  se <- chromVAR::addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)

  #se <- chromVAR::filterSamples(se, min_depth = 2100)
  se <- chromVAR::filterPeaks(se)

  cat("\nAfter filtering:\n")
  cat("Peaks:", nrow(se), "Cells:", ncol(se), "\n")

  if (nrow(se) < 5000) {
    cat("Too few peaks after filtering\n")
    next
  }


  cat("\n Matching motifs...\n")

  motif_ix <- motifmatchr::matchMotifs(
    pfm,
    se,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )

  cat(" Motifs matched\n")

  # ------------------------------
  # chromVAR
  # ------------------------------
  cat("\n Running chromVAR...\n")

  dev <- chromVAR::computeDeviations(se, annotations = motif_ix)
  tf_mat <- chromVAR::deviationScores(dev)

  cat(" chromVAR completed\n")
  cat("TF matrix:", dim(tf_mat), "\n")


  atac[["chromvar"]] <- CreateAssayObject(counts = tf_mat)
  atac <- NormalizeData(atac, assay = "chromvar", verbose = FALSE)


  saveRDS(atac, out_file)

  cat(" Saved:", sample_id, "\n")

  rm(atac, se, dev, tf_mat, peak_mat, peak_ranges, motif_ix)
  gc()
}

cat("\n chromVAR COMPLETED SUCCESSFULLY\n")



##Ceiero analysis
####################################################################cicero.R
#!/usr/bin/env Rscript

rm(list = ls())

library(Signac)
library(Seurat)
library(Matrix)
library(monocle)     
library(cicero)
library(GenomicRanges)
library(AnnotationHub)

setwd("~/Documents/Milan/Epigenomics/data/scatac")

input_dir  <- "processed_seurat_objects/"
output_dir <- "processed_seurat_objects/cicero/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\nCICERO PIPELINE (FINAL - YOUR LOGIC)\n")


cat("\n Loading hg38 genome info...\n")

ah <- AnnotationHub()
ensdb_v98 <- ah[["AH75011"]]

annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
annotations <- keepStandardChromosomes(annotations, pruning.mode = "coarse")

genome.df <- data.frame(
  chr = seqlevels(annotations),
  length = seqlengths(annotations)
)

genome.df <- genome.df[!is.na(genome.df$length), ]


print(head(genome.df))


files <- list.files(
  input_dir,
  pattern = "_filtered1\\.rds$",
  full.names = TRUE
)

cat("\nTotal samples:", length(files), "\n")

# ==============================
# MAIN LOOP
# ==============================
for (file in files) {

  sample_id <- gsub(".*/(GBMx_.*)_filtered1\\.rds$", "\\1", file)
  out_file <- file.path(output_dir, paste0(sample_id, "_cicero.rds"))

  if (file.exists(out_file)) {
    cat("\nSkipping:", sample_id, "\n")
    next
  }

  cat("\n==============================\n")
  cat("Sample:", sample_id, "\n")
  cat("==============================\n")


  atac <- readRDS(file)
  DefaultAssay(atac) <- "peaks"

  cat("\n Loaded object\n")
  print(atac)

  peak_mat <- GetAssayData(atac, assay = "peaks", layer = "counts")

  cat("\n Raw matrix:", dim(peak_mat), "\n")

  # ------------------------------
  # Q1 FILTER 
  # ------------------------------
  peak_sums <- Matrix::rowSums(peak_mat)
  cutoff <- quantile(peak_sums, 0.25)

  peak_mat <- peak_mat[peak_sums > cutoff, , drop = FALSE]
  atac <- atac[rownames(peak_mat), ]

  cat(" Peaks after Q1 filtering:", dim(peak_mat), "\n")

  if (nrow(peak_mat) < 5000) {
    cat("Too few peaks — skipping\n")
    next
  }


  peak_mat@x[peak_mat@x > 0] <- 1

  # ------------------------------
  # SUBSAMPLE CELLS 
  # ------------------------------
  if (ncol(peak_mat) > 20000) {
    set.seed(123)
    keep_cells <- sample(colnames(peak_mat), 20000)
    peak_mat <- peak_mat[, keep_cells]
    atac <- subset(atac, cells = keep_cells)
    cat(" Subsampled cells:", ncol(peak_mat), "\n")
  }

  cat("Cells used:", ncol(peak_mat), "\n")


  cat("\n Creating CellDataSet...\n")

  cellinfo <- data.frame(row.names = colnames(peak_mat))
  peakinfo <- data.frame(
    gene_short_name = rownames(peak_mat),
    row.names = rownames(peak_mat)
  )

  pd <- new("AnnotatedDataFrame", data = cellinfo)
  fd <- new("AnnotatedDataFrame", data = peakinfo)

  cds <- newCellDataSet(
    peak_mat,
    phenoData = pd,
    featureData = fd,
    expressionFamily = binomialff()
  )

  cat(" CDS created:", dim(exprs(cds)), "\n")

  # ------------------------------
  # LSI 
  # ------------------------------
  cat("\n Running LSI...\n")

  atac <- RunTFIDF(atac)
  atac <- FindTopFeatures(atac, min.cutoff = 'q0')
  atac <- RunSVD(atac)

  reduced_coords <- Embeddings(atac, "lsi")[colnames(peak_mat), 1:30]

  cat(" LSI done:", dim(reduced_coords), "\n")

  # ------------------------------
  #CICERO CDS
  # ------------------------------
  cat("\nCreating Cicero CDS...\n")

  cicero_cds <- make_cicero_cds(
    cds,
    reduced_coordinates = reduced_coords
  )

  cat(" Cicero CDS ready\n")


  cat("\n Running Cicero...\n")

  conns <- run_cicero(
    cicero_cds,
    genomic_coords = genome.df,
    sample_num = 100
  )

  cat("Cicero completed\n")

  # ------------------------------
  # DEBUG OUTPUT
  # ------------------------------
  cat("\n Connections preview:\n")
  print(head(conns))

  cat("\n Total connections:", nrow(conns), "\n")
  cat("\n Co-accessibility summary:\n")
  print(summary(conns$coaccess))


  saveRDS(conns, out_file)

  cat("Saved:", sample_id, "\n")

  rm(atac, cds, cicero_cds, conns, peak_mat)
  gc()
}

cat("\nCICERO COMPLETED SUCCESSFULLY\n")

#####################################cicero visualization
#!/usr/bin/env Rscript

rm(list = ls())


library(cicero)
library(GenomicRanges)
library(AnnotationHub)
library(rtracklayer)
library(Signac)
library(Seurat)
library(Matrix)
library(monocle)     
library(cicero)
library(GenomicRanges)
library(AnnotationHub)

cat("\nCICERO VISUALIZATION PIPELINE\n")


setwd("~/Documents/Milan/Epigenomics/data/scatac")


cicero_file <-"processed_seurat_objects/cicero/GBMx_09C0DCE7_D669_4D28_980D_BF71179116A4_X005_S04_B1_T1_cicero.rds"

conns <- readRDS(cicero_file)

cat("\n Loaded Cicero connections\n")
print(head(conns))


cat("\n Loading hg38 gene annotations...\n")

ah <- AnnotationHub()
ensdb_v98 <- ah[["AH75011"]]   # EnsDb v98 (hg38)

annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)

# standardize to UCSC style (chr1, chr2...)
seqlevelsStyle(annotations) <- "UCSC"
annotations <- keepStandardChromosomes(annotations, pruning.mode = "coarse")

genome(annotations) <- "hg38"

cat(" Annotation loaded\n")


cat("\n Building clean gene model...\n")
gene_anno <- data.frame(
  chromosome = as.character(seqnames(annotations)),
  start = start(annotations),
  end = end(annotations),
  gene = annotations$gene_id,
  symbol = annotations$gene_name,
  strand = as.character(strand(annotations)),
  transcript = annotations$tx_id   # <-- REQUIRED FIX
)

# replace missing transcript IDs if any
gene_anno$transcript[is.na(gene_anno$transcript)] <- "NA"

# remove duplicates (important)
gene_anno <- gene_anno[!duplicated(gene_anno$gene), ]


cat(" Clean gene model ready\n")
print(head(gene_anno))
cat("Gene model ready\n")



# ==============================
# COORDINATES (±1MB WINDOW)
# ==============================

chr3_conns <- conns[
  grepl("^chr3-", conns$Peak1) | grepl("^chr3-", conns$Peak2),
]

chr <- "chr3"
minbp <- 46089058 - 1e6
maxbp <- 46090453 + 1e6


get_midpoint <- function(x) {
  parts <- strsplit(gsub("^chr3-", "", x), "-")
  sapply(parts, function(p) mean(as.numeric(p)))
}

mid1 <- get_midpoint(chr3_conns$Peak1)
mid2 <- get_midpoint(chr3_conns$Peak2)


region_idx <- (mid1 >= minbp & mid1 <= maxbp) |
              (mid2 >= minbp & mid2 <= maxbp)

chr3_local <- chr3_conns[region_idx, ]

top_hits <- chr3_local[order(chr3_local$coaccess, decreasing = TRUE), ]

viewpoint <- top_hits$Peak1[1]
viewpoint

chr3_local$Peak2 <- as.character(chr3_local$Peak2)

cons_4c <- chr3_local[
  chr3_local$Peak1 == viewpoint |
  chr3_local$Peak2 == viewpoint,
]


tiff('ccr1_scatac.tiff', units="in", width=12, height=10, res=1200, compression = 'lzw')

plot_connections(
  conns,
  chr = "chr3",
  minbp = minbp,
  maxbp = maxbp,
  viewpoint = viewpoint,
  gene_model = gene_anno,
  coaccess_cutoff = 0,
  comparison_track = cons_4c,
  connection_width = 0.5,
  include_axis_track = FALSE
)
dev.off()

cat("\n DONE: chr3 4C plot generated\n")

#######LSI and Harmony integration and custers


#!/usr/bin/env Rscript

rm(list = ls())

library(Signac)
library(Seurat)
library(future)


plan(multisession, workers = 2)

Sys.setenv(
  OMP_NUM_THREADS = 2,
  MKL_NUM_THREADS = 2,
  OPENBLAS_NUM_THREADS = 2
)

options(future.globals.maxSize = 80 * 1024^3)

setwd("~/Documents/Milan/Epigenomics/data/scatac")

input_file  <- "Cancer_scATACseq_data/GBMx_integrated_scATAC_with_labels.rds"
output_file <- "Cancer_scATACseq_data/GBMx_GLIOMA_LSI.rds"


cat("\n Loading object...\n")
combined <- readRDS(input_file)

cat("\nFiltering Glioma...\n")

cells_use <- colnames(combined)[combined$predicted.id == "Glioma"]
combined <- subset(combined, cells = cells_use)

rm(cells_use)
gc()

cat("Glioma cells:", ncol(combined), "\n")


combined@reductions <- list()
combined@graphs <- list()

DefaultAssay(combined) <- "peaks"

# ------------------------------
# LSI
# ------------------------------
cat("\n Running LSI...\n")

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = "q0")

VariableFeatures(combined) <- head(VariableFeatures(combined), 30000)

combined <- RunSVD(combined, n = 30)

gc()


saveRDS(combined, output_file)

cat("\nDONE: LSI saved\n")

######################################################################################

#!/usr/bin/env Rscript

rm(list = ls())

library(Signac)
library(Seurat)
library(harmony)
library(future)


plan(multisession, workers = 2)

Sys.setenv(
  OMP_NUM_THREADS = 2,
  MKL_NUM_THREADS = 2,
  OPENBLAS_NUM_THREADS = 2
)

options(future.globals.maxSize = 80 * 1024^3)

setwd("~/Documents/Milan/Epigenomics/data/scatac")

input_file  <- "Cancer_scATACseq_data/GBMx_GLIOMA_LSI.rds"
output_file <- "Cancer_scATACseq_data/GBMx_GLIOMA_HARMONY.rds"


cat("\nLoading LSI object...\n")
combined <- readRDS(input_file)


# ------------------------------
# STEP 2: EXTRACT LSI 
# ------------------------------
emb <- Embeddings(combined, "lsi")
emb <- emb[, 2:15]   

harmony_embeddings <- harmony::HarmonyMatrix(
  data_mat = emb,
  meta_data = combined@meta.data,
  vars_use = "sample",
  do_pca = FALSE,
  theta = 2,
  lambda = 0.1,
  max.iter.harmony = 10,
  verbose = TRUE
)

gc()


cat("\n Storing Harmony embeddings...\n")

combined[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embeddings,
  key = "harmony_",
  assay = DefaultAssay(combined)
)

gc()

# ------------------------------
# QUICK QC
# ------------------------------
cat("\n QC check...\n")

print(dim(harmony_embeddings))

cat("\nCells per sample after balancing:\n")
print(table(combined$sample))


cat("\nSaving...\n")

saveRDS(combined, output_file)

cat("\nDONE: Harmony saved\n")

####################################################################################
#!/usr/bin/env Rscript

rm(list = ls())

library(Signac)
library(Seurat)
library(harmony)
library(future)


plan(multisession, workers = 2)

Sys.setenv(
  OMP_NUM_THREADS = 2,
  MKL_NUM_THREADS = 2,
  OPENBLAS_NUM_THREADS = 2
)

options(future.globals.maxSize = 80 * 1024^3)


setwd("~/Documents/Milan/Epigenomics/data/scatac")

input_file  <- "Cancer_scATACseq_data/GBMx_GLIOMA_HARMONY.rds"
output_file <- "Cancer_scATACseq_data/GBMx_GLIOMA_CLUSTERED.rds"


cat("\n[1] Loading Harmony object...\n")
combined <- readRDS(input_file)

cat("Cells:", ncol(combined), "\n")

# ------------------------------
# SANITY CHECK
# ------------------------------
if (!"harmony" %in% names(combined@reductions)) {
  stop("ERROR: Harmony reduction not found")
}

harmony_dims <- ncol(Embeddings(combined, "harmony"))
cat("Harmony dims:", harmony_dims, "\n")

if (harmony_dims < 10) {
  stop("ERROR: Not enough Harmony dimensions")
}

# ------------------------------
# NEIGHBORS 
# ------------------------------

cat("\nFinding neighbors...\n")

combined <- FindNeighbors(
  combined,
  reduction = "harmony",
  dims = 1:10,     
  k.param = 100     
)



gc()


cat("\nClustering...\n")


combined <- FindClusters(
  combined,
  resolution = 0.03,   
  algorithm = 4,
  random.seed = 1234
)
gc()

# ------------------------------
# QC: cluster numbers
# ------------------------------
cat("\nCluster distribution:\n")
print(table(combined$seurat_clusters))

n_clusters <- length(unique(combined$seurat_clusters))
cat("\nTotal clusters:", n_clusters, "\n")

# ------------------------------
# QC: sample mixing
# ------------------------------
if ("sample" %in% colnames(combined@meta.data)) {
  cat("\n Sample distribution per cluster:\n")
  print(table(combined$sample, combined$seurat_clusters))
} else {
  cat("\n WARNING: 'sample' column not found\n")
}




saveRDS(combined, output_file)

cat("\nDONE: Coarse clustering complete\n")
########################################################################
#!/usr/bin/env Rscript

rm(list = ls())

library(Signac)
library(Seurat)
library(harmony)
library(future)


plan(multisession, workers = 2)

Sys.setenv(
  OMP_NUM_THREADS = 2,
  MKL_NUM_THREADS = 2,
  OPENBLAS_NUM_THREADS = 2
)

options(future.globals.maxSize = 80 * 1024^3)


setwd("~/Documents/Milan/Epigenomics/data/scatac")

input_file  <- "Cancer_scATACseq_data/GBMx_GLIOMA_CLUSTERED.rds"
output_file <- "Cancer_scATACseq_data/GBMx_GLIOMA_FINAL.rds"


cat("\nLoading clustered object...\n")
combined <- readRDS(input_file)

cat("Cells:", ncol(combined), "\n")

# ------------------------------
# SAFETY CHECK
# ------------------------------
if (!"harmony" %in% names(combined@reductions)) {
  stop("ERROR: Harmony reduction not found!")
}

cat("Harmony dims:", dim(Embeddings(combined, "harmony")), "\n")

# ------------------------------
# UMAP 
# ------------------------------
cat("\nRunning UMAP...\n")

set.seed(1234)

combined <- RunUMAP(
  combined,
  reduction = "harmony",
  dims = 1:10,        
  n.neighbors = 50,   
  min.dist = 0.6,     
  spread = 1.2,
  verbose = TRUE
)

gc()

# ------------------------------
# QUICK QC
# ------------------------------
cat("\nCluster sizes:\n")
print(table(combined$seurat_clusters))



saveRDS(combined, output_file)

cat("\n DONE: Final object saved\n")

