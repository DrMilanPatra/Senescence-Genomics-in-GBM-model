
## =========================================
## chromVAR TF ~ Signature (Linear Model)
## =========================================

library(Seurat)
library(Signac)
library(Matrix)
library(matrixStats)
library(limma)
library(dplyr)

setwd("~/Documents/Milan/Epigenomics/data/scatac")

obj <- readRDS("Cancer_scATACseq_data/GBMx_all_with_chromvar.rds")

# ==============================
#KEEP VALID chromVAR CELLS
# ==============================
cat("\n[1] Subsetting chromVAR cells...\n")

DefaultAssay(obj) <- "chromvar"

tf_counts <- GetAssayData(obj, layer = "counts")
valid_cells <- colnames(obj)[Matrix::colSums(tf_counts) != 0]

obj_sub <- subset(obj, cells = valid_cells)

cat("Cells with chromVAR:", ncol(obj_sub), "\n")

# ==============================
#CLEAN TF MATRIX
# ==============================
cat("\n[2] Cleaning TF matrix...\n")

tf_mat <- GetAssayData(obj_sub, layer = "data")
tf_mat <- as.matrix(tf_mat)

# --- remove TFs with too many NA (>50%)
tf_keep <- rowMeans(is.na(tf_mat)) < 0.5

# --- remove cells with too many NA (>50%)
cell_keep <- colMeans(is.na(tf_mat)) < 0.5

tf_mat <- tf_mat[tf_keep, cell_keep]

# --- remove zero variance TFs
tf_sd <- matrixStats::rowSds(tf_mat, na.rm = TRUE)
tf_mat <- tf_mat[tf_sd > 0, ]

# --- replace remaining NA with 0
tf_mat[is.na(tf_mat)] <- 0

cat("Final TF matrix:", dim(tf_mat), "\n")

# update Seurat object to match filtered cells
obj_sub <- subset(obj_sub, cells = colnames(tf_mat))


cat("\n[3] Defining signatures...\n")

signatures <- c(
  "p16_top150_GA_score",
  "senmayo_GA_score",
  "state5_2fc_GA_score",
  "state1_2fc_GA_score",
  "state2_2fc_GA_score"
)

signatures <- signatures[signatures %in% colnames(obj_sub@meta.data)]

cat("Using signatures:", length(signatures), "\n")

# ==============================
#LINEAR MODEL LOOP
# ==============================
cat("\n[4] Running limma models...\n")

all_results <- list()

for (sig in signatures) {

  cat("\n==============================\n")
  cat("Processing:", sig, "\n")

  score <- obj_sub[[sig]][,1]

  
  valid <- !is.na(score)

  tf_use <- tf_mat[, valid]
  score_use <- score[valid]

  cat("Cells used:", length(score_use), "\n")

  
  if (length(score_use) < 100) {
    cat("Skipping (too few cells)\n")
    next
  }

  # ------------------------------
  # Design matrix
  # ------------------------------
  design <- model.matrix(~ scale(score_use))

  # ------------------------------
  # limma model
  # ------------------------------
  fit <- lmFit(tf_use, design)
  fit <- eBayes(fit)

  res <- topTable(fit, coef = 2, number = Inf)

  if (nrow(res) == 0) {
    cat(" No results for:", sig, "\n")
    next
  }

  res$TF <- rownames(res)
  res$signature <- sig


  write.csv(
    res,
    paste0("chromVAR_lm_", sig, ".csv"),
    row.names = FALSE
  )

  all_results[[sig]] <- res
  
}


cat("\n[5] Merging results...\n")

all_df <- bind_rows(all_results)

write.csv(all_df, "chromVAR_lm_all_signatures.csv", row.names = FALSE)




