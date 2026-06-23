#!/usr/bin/env Rscript
library(Seurat)
library(nebula)
library(RColorBrewer)
library(gplots)
library("ggplot2")
library("gridExtra")
library(ggplot2)
library(dplyr)
library(readr)

load("/home/user/Documents/Milan/Epigenomics/data/scRNA_seq_GBM4set_mice2set/GSE182109/gbm_GSE182109_filtered_scores.RData")

# Define your gene score columns
score_columns <- c(
  "p16_top150.score1",
  "senmayo.score1",
  "IR_top150.score1",
)


dir <- "~/nebula_gbm/"
if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)


for (col in score_columns) {
  message(">>> Running NEBULA for: ", col)
  t1 <- Sys.time()
  
  
  classif_col <- paste0("classif_", col)
  q25 <- quantile(se@meta.data[[col]], probs = 0.3, na.rm = TRUE)
  q75 <- quantile(se@meta.data[[col]], probs = 0.7, na.rm = TRUE)
  se@meta.data[[classif_col]] <- ifelse(
    se@meta.data[[col]] > q75, "high",
    ifelse(se@meta.data[[col]] < q25, "low", NA)
  )
  
  
  subset_cells <- !is.na(se@meta.data[[classif_col]])
  se_sub <- subset(se, cells = colnames(se)[subset_cells])
  se_sub@meta.data[[classif_col]] <- factor(se_sub@meta.data[[classif_col]], levels = c("low", "high"))
  
 
  data_neb <- scToNeb(obj = se_sub, assay = "RNA", id = "Patient", pred = classif_col)
  
 
  pred_df <- data_neb$pred
  colnames(pred_df) <- "classif"
  pred_matrix <- model.matrix(~ classif, data = pred_df)
  
  
  count_mat <- data_neb$count
  id_vec <- as.factor(data_neb$id)
  offset_vec <- log1p(colSums(count_mat))
  names(offset_vec) <- NULL

  # Run NEBULA safely
  neb_res <- try(
    nebula(
      count = count_mat,
      id = id_vec,
      pred = pred_matrix,
      offset = offset_vec,
      verbose = TRUE
    ),
    silent = TRUE
  )
  
  #Check NEBULA result
  if (inherits(neb_res, "try-error")) {
    message("NEBULA failed for ", col, ": ", conditionMessage(attr(neb_res, "condition")))
    next
  }

  if (!is.null(neb_res$summary) && nrow(neb_res$summary) > 0) {
    neb_summary <- neb_res$summary
    
    # Save full result
    full_path <- file.path(dir, paste0("NEBULA_gbm_", col, ".txt"))
    write.table(neb_summary, file = full_path, sep = "\t", quote = FALSE, row.names = FALSE)

    # Filter significant genes
    sig_genes <- subset(neb_summary, p_classifhigh < 0.001 & abs(logFC_classifhigh) > 0.5)

    if (nrow(sig_genes) > 0) {
      sorted_sig_genes <- sig_genes[order(-sig_genes$logFC_classifhigh), ]
      
     
      filter_path <- file.path(dir, paste0("NEBULA_gbm_filter_", col, ".txt"))
      write.table(sorted_sig_genes, file = filter_path, sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
      message("NEBULA finished but no significant genes for ", col)
    }

    message("Finished ", col, " in ", round(Sys.time() - t1, 2), " seconds.")

  } else {
    message("NEBULA returned no usable results for ", col)
  }
}


#################testing padj values####

library(dplyr)

# -----------------------------------------
# directory containing NEBULA outputs
# -----------------------------------------
dir <- "C:/Users/mlnpa/Documents/nebula_gbm/"


# -----------------------------------------
# score names
# -----------------------------------------
score_columns <- c(
  "p16_top150.score1",
  "senmayo.score1",
  "IR_top150.score1"
)

# -----------------------------------------
# loop through analyses
# -----------------------------------------
for (col in score_columns) {

  message("\n==============================")
  message("Checking: ", col)

  # -----------------------------------------
  # load full result
  # -----------------------------------------
  full_file <- file.path(
    dir,
    paste0("NEBULA_gbm_", col, ".txt")
  )

  neb_summary <- read.table(
    full_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  # -----------------------------------------
  # calculate adjusted p-values
  # -----------------------------------------
  neb_summary$padj_classifhigh <- p.adjust(
    neb_summary$p_classifhigh,
    method = "BH"
  )

  # -----------------------------------------
  # original filtering
  # -----------------------------------------
  old_filter <- subset(
    neb_summary,
    p_classifhigh < 0.001 &
    abs(logFC_classifhigh) > 0.5
  )


  old_filter$pass_fdr <- old_filter$padj_classifhigh < 0.05

  # -----------------------------------------
  # summary statistics
  # -----------------------------------------
  n_total <- nrow(old_filter)
  n_fdr <- sum(old_filter$pass_fdr)

  message("Genes passing original filter: ", n_total)
  message("Genes also passing padj < 0.05: ", n_fdr)

  message(
    "Percent surviving FDR: ",
    round(100 * n_fdr / n_total, 2),
    "%"
  )


  out_file <- file.path(
    dir,
    paste0("NEBULA_gbm_filter_FDRcheck_", col, ".txt")
  )

  write.table(
    old_filter,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

}