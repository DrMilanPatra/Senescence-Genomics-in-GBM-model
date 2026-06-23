#!/usr/bin/env Rscript

library(Seurat)
library(nebula)

# Load the precomputed Seurat object with scores
load("/home/user/Documents/Milan/Epigenomics/data/scRNA_seq_GBM4set_mice2set/GSE182109/gbm_GSE182109_filtered_scores.RData")


# Defining gene score columns
score_columns <- c(
  "p16_top150.score1",
  "senmayo.score1",
  "IR_top150.score1",

)



p53_mut_map <- c(
  "CNSTM-068" = "mut",
  "CNSTM-096" = "mut",
  "CNSTM-390" = "mut",
  "CNSTM-397" = "mut",
  "MDAG-9"    = "mut"
)


tert_mut_map <- c(
  "CNSTM-394" = "mut",
  "CNSTM-397" = "mut",
  "MDAG-4" = "mut"
)

pten_mut_map <- c(
  "CNSTM-068" = "mut",
  "CNSTM-070" = "mut",
  "CNSTM-394" = "mut",
  "MDAG-9"    = "mut"
)



mgmt_map <- c(
  "CNSTM-068" = "methylated",
  "CNSTM-397" = "methylated",
  "MDAG-4" = "methylated",
  "MDAG-12" = "methylated",
  "MDAG-1" = "methylated",
  "CNSTM-390" = "unmethylated",
  "CNSTM-394" = "unmethylated"
)



se@meta.data$p53_mut <- ifelse(
  se@meta.data$orig_ident %in% names(p53_mut_map),
  "mut", "wt"
)

se@meta.data$tert_mut <- ifelse(
  se@meta.data$orig_ident %in% names(tert_mut_map),
  "mut", "wt"
)

se@meta.data$pten_mut <- ifelse(
  se@meta.data$orig_ident %in% names(pten_mut_map),
  "mut", "wt"
)


se@meta.data$mgmt_status <- mgmt_map[as.character(se@meta.data$orig_ident)]
se@meta.data$mgmt_status[is.na(se@meta.data$mgmt_status)] <- "not_obtained"


dir <- "~/nebula_gbm1/"
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
se_sub@meta.data$p53_mut <- factor(se_sub@meta.data$p53_mut, levels = c("wt", "mut"))
se_sub@meta.data$tert_mut <- factor(se_sub@meta.data$tert_mut, levels = c("wt", "mut"))
se_sub@meta.data$pten_mut <- factor(se_sub@meta.data$pten_mut, levels = c("wt", "mut"))

se_sub@meta.data$mgmt_status <- factor(
  se_sub@meta.data$mgmt_status,
  levels = c("unmethylated", "methylated", "not_obtained")
)
  

  data_neb <- scToNeb(
    obj = se_sub,
    assay = "RNA",
    id = "Patient",
    pred = c(classif_col, "mgmt_status", "p53_mut", "pten_mut", "tert_mut")
  )

  
  pred_df <- data_neb$pred
  colnames(pred_df)[colnames(pred_df) == classif_col] <- "classif"  # for formula compatibility

  # Model: gene expression ~ score class + fixed effects
  pred_matrix <- model.matrix(~ classif + mgmt_status + p53_mut + pten_mut + tert_mut, data = pred_df)
  
 
  count_mat <- data_neb$count
  id_vec <- as.factor(data_neb$id)
  offset_vec <- log1p(colSums(count_mat))
  names(offset_vec) <- NULL


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

  
  if (inherits(neb_res, "try-error")) {
    message("NEBULA failed for ", col, ": ", conditionMessage(attr(neb_res, "condition")))
    next
  }

  if (!is.null(neb_res$summary) && nrow(neb_res$summary) > 0) {
    neb_summary <- neb_res$summary
    
    
    full_path <- file.path(dir, paste0("NEBULA_gbm_", col, ".txt"))
    write.table(neb_summary, file = full_path, sep = "\t", quote = FALSE, row.names = FALSE)

    
    if ("p_classifhigh" %in% colnames(neb_summary)) {
      sig_genes <- subset(neb_summary, p_classifhigh < 0.001 & abs(logFC_classifhigh) > 0.5)
      if (nrow(sig_genes) > 0) {
        sorted_sig_genes <- sig_genes[order(-sig_genes$logFC_classifhigh), ]
        filter_path <- file.path(dir, paste0("NEBULA_gbm_filter_", col, ".txt"))
        write.table(sorted_sig_genes, file = filter_path, sep = "\t", quote = FALSE, row.names = FALSE)
      } else {
        message("NEBULA finished but no significant genes for ", col)
      }
    } else {
      message("p_classifhigh not found in NEBULA result for ", col)
    }

    message("Finished ", col, " in ", round(Sys.time() - t1, 2), " seconds.")
  } else {
    message("NEBULA returned no usable results for ", col)
  }
}


#################testing padj values####

library(dplyr)

dir <- "C:/Users/mlnpa/Documents/nebula_gbm1/"



score_columns <- c(
  "p16_top150.score1",
  "senmayo.score1",
  "IR_top150.score1"
)

for (col in score_columns) {

  message("\n==============================")
  message("Checking: ", col)


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

  # -----------------------------------------
  # check FDR survival
  # -----------------------------------------
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

  # -----------------------------------------
  # save checked table
  # -----------------------------------------
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