################################################## permutation test for crispr, chromatin region and cell state###########

##State2,5 upregulated chromatin wilcox test
library(Seurat)
library(Signac)
library(dplyr)
library(lme4)
library(Matrix)
library(GenomicRanges)
library(rtracklayer)


library(future)

plan("multicore", workers = 8) 
options(future.globals.maxSize = 50 * 1024^3)



setwd("~/Documents/Milan/Epigenomics/data/scatac")
sce_file <- "Cancer_scATACseq_data/GBMx_all.rds"
cicero_file <- "Cancer_scATACseq_data/GBMx_CICERO_GLOBAL.rds"

bulk_senescence_peaks_file <- "Cancer_scATACseq_data/atac_up.bed"
anchor1_file <- "anchor1_scatac.bed"
anchor2_file <- "anchor2_scatac.bed"


cat("\n[1] Loading data...\n")

sce <- readRDS(sce_file)
cicero <- readRDS(cicero_file)


cat("\n[7] ATAC_up integration...\n")

bulk_peaks <- import(bulk_senescence_peaks_file)
sc_peaks <- StringToGRanges(rownames(sce))

hits <- findOverlaps(sc_peaks, bulk_peaks)
length(hits)

idx <- unique(queryHits(hits))

# Extract peak names
peak_names <- rownames(sce[["peaks"]])

# Safety check
if (length(idx) == 0) stop("No overlapping peaks found!")
if (max(idx) > length(peak_names)) stop("Index exceeds peak names!")

# Add module score
sce <- AddModuleScore(
  object = sce,
  features = list(ATAC_up = peak_names[idx]),
  assay = "peaks",
  name = "ATAC_up"
)


# Rename for clarity
sce$ATAC_up_score <- sce$ATAC_up1

sce$ATAC_up1 <- NULL


score_columns <- c(
  "p16_top150_GA_score",
  "senmayo_GA_score",
  "state5_2fc_GA_score",
  "state2_2fc_GA_score"
)


all_da <- list()


for (col in score_columns) {

  cat("\n==============================\n")
  cat("Processing:", col, "\n")
  cat("==============================\n")


  q30 <- quantile(sce[[col]][,1], 0.30, na.rm = TRUE)
  q70 <- quantile(sce[[col]][,1], 0.70, na.rm = TRUE)

  class_col <- paste0("classif_", col)

  sce[[class_col]] <- ifelse(
    sce[[col]][,1] > q70, "High",
    ifelse(sce[[col]][,1] < q30, "Low", "Mid")
  )

  cat("High cells:", sum(sce[[class_col]] == "High"), "\n")
  cat("Low cells:", sum(sce[[class_col]] == "Low"), "\n")

  # ------------------------------
  # Subset only High/Low
  # ------------------------------

cells_use <- colnames(sce)[sce@meta.data[[class_col]] %in% c("High", "Low")]

sub <- subset(sce, cells = cells_use)


Idents(sub) <- factor(sub@meta.data[[class_col]])

# ------------------------------
# Differential accessibility
# ------------------------------
DefaultAssay(sub) <- "peaks"

da <- FindMarkers(
  sub,
  ident.1 = "High",
  ident.2 = "Low",
  test.use = "wilcox"
)

  da$peak <- rownames(da)
  da$signature <- col


  out_file <- paste0("DA_", col, ".csv")
  write.csv(da, out_file, row.names = FALSE)

  all_da[[col]] <- da
}


da_all <- bind_rows(all_da)

write.csv(da_all, "DA_all_signatures.csv", row.names = FALSE)

cat("\nDONE: All DA results saved\n")

#######################################################################################
##Enrichments
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(rtracklayer)
library(Matrix)
library(dplyr)
library(lme4)
library(future)
library(data.table)

plan("multicore", workers = 8)
options(future.globals.maxSize = 50 * 1024^3)

set.seed(1)

setwd("~/Documents/Milan/Epigenomics/data/scatac")



obj <- readRDS("Cancer_scATACseq_data/GBMx_all.rds")

cicero <- readRDS("Cancer_scATACseq_data/GBMx_CICERO_GLOBAL.rds")

da_state2 <- read.csv("DA_state2_2fc_GA_score.csv")
da_state5 <- read.csv("DA_state5_2fc_GA_score.csv")

state2_sig <- da_state2 %>%
  dplyr::filter(p_val_adj < 0.05)
  
state5_sig <- da_state5 %>%
  dplyr::filter(p_val_adj < 0.05)



bulk_gr  <- import("Cancer_scATACseq_data/atac_up.bed")
ctcf_gr  <- import("all_CTCF_raw.bed")

state2_gr <- StringToGRanges(state2_sig$peak, sep = c("-", "-"))
state5_gr <- StringToGRanges(state5_sig$peak, sep = c("-", "-"))



crispr_df <- data.frame(
  chr = c("chr3", "chr15"),
  start = c(46089058, 46090453),
  end = c(73955958, 73957981),
  gene = c("CCR1", "CD276")
)



crispr_gr <- GRanges(
  seqnames = crispr_df$chr,
  ranges   = IRanges(crispr_df$start, crispr_df$end),
  gene     = crispr_df$gene
)



peak_names <- rownames(obj[["peaks"]])
peak_gr <- StringToGRanges(peak_names, sep = c("-", "-"))

hits_bulk <- findOverlaps(peak_gr, bulk_gr)

bulk_idx <- unique(queryHits(hits_bulk))

if (length(bulk_idx) == 0) stop("No overlap between scATAC peaks and bulk ATAC_up")

obj <- AddModuleScore(
  obj,
  features = list(bulk_ATAC = peak_names[bulk_idx]),
  assay = "peaks",
  name = "ATAC_up"
)

obj$ATAC_up_score <- obj$ATAC_up1
obj$ATAC_up1 <- NULL

meta <- obj@meta.data

############################################################
#CTCF ENRICHMENT
############################################################


# ==============================
# BACKGROUND PEAK SET
# ==============================
all_peaks <- rownames(obj[["peaks"]])

all_gr <- StringToGRanges(all_peaks, sep = c("-", "-"))

# ==============================
# OBSERVED OVERLAP
# ==============================
hits <- findOverlaps(crispr_gr, ctcf_gr)
obs <- length(unique(queryHits(hits)))

cat("\nObserved overlaps:", obs, "\n")

# ==============================
# PERMUTATION TEST
# ==============================
set.seed(1)
n_perm <- 10000

perm_dist <- replicate(n_perm, {

  idx <- sample(seq_along(all_gr), length(crispr_gr))
  shuf <- all_gr[idx]

  length(unique(queryHits(
    findOverlaps(shuf, ctcf_gr)
  )))
})


emp_p <- mean(perm_dist >= obs)
fold_enrichment <- obs / mean(perm_dist)


cat("\n========== RESULTS ==========\n")
cat("Observed overlap:", obs, "\n")
cat("Mean permuted overlap:", mean(perm_dist), "\n")
cat("Fold enrichment:", fold_enrichment, "\n")
cat("Empirical P-value:", emp_p, "\n")


results <- data.frame(
  observed = obs,
  mean_perm = mean(perm_dist),
  fold_enrichment = fold_enrichment,
  empirical_p = emp_p
)

write.csv(results, "CRISPR_CTCF_permutation_results.csv", row.names = FALSE)

write.csv(
  data.frame(null = perm_dist),
  "CRISPR_CTCF_null_distribution.csv",
  row.names = FALSE
)

cat("\nDONE\n")



############################################################
#STATE5 ENRICHMENT vs CRISPR LOCI
############################################################

library(GenomicRanges)
library(Signac)

# ------------------------------
# BACKGROUND PEAK SPACE
# ------------------------------
all_peaks <- rownames(obj[["peaks"]])

all_gr <- StringToGRanges(
  all_peaks,
  sep = c("-", "-")
)



# ------------------------------
# OBSERVED OVERLAP
# ------------------------------
hits <- findOverlaps(crispr_gr, state5_gr)
obs <- length(unique(queryHits(hits)))

cat("\nObserved state5 overlaps:", obs, "\n")

# ------------------------------
# PERMUTATION TEST (10,000)
# ------------------------------
set.seed(1)
n_perm <- 10000

perm_state5 <- replicate(n_perm, {

  shuf <- sample(all_gr, length(crispr_gr))

  length(unique(queryHits(
    findOverlaps(shuf, state5_gr)
  )))
})


mean_perm <- mean(perm_state5)
emp_p <- mean(perm_state5 >= obs)
fold_enrichment <- obs / mean_perm


cat("\n========== STATE5 ENRICHMENT ==========\n")
cat("Observed overlap:", obs, "\n")
cat("Mean permuted overlap:", mean_perm, "\n")
cat("Fold enrichment:", fold_enrichment, "\n")
cat("Empirical P-value:", emp_p, "\n")


write.csv(
  data.frame(
    observed = obs,
    mean_perm = mean_perm,
    fold_enrichment = fold_enrichment,
    empirical_p = emp_p
  ),
  "CRISPR_state5_permutation_results.csv",
  row.names = FALSE
)

write.csv(
  data.frame(null = perm_state5),
  "CRISPR_state5_null_distribution.csv",
  row.names = FALSE
)

cat("\nDONE\n")
###################################################state2 permutation Loop
library(GenomicRanges)
library(Signac)

# ------------------------------
# BACKGROUND PEAK SPACE
# ------------------------------
all_peaks <- rownames(obj[["peaks"]])

all_gr <- StringToGRanges(
  all_peaks,
  sep = c("-", "-")
)



# ------------------------------
# OBSERVED OVERLAP
# ------------------------------
hits <- findOverlaps(crispr_gr, state2_gr)
obs <- length(unique(queryHits(hits)))

cat("\nObserved state2 overlaps:", obs, "\n")

# ------------------------------
# PERMUTATION TEST (10,000)
# ------------------------------
set.seed(1)
n_perm <- 10000

perm_state2 <- replicate(n_perm, {

  # random CRISPR-sized sampling from peak universe
  shuf <- sample(all_gr, length(crispr_gr))

  length(unique(queryHits(
    findOverlaps(shuf, state2_gr)
  )))
})


mean_perm <- mean(perm_state2)
emp_p <- mean(perm_state2 >= obs)
fold_enrichment <- obs / mean_perm


cat("\n========== STATE2 ENRICHMENT ==========\n")
cat("Observed overlap:", obs, "\n")
cat("Mean permuted overlap:", mean_perm, "\n")
cat("Fold enrichment:", fold_enrichment, "\n")
cat("Empirical P-value:", emp_p, "\n")


write.csv(
  data.frame(
    observed = obs,
    mean_perm = mean_perm,
    fold_enrichment = fold_enrichment,
    empirical_p = emp_p
  ),
  "CRISPR_state2_permutation_results.csv",
  row.names = FALSE
)

write.csv(
  data.frame(null = perm_state2),
  "CRISPR_state2_null_distribution.csv",
  row.names = FALSE
)

cat("\nDONE\n")

