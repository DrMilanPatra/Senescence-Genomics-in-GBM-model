

#!/usr/bin/env Rscript

rm(list = ls())


library(Signac)
library(Seurat)
library(Matrix)
library(monocle)
library(cicero)
library(GenomicRanges)
library(AnnotationHub)
library(future)



options(mc.cores = 8)

Sys.setenv(
  OMP_NUM_THREADS = 8,
  MKL_NUM_THREADS = 8,
  OPENBLAS_NUM_THREADS = 8
)

setwd("~/Documents/Milan/Epigenomics/data/scatac")

input_file <- "Cancer_scATACseq_data/GBMx_all.rds"
output_dir <- "processed_seurat_objects/cicero/gliomacell_ciciero"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\nCICERO PIPELINE (GBMx_all SPLIT BY SAMPLE)\n")


cat("\nLoading GBMx_all...\n")
obj <- readRDS(input_file)

DefaultAssay(obj) <- "peaks"

cat("Total cells:", ncol(obj), "\n")

# ==============================
# SPLIT BY SAMPLE
# ==============================
cat("\nSplitting by sample...\n")

obj_list <- SplitObject(obj, split.by = "sample")

cat("Total samples:", length(obj_list), "\n")


cat("\nLoading hg38 genome...\n")

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


# ==============================
# MAIN LOOP (PER SAMPLE)
# ==============================
for (sample_id in names(obj_list)) {

  out_file <- file.path(output_dir, paste0(sample_id, "_cicero.rds"))

  if (file.exists(out_file)) {
    cat("\nSkipping:", sample_id, "\n")
    next
  }

  cat("\n==============================\n")
  cat("Sample:", sample_id, "\n")
  cat("==============================\n")

  atac <- obj_list[[sample_id]]
  DefaultAssay(atac) <- "peaks"

  peak_mat <- GetAssayData(atac, assay = "peaks", layer = "counts")

  cat("Raw matrix:", dim(peak_mat), "\n")

  # ------------------------------
  # Q1 FILTER
  # ------------------------------
  peak_sums <- Matrix::rowSums(peak_mat)
  cutoff <- quantile(peak_sums, 0.25)

  peak_mat <- peak_mat[peak_sums > cutoff, , drop = FALSE]
  atac <- atac[rownames(peak_mat), ]

  cat("Peaks after filtering:", dim(peak_mat), "\n")

  if (nrow(peak_mat) < 5000) {
    cat("Too few peaks — skipping\n")
    next
  }


  peak_mat@x[peak_mat@x > 0] <- 1

  cat("Cells used:", ncol(peak_mat), "\n")


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

  cat("CDS created:", dim(exprs(cds)), "\n")

  # ------------------------------
  # LSI (PER SAMPLE)
  # ------------------------------
  cat("\nRunning LSI...\n")

  atac <- RunTFIDF(atac)
  atac <- FindTopFeatures(atac, min.cutoff = 'q0')
  atac <- RunSVD(atac)

  reduced_coords <- Embeddings(atac, "lsi")[colnames(peak_mat), 1:30]


  cicero_cds <- make_cicero_cds(
    cds,
    reduced_coordinates = reduced_coords
  )

  cat("Cicero CDS ready\n")


  cat("\nRunning Cicero...\n")

  conns <- run_cicero(
    cicero_cds,
    genomic_coords = genome.df,
    sample_num = 200
  )

  cat("Cicero completed\n")

  # ------------------------------
  # DEBUG OUTPUT
  # ------------------------------
  cat("\n Connections preview:\n")
  print(head(conns))

  cat("\nTotal connections:", nrow(conns), "\n")
  cat("\nCo-accessibility summary:\n")
  print(summary(conns$coaccess))


  saveRDS(conns, out_file)

  cat("Saved:", sample_id, "\n")

  rm(atac, cds, cicero_cds, conns, peak_mat)
  gc()
}

cat("\n CICERO COMPLETED SUCCESSFULLY\n")

######################################################################
#!/usr/bin/env Rscript

rm(list = ls())

library(Signac)
library(Seurat)
library(Matrix)
library(monocle)
library(cicero)
library(GenomicRanges)
library(AnnotationHub)
library(future)



options(mc.cores = 8)

Sys.setenv(
  OMP_NUM_THREADS = 8,
  MKL_NUM_THREADS = 8,
  OPENBLAS_NUM_THREADS = 8
)


setwd("~/Documents/Milan/Epigenomics/data/scatac")

cicero_dir <- "processed_seurat_objects/cicero/gliomacell_ciciero"
out_file   <- "Cancer_scATACseq_data/GBMx_CICERO_GLOBAL.rds"


files <- list.files(
  cicero_dir,
  pattern = "_cicero\\.rds$",
  full.names = TRUE
)

cat("Total Cicero files:", length(files), "\n")


all_conns <- list()


for (f in files) {

  cat("\nProcessing:", basename(f), "\n")

  conns <- readRDS(f)

  if (!all(c("Peak1", "Peak2", "coaccess") %in% colnames(conns))) {
    cat("Invalid format. Skipping\n")
    next
  }


  conns$Peak1 <- as.character(conns$Peak1)
  conns$Peak2 <- as.character(conns$Peak2)

  # ------------------------------
  # ROBUST PAIR HARMONIZATION
  # ------------------------------
  conns$pair_id <- apply(
    conns[, c("Peak1", "Peak2")],
    1,
    function(x) paste(sort(x), collapse = "_")
  )


  sample_id <- gsub("_cicero\\.rds$", "", basename(f))
  conns$sample <- sample_id




conns <- conns %>% dplyr::filter(!is.na(coaccess) & coaccess > 0) %>% dplyr::select(pair_id, coaccess, sample)


  all_conns[[sample_id]] <- conns[, c("pair_id", "coaccess", "sample")]

  rm(conns)
  gc()
}


cat("\n[MERGE]\n")

all_conns <- bind_rows(all_conns)

cat("Total rows:", nrow(all_conns), "\n")

# ==============================
# CONSENSUS NETWORK
# ==============================
cat("\n[CONSENSUS NETWORK]\n")

cicero_global <- all_conns %>%
  group_by(pair_id) %>%
  summarise(
    mean_coaccess   = mean(coaccess, na.rm = TRUE),
    median_coaccess = median(coaccess, na.rm = TRUE),
    sd_coaccess     = sd(coaccess, na.rm = TRUE),
    n_samples       = n_distinct(sample),
    .groups = "drop"
  )

# ==============================
# REPRODUCIBILITY
# ==============================
total_samples <- length(files)

cicero_global <- cicero_global %>%
  mutate(
    reproducibility = n_samples / total_samples
  )

cat("Unique pairs:", nrow(cicero_global), "\n")

# ==============================
# SPLIT BACK PEAKS
# ==============================
split_pairs <- strsplit(cicero_global$pair_id, "_")

cicero_global$Peak1 <- sapply(split_pairs, `[`, 1)
cicero_global$Peak2 <- sapply(split_pairs, `[`, 2)

# ==============================
# DISTANCE CALCULATION (≤1 Mb)
# ==============================
get_dist <- function(p1, p2) {
  s1 <- as.numeric(strsplit(p1, "-")[[1]][2])
  s2 <- as.numeric(strsplit(p2, "-")[[1]][2])
  abs(s1 - s2)
}

cicero_global$distance <- mapply(
  get_dist,
  cicero_global$Peak1,
  cicero_global$Peak2
)

saveRDS(cicero_global, out_file)

  
cicero_filtered <- cicero_global %>% dplyr::filter(reproducibility >= 0.5 & distance <= 1e6 & median_coaccess > 0.05)

saveRDS(cicero_filtered, "Cancer_scATACseq_data/GBMx_CICERO_FILTERED.rds")



cat("\nDONE: Global Cicero network done!)\n")
###################cicero part result analysis for integration of cicero DE to SE-MP###########

######################################################################

rm(list = ls())

library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(GenomicRanges)
library(IRanges)
library(Matrix)
library(limma)


setwd("~/Documents/Milan/Epigenomics/data/scatac")

cicero_file <- "Cancer_scATACseq_data/GBMx_CICERO_GLOBAL.rds"
obj_file    <- "Cancer_scATACseq_data/GBMx_all.rds"
mp_rds      <- "/home/user/Documents/Milan/Epigenomics/MP_downstream.rds"
hic_file    <- "GBM_HiChIP_loops_mean1.txt"

crispr_gr <- GRanges("chr5", IRanges(132202596, 132228449))


n_perm <- 2000
set.seed(1)

# ==============================
# LOAD CICERO EDGES
# ==============================
cat("\nLoading Cicero...\n")

cicero_global <- readRDS(cicero_file) %>%
  filter(
    reproducibility >= 0.5,
    distance <= 1e6,
    median_coaccess > 0.01
  )

cat("Edges after filtering:", nrow(cicero_global), "\n")

# ------------------------------
# PRESERVE ORIGINAL pair_id
# ------------------------------
cicero_global$pair_id_original <- cicero_global$pair_id

# ------------------------------
# DERIVE CANONICAL pair_id (FROM PEAKS)
# ------------------------------
cicero_global$pair_id_canonical <- paste(
  pmin(cicero_global$Peak1, cicero_global$Peak2),
  pmax(cicero_global$Peak1, cicero_global$Peak2),
  sep = "_"
)

# ==============================
# QC: VALIDATE pair_id INTEGRITY
# ==============================
cat("\n[QC] Validating pair_id integrity (GLOBAL)...\n")

mismatch_idx <- which(
  cicero_global$pair_id_original != cicero_global$pair_id_canonical
)

cat("Total edges:", nrow(cicero_global), "\n")
cat("Mismatches:", length(mismatch_idx), "\n")

if (length(mismatch_idx) > 0) {
  cat("ERROR: pair_id mismatch detected!\n")
  print(head(cicero_global[mismatch_idx, ], 5))
  stop("pair_id validation failed (GLOBAL)")
} else {
  cat("pair_id integrity verified (GLOBAL)\n")
}


cicero_global$pair_id <- cicero_global$pair_id_canonical

# ==============================
# CREATE EDGE TABLE (POST-QC ONLY)
# ==============================
edge_table <- cicero_global %>%
  select(pair_id, Peak1, Peak2, median_coaccess)

  
  
# ==============================
#  SE–MP MODEL
# ==============================
cat("\n SE–MP mapping...\n")

mp_df <- readRDS(mp_rds)$SE_MP_coords

mp_gr <- GRanges(
  seqnames = mp_df$chr_se,
  ranges   = IRanges(mp_df$start_se, mp_df$end_se),
  MP       = mp_df$meta_program
)

# ==============================
# ANCHOR MAPPING
# ==============================
cat("\n Anchor mapping...\n")

gr1 <- StringToGRanges(edge_table$Peak1)
gr2 <- StringToGRanges(edge_table$Peak2)

mcols(gr1)$pair_id <- edge_table$pair_id
mcols(gr2)$pair_id <- edge_table$pair_id

ov1 <- findOverlaps(gr1, mp_gr)
ov2 <- findOverlaps(gr2, mp_gr)

anchor_map <- bind_rows(
  data.frame(pair_id=mcols(gr1)$pair_id[queryHits(ov1)],
             anchor="P1",
             MP=mp_gr$MP[subjectHits(ov1)]),
  data.frame(pair_id=mcols(gr2)$pair_id[queryHits(ov2)],
             anchor="P2",
             MP=mp_gr$MP[subjectHits(ov2)])
)

# ==============================
# PROBABILISTIC SE–MP EDGE MODEL
# ==============================
cat("\n Building SE–MP edge model...\n")

# anchor-level probabilities
anchor_map <- anchor_map %>%
  group_by(pair_id, anchor, MP) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(pair_id, anchor) %>%
  mutate(p_anchor = count / sum(count)) %>%
  ungroup()


edge_mp <- anchor_map %>%
  group_by(pair_id, MP) %>%
  summarise(weight = mean(p_anchor), .groups = "drop") %>%
  group_by(pair_id) %>%
  mutate(weight = weight / sum(weight)) %>%
  ungroup()

edge_summary <- edge_mp %>%
  group_by(pair_id) %>%
  summarise(
    MP_entropy = -sum(weight * log(weight + 1e-10)),
    MP_max = max(weight),
    MP_diversity = n_distinct(MP),
    .groups="drop"
  )
  
  
  edge_dominant_mp <- edge_mp %>%
  group_by(pair_id) %>%
  slice_max(weight, n = 1, with_ties = FALSE) %>%
  select(pair_id, dominant_MP = MP)
  
  
edge_table$SE_flag <- edge_table$pair_id %in% edge_mp$pair_id

cat("SE-linked edges:", sum(edge_table$SE_flag), "\n")
cat("Non-SE edges:", sum(!edge_table$SE_flag), "\n")

edge_table <- left_join(edge_table, edge_summary, by="pair_id")
edge_table <- left_join(edge_table, edge_dominant_mp, by = "pair_id")


edge_table$dominant_MP[is.na(edge_table$dominant_MP)] <- "None"
edge_table$MP_entropy[is.na(edge_table$MP_entropy)] <- 0
edge_table$MP_max[is.na(edge_table$MP_max)] <- 0
edge_table$MP_diversity[is.na(edge_table$MP_diversity)] <- 0

# ==============================
# HiChIP AS PRIOR (LOOP-AWARE, CLEAN)
# ==============================
cat("\n HiChIP prior (loop-aware)...\n")

hic <- fread(hic_file)

# assign loop IDs (stable)
hic$loop_id <- if ("V8" %in% colnames(hic)) hic$V8 else paste0("loop_", seq_len(nrow(hic)))

# build GRanges with loop identity
hic_gr1 <- GRanges(
  seqnames = hic$V1,
  ranges   = IRanges(hic$V2, hic$V3),
  loop_id  = hic$loop_id
)

hic_gr2 <- GRanges(
  seqnames = hic$V4,
  ranges   = IRanges(hic$V5, hic$V6),
  loop_id  = hic$loop_id
)

# Cicero anchors (DO THIS ONCE)
gr1_full <- StringToGRanges(edge_table$Peak1)
gr2_full <- StringToGRanges(edge_table$Peak2)

mcols(gr1_full)$edge_idx <- seq_len(nrow(edge_table))
mcols(gr2_full)$edge_idx <- seq_len(nrow(edge_table))

# overlaps
ov1 <- findOverlaps(gr1_full, hic_gr1)
ov2 <- findOverlaps(gr2_full, hic_gr2)

# Peak1 → loop_id
map1 <- data.frame(
  edge_idx = mcols(gr1_full)$edge_idx[queryHits(ov1)],
  loop_id  = mcols(hic_gr1)$loop_id[subjectHits(ov1)]
)

# Peak2 → loop_id
map2 <- data.frame(
  edge_idx = mcols(gr2_full)$edge_idx[queryHits(ov2)],
  loop_id  = mcols(hic_gr2)$loop_id[subjectHits(ov2)]
)

# forward orientation
loop_support_fwd <- merge(map1, map2, by = c("edge_idx", "loop_id"))

# =========================
# ADD SWAP ORIENTATION HERE
# =========================

ov1_swap <- findOverlaps(gr1_full, hic_gr2)
ov2_swap <- findOverlaps(gr2_full, hic_gr1)

map1_swap <- data.frame(
  edge_idx = mcols(gr1_full)$edge_idx[queryHits(ov1_swap)],
  loop_id  = mcols(hic_gr2)$loop_id[subjectHits(ov1_swap)]
)

map2_swap <- data.frame(
  edge_idx = mcols(gr2_full)$edge_idx[queryHits(ov2_swap)],
  loop_id  = mcols(hic_gr1)$loop_id[subjectHits(ov2_swap)]
)

loop_support_swap <- merge(map1_swap, map2_swap, by = c("edge_idx", "loop_id"))

# =========================
# COMBINE BOTH DIRECTIONS
# =========================

loop_support <- rbind(loop_support_fwd, loop_support_swap)

# now define supported edges
supported_edges <- unique(loop_support$edge_idx)

# assign HiChIP support
edge_table$HiChIP <- 0L
edge_table$HiChIP[supported_edges] <- 1L

cat("HiChIP-supported edges:", sum(edge_table$HiChIP), "\n")


edge_table$regulatory_prior <- edge_table$MP_max * (1 + edge_table$HiChIP)

# ==============================
# CRISPR LOCUS OVERLAY (CONSISTENT)
# ==============================
cat("\nCRISPR locus mapping...\n")



# attach pair_id once
mcols(gr1_full)$pair_id <- edge_table$pair_id
mcols(gr2_full)$pair_id <- edge_table$pair_id

# overlaps
ov_c1 <- findOverlaps(gr1_full, crispr_gr)
ov_c2 <- findOverlaps(gr2_full, crispr_gr)

# collect edges touching CRISPR locus
crispr_edges <- unique(c(
  mcols(gr1_full)$pair_id[queryHits(ov_c1)],
  mcols(gr2_full)$pair_id[queryHits(ov_c2)]
))


edge_table$crispr_edge <- edge_table$pair_id %in% crispr_edges

cat("CRISPR-overlapping edges:", sum(edge_table$crispr_edge), "\n")




saveRDS(edge_table, "edge_table_annotated.rds")


saveRDS(edge_mp, "edge_mp_weights.rds")
saveRDS(mp_gr, "mp_granges.rds")

cat("Saved:\n")
cat("- edge_table_annotated.rds\n")
cat("- edge_mp_weights.rds\n")
cat("- mp_granges.rds\n")



################ limm model for CICERO CO-ACCESSIBILITY and cell state signature################
#!/usr/bin/env Rscript

rm(list = ls())
gc()


suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(Matrix)
  library(matrixStats)
  library(limma)
  library(Seurat)
})

set.seed(1)


setwd("~/Documents/Milan/Epigenomics/data/scatac")

cicero_file <- "Cancer_scATACseq_data/GBMx_CICERO_GLOBAL.rds"
obj_file    <- "Cancer_scATACseq_data/GBMx_all_with_chromvar.rds"

cicero_dir  <- "processed_seurat_objects/cicero/gliomacell_ciciero"
out_dir     <- "processed_seurat_objects/cicero/diff_coaccess_no_distance_model"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


cat("\nLoading data...\n")

cicero_global <- readRDS(cicero_file)
obj <- readRDS(obj_file)

cat("Total global edges:", nrow(cicero_global), "\n")


cicero_global <- cicero_global %>%
  filter(
    reproducibility >= 0.5,
    distance <= 1e6,
    median_coaccess > 0.01
  )

cat("Filtered edges:", nrow(cicero_global), "\n")


cat("\n[QC] Validating pair_id integrity (GLOBAL)...\n")

pair_id_check <- paste(
  pmin(cicero_global$Peak1, cicero_global$Peak2),
  pmax(cicero_global$Peak1, cicero_global$Peak2),
  sep = "_"
)

mismatch_idx <- which(cicero_global$pair_id != pair_id_check)

cat("Total edges:", nrow(cicero_global), "\n")
cat("Mismatches:", length(mismatch_idx), "\n")

if (length(mismatch_idx) > 0) {
  print(head(cicero_global[mismatch_idx, ], 5))
  stop("pair_id validation failed (GLOBAL)")
} else {
  cat("pair_id integrity verified (GLOBAL)\n")
}


keep_pairs <- unique(cicero_global$pair_id)
keep_hash  <- setNames(rep(TRUE, length(keep_pairs)), keep_pairs)

# ==============================
# BUILD SAMPLE-WISE EDGE TABLE
# ==============================
cat("\nBuilding per-sample edge table...\n")

files <- list.files(
  cicero_dir,
  pattern = "_cicero\\.rds$",
  full.names = TRUE
)

cat("Total samples:", length(files), "\n")

all_list <- vector("list", length(files))
names(all_list) <- basename(files)

for (i in seq_along(files)) {

  f <- files[i]
  sample_id <- sub("_cicero\\.rds$", "", basename(f))

  conns <- readRDS(f)

  if (!all(c("Peak1","Peak2","coaccess") %in% colnames(conns))) {
    cat("Skipping malformed:", sample_id, "\n")
    next
  }


  p1 <- as.character(conns$Peak1)
  p2 <- as.character(conns$Peak2)


  pair_id <- paste(
    pmin(p1, p2),
    pmax(p1, p2),
    sep = "_"
  )


  if ("pair_id" %in% colnames(conns)) {

    cat("[QC] Checking pair_id in:", sample_id, "\n")

    mismatch_idx <- which(conns$pair_id != pair_id)

    cat("Mismatches:", length(mismatch_idx), "\n")

    if (length(mismatch_idx) > 0) {
      print(head(conns[mismatch_idx, ], 5))
      stop(paste("pair_id mismatch detected in", sample_id))
    } else {
      cat("pair_id integrity OK\n")
    }
  }


  conns <- data.table(
    pair_id = pair_id,
    coaccess = as.numeric(conns$coaccess),
    sample = sample_id
  )


  conns <- conns[keep_hash[pair_id] == TRUE]

  if (nrow(conns) == 0) next

  all_list[[i]] <- conns
}

# ==============================
# MERGE + AGGREGATE
# ==============================
cat("\n Aggregating...\n")

all_dt <- rbindlist(all_list, fill = TRUE, use.names = TRUE)

cat("Total rows:", nrow(all_dt), "\n")

# collapse duplicates safely
all_dt <- all_dt[
  , .(coaccess = mean(coaccess, na.rm = TRUE)),
  by = .(pair_id, sample)
]

cat("Unique edges:", uniqueN(all_dt$pair_id), "\n")

# ==============================
# MATRIX CONSTRUCTION
# ==============================
cat("\nMatrix build...\n")

mat_df <- dcast(
  all_dt,
  pair_id ~ sample,
  value.var = "coaccess",
  fill = 0
)

pair_ids <- mat_df$pair_id
mat <- as.matrix(mat_df[, -1, with = FALSE])
rownames(mat) <- pair_ids

rm(mat_df, all_dt)
gc()

cat("Matrix dim:", dim(mat), "\n")

# ==============================
# FILTER LOW VARIANCE
# ==============================


pair_sd <- matrixStats::rowSds(mat)

keep_idx <- pair_sd > 0.001

mat <- mat[keep_idx, , drop = FALSE]
pair_ids <- pair_ids[keep_idx]

cat("Edges after filtering:", nrow(mat), "\n")

# ==============================
# SIGNATURE SCORING
# ==============================
cat("\n[6] Signature scoring...\n")

meta <- obj@meta.data

signatures <- c(
  "state5_2fc_GA_score",
  "state1_2fc_GA_score",
  "state2_2fc_GA_score"
)

signatures <- signatures[signatures %in% colnames(meta)]

sample_scores <- meta %>%
  group_by(sample) %>%
  summarise(across(all_of(signatures), ~ mean(.x, na.rm = TRUE)))

# ==============================
#ALIGN SAMPLES
# ==============================
cat("\nAligning samples...\n")

common <- intersect(colnames(mat), sample_scores$sample)

if (length(common) < 3) {
  stop("Too few samples after alignment — check data consistency")
}

mat <- mat[, common, drop = FALSE]
sample_scores <- sample_scores[match(common, sample_scores$sample), ]

stopifnot(all(colnames(mat) == sample_scores$sample))

cat("Samples used:", length(common), "\n")

# ==============================
# LIMMA MODEL 
# ==============================
cat("\nRunning limma...\n")

all_results <- list()

for (sig in signatures) {

  cat("\nProcessing:", sig, "\n")

  signature_vec <- sample_scores[[sig]]

  # guard against constant vector
  if (sd(signature_vec, na.rm = TRUE) == 0) {
    cat("Skipping (no variance):", sig, "\n")
    next
  }

  design <- model.matrix(~ scale(signature_vec))

  fit <- lmFit(mat, design)
  fit <- eBayes(fit)

  res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

 
  res$pair_id <- pair_ids
  res$signature <- sig

  # sanity check
  stopifnot(nrow(res) == length(pair_ids))

  fwrite(
    res,
    file.path(out_dir, paste0("diff_coaccess_", sig, ".csv"))
  )

  all_results[[sig]] <- res
}

# ==============================
# MERGE RESULTS
# ==============================
cat("\nMerging results...\n")

final <- rbindlist(all_results, use.names = TRUE)

fwrite(
  final,
  file.path(out_dir, "diff_coaccess_all_signatures.csv")
)

saveRDS(
  final,
  file.path(out_dir, "diff_coaccess_all_signatures.rds")
)

# ==============================
# SUMMARY
# ==============================
cat("\n====================================\n")
cat(" limma model for CICERO COMPLETE\n")
cat("Total tested edges:", nrow(final), "\n")
cat("Unique edges:", uniqueN(final$pair_id), "\n")
cat("====================================\n")
#####################################################################################3


#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(Matrix)
  library(limma)
})

set.seed(1)

setwd("~/Documents/Milan/Epigenomics/data/scatac")

edge_file <- "edge_table_annotated.rds"
diff_file <- "/home/user/Documents/Milan/Epigenomics/data/scatac/processed_seurat_objects/cicero/diff_coaccess_no_distance_model/diff_coaccess_all_signatures.rds"
out_dir   <- "processed_seurat_objects/cicero/diff_coaccess_no_distance_model"



# ==============================
# LOAD EDGE ANNOTATION
# ==============================


edge_table <- readRDS(edge_file)
edge_table$pair_id <- as.character(edge_table$pair_id)

cat("Edges in annotation:", nrow(edge_table), "\n")

# ==============================
#  LOAD DIFFERENTIAL RESULTS
# ==============================
cat("\n Loading differential Cicero...\n")

diff_res <- readRDS(diff_file)
diff_res$pair_id <- as.character(diff_res$pair_id)

# keep only required columns
diff_res <- diff_res %>%
  dplyr::select(pair_id, logFC, P.Value, adj.P.Val, signature)

cat("Total diff rows:", nrow(diff_res), "\n")
cat("Unique diff edges:", length(unique(diff_res$pair_id)), "\n")

# ==============================
# STRICT QC: EDGE CONSISTENCY
# ==============================
cat("\n QC: edge consistency check...\n")

common_edges <- intersect(edge_table$pair_id, diff_res$pair_id)

cat("Common edges:", length(common_edges), "\n")


if (length(common_edges) < 0.8 * length(unique(diff_res$pair_id))) {
  stop("Major mismatch: annotation vs differential space")
}

# ==============================
# ALIGN BOTH TABLES (STRICT ORDERING)
# ==============================
cat("\nAligning datasets...\n")

edge_annot <- edge_table %>%
  dplyr::select(pair_id, SE_flag, dominant_MP, MP_max,
                MP_entropy, MP_diversity, HiChIP,
                regulatory_prior, crispr_edge)

# filter both
diff_res   <- diff_res %>% filter(pair_id %in% common_edges)
edge_annot <- edge_annot %>% filter(pair_id %in% common_edges)

# enforce identical ordering
edge_annot <- edge_annot[match(diff_res$pair_id, edge_annot$pair_id), ]

stopifnot(all(edge_annot$pair_id == diff_res$pair_id))

cat("Aligned edges:", nrow(diff_res), "\n")

# ==============================
# FINAL MERGE (SAFE)
# ==============================


final_annot <- bind_cols(diff_res, edge_annot %>% select(-pair_id))

cat("Final rows:", nrow(final_annot), "\n")


cat("\nCleaning NA values...\n")

final_annot <- final_annot %>%
  mutate(
    dominant_MP   = ifelse(is.na(dominant_MP), "None", dominant_MP),
    SE_flag       = ifelse(is.na(SE_flag), FALSE, SE_flag),
    HiChIP        = ifelse(is.na(HiChIP), 0, HiChIP),
    crispr_edge   = ifelse(is.na(crispr_edge), FALSE, crispr_edge),
    MP_max        = ifelse(is.na(MP_max), 0, MP_max)
  )

wilcox.test(MP_max ~ SE_flag, data = final_annot)
  
# ==============================
# DEFINE REWIRED EDGES
# ==============================
cat("\nDefining rewired edges...\n")

final_annot$is_rewired <- final_annot$adj.P.Val < 0.05

cat("Rewired edges:", sum(final_annot$is_rewired), "\n")

cat("\nState-specific enrichment analysis...\n")

state_results <- lapply(state_list, function(st) {

  df <- final_annot %>% filter(signature == st)

  if (nrow(df) < 1000) return(NULL)

  safe_fisher <- function(x, y) {

    tab <- table(x, y)

    # FORCE 2×2 completeness
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      return(list(p.value = NA, estimate = NA))
    }

    ft <- fisher.test(tab)

    list(p.value = ft$p.value, estimate = unname(ft$estimate))
  }

  se_res <- safe_fisher(df$SE_flag, df$is_rewired)
  hic_res <- safe_fisher(df$HiChIP, df$is_rewired)
  cr_res <- safe_fisher(df$crispr_edge, df$is_rewired)

  data.frame(
    signature = st,

    SE_p = se_res$p.value,
    SE_OR = se_res$estimate,

    HiChIP_p = hic_res$p.value,
    HiChIP_OR = hic_res$estimate,

    CRISPR_p = cr_res$p.value,
    CRISPR_OR = cr_res$estimate,

    n_edges = nrow(df),
    n_rewired = sum(df$is_rewired)
  )
})

state_results <- bind_rows(state_results)

write.csv(state_results,
          file.path(out_dir, "STATE_SPECIFIC_ENRICHMENT.csv"),
          row.names = FALSE)

print(state_results)

cat("\nState interaction model...\n")

interaction_model <- glm(
  is_rewired ~ SE_flag * signature,
  data = final_annot,
  family = binomial()
)

anova_res <- anova(interaction_model, test = "Chisq")

write.csv(as.data.frame(anova_res),
          file.path(out_dir, "STATE_INTERACTION_TEST.csv"))
		  
		  

# ==============================
# MP ENRICHMENT
# ==============================
cat("\nMP enrichment...\n")

mp_enrichment <- final_annot %>%
  filter(is_rewired) %>%
  group_by(signature, dominant_MP) %>%
  summarise(
    mean_logFC = mean(logFC, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(signature, desc(mean_logFC))

write.csv(mp_enrichment,
          file.path(out_dir, "MP_enrichment_rewired_edges.csv"),
          row.names = FALSE)

# ==============================
# FISHER TESTS
# ==============================
cat("\nFisher tests...\n")

semp_fisher   <- fisher.test(table(final_annot$SE_flag, final_annot$is_rewired))
hic_fisher    <- fisher.test(table(final_annot$HiChIP, final_annot$is_rewired))
crispr_fisher <- fisher.test(table(final_annot$crispr_edge, final_annot$is_rewired))

cat("SE–MP P:", semp_fisher$p.value, "\n")
cat("HiChIP P:", hic_fisher$p.value, "\n")
cat("CRISPR P:", crispr_fisher$p.value, "\n")

# ==============================
# MP DIRECTIONALITY
# ==============================
cat("\nMP directionality...\n")

mp_direction <- final_annot %>%
  group_by(signature, dominant_MP) %>%
  summarise(
    mean_logFC = mean(logFC, na.rm = TRUE),
    frac_positive = mean(logFC > 0, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(signature, desc(mean_logFC))

write.csv(mp_direction,
          file.path(out_dir, "MP_directionality_all_edges.csv"),
          row.names = FALSE)
		  
		  
# ==============================
# SAVE FINAL OBJECT
# ==============================
saveRDS(final_annot,
        file.path(out_dir, "FINAL_INTEGRATED_EDGE_MODEL.rds"))
		

# ==============================
# PERMUTATION TESTS
# ==============================
cat("\nPermutation tests...\n")

perm_test <- function(feature_vec, label_vec, n_perm = 1000) {

  feature_vec <- as.numeric(feature_vec)
  label_vec   <- as.logical(label_vec)

  obs <- mean(feature_vec[label_vec], na.rm = TRUE)

  perm <- replicate(n_perm, {
    shuf <- sample(label_vec)
    mean(feature_vec[shuf], na.rm = TRUE)
  })

  mean(perm >= obs)
}

perm_p_semp <- perm_test(final_annot$SE_flag, final_annot$is_rewired)
perm_p_hic  <- perm_test(final_annot$HiChIP, final_annot$is_rewired)
perm_p_cr   <- perm_test(final_annot$crispr_edge, final_annot$is_rewired)

cat("Permutation SE–MP P:", perm_p_semp, "\n")
cat("Permutation HiChIP P:", perm_p_hic, "\n")
cat("Permutation CRISPR P:", perm_p_cr, "\n")

# ==============================
# DONE
# ==============================
cat("\n====================================\n")
cat("FULL INTEGRATION COMPLETE\n")
cat("SE–MP P:", semp_fisher$p.value, "\n")
cat("HiChIP P:", hic_fisher$p.value, "\n")
cat("CRISPR P:", crispr_fisher$p.value, "\n")
cat("Permutation SE–MP:", perm_p_semp, "\n")
cat("====================================\n")



crispr_hub <- edge_table %>% filter(crispr_edge == TRUE)

nrow(crispr_hub)


t.test(
  final_annot$regulatory_prior[final_annot$crispr_edge],
  final_annot$regulatory_prior[!final_annot$crispr_edge]
)

tab <- table(final_annot$crispr_edge, final_annot$SE_flag)
fisher.test(tab)

mean(final_annot$regulatory_prior[final_annot$crispr_edge])
mean(final_annot$regulatory_prior[!final_annot$crispr_edge])


final_annot <- diff_res %>%
  inner_join(edge_table, by = "pair_id")
  
  final_annot$is_rewired <- final_annot$adj.P.Val < 0.05
  
  table(final_annot$is_rewired, useNA = "ifany")
table(final_annot$HiChIP, useNA = "ifany")

tab <- table(final_annot$HiChIP, final_annot$is_rewired)
fisher.test(tab)



tab_global <- table(final_annot$HiChIP, final_annot$is_rewired)

ft_global <- fisher.test(tab_global)

OR_global <- as.numeric(ft_global$estimate)
CI_global <- ft_global$conf.int
P_global  <- ft_global$p.value

OR_global
CI_global
P_global

#################################################

mp_list <- unique(final_annot$dominant_MP)

mp_res <- lapply(mp_list, function(mp) {

  df <- final_annot %>%
    filter(dominant_MP == mp)

  # skip small MPs
  if (nrow(df) < 50) return(NULL)

  tab <- table(df$HiChIP, df$is_rewired)

  if (any(dim(tab) < 2)) return(NULL)

  ft <- fisher.test(tab)

  data.frame(
    MP = mp,
    OR = as.numeric(ft$estimate),
    CI_low = ft$conf.int[1],
    CI_high = ft$conf.int[2],
    p = ft$p.value,
    n = nrow(df)
  )
})

mp_res <- bind_rows(mp_res) %>%
  arrange(desc(OR))
final_annot %>%
  group_by(dominant_MP) %>%
  summarise(
    OR = {
      tab <- table(HiChIP, is_rewired)
      if (all(dim(tab) == c(2,2))) {
        fisher.test(tab)$estimate
      } else {
        NA
      }
    },
    p = {
      tab <- table(HiChIP, is_rewired)
      if (all(dim(tab) == c(2,2))) {
        fisher.test(tab)$p.value
      } else {
        NA
      }
    }
  )
  
  
####################################################################################

