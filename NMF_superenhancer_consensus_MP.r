
library(tidyverse)
library(DESeq2)
library(tibble)
library(NMF)  

##Super enhancer library size normalization
read.mat <- read.table("/home/user/genomes/se_counts1.txt", sep="\t", header=TRUE)


df <- read.mat %>% 
  select(-Chr, -Start, -End, -Strand, -Length) %>% 
  column_to_rownames("Geneid")


df <- round(df)

metadata <- read.table("/home/user/Documents/Milan/Epigenomics/SE_metadata.txt",
                       sep="\t", header=TRUE)


dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = metadata,
                              design = ~ 1)  # no design needed
dds <- estimateSizeFactors(dds)  # normalizes library size

# Get normalized counts
norm_counts <- counts(dds, normalized=TRUE)
rpk_log <- log2(norm_counts + 1)

write.table(rpk_log, file = "/home/user/Documents/Milan/Epigenomics/SE_normalized.txt", sep = "\t", quote = FALSE)

###########################################################
### Evaluate NMF quality metrics
library(NMF)
library(Matrix)
library(stats)
library(cluster)
library(igraph)
library(dplyr)
library(ggplot2)


# Read normalized log2 RPK matrix
rpk_log <- read.table(
  "/home/user/Documents/Milan/Epigenomics/SE_normalized.txt",
  header = TRUE,       # first row has column names (samples)
  row.names = 1,       # first column has SE identifiers
  sep = "\t",          # tab-delimited file
  check.names = FALSE  # keep column names exactly as they appear
)

sum(is.na(rpk_log))

rpk_log_mat <- as.matrix(rpk_log)
storage.mode(rpk_log_mat) <- "numeric"


rpk_log_mat <- rpk_log_mat[rowSums(rpk_log_mat, na.rm = TRUE) > 0, ]
rpk_log_mat <- rpk_log_mat[complete.cases(rpk_log_mat), ]

head(rownames(rpk_log_mat))
head(colnames(rpk_log_mat))

summary(as.vector(rpk_log_mat))
hist(as.vector(rpk_log_mat), breaks = 100, main = "Distribution of log-normalized SE counts", xlab = "log2 normalized count")
boxplot(rpk_log_mat[, 1:25], main="First 5 samples", las=2)

set.seed(12345)

rank_range <- 4:12  
nmf_rank_eval <- nmfEstimateRank(rpk_log_mat, range = rank_range, method = "brunet", nrun = 10, seed = 12345)

# visualize cophenetic and residual errors
png("estimateranks_nmf_programs.png",
    width = 3000, height = 3000, res = 300)
    
    plot(nmf_rank_eval)

dev.off()


plot(nmf_rank_eval)



set.seed(123)

# Evaluate NMF quality metrics (cophenetic, dispersion, residuals)
rank_range <- 4:12
nmf_rank_eval <- nmfEstimateRank(rpk_log, range = rank_range, method = "brunet", nrun = 10, seed = 123456)

# Visualize the stability of ranks
plot(nmf_rank_eval)


###################################################
##NMF decomposition seperately for each rank
library(NMF)
library(FNN)
library(igraph)
library(matrixStats)
#library(pheatmap)
library(doParallel)
library(foreach)

# ==============================
# Parallelized Subsampled NMF
# ==============================


# ------------------------------
#Setup parameters
# ------------------------------

set.seed(2025)
k_best   <- 10          # chosen rank
n_runs   <- 100        # number of subsampled runs
frac_SE  <- 0.7        # fraction of SEs per subsample
frac_smp <- 1          # fraction of samples per subsample
nrun_sub <- 3          # internal NMF repeats per subsample


SE_ids <- rownames(rpk_log_mat)
S <- nrow(rpk_log_mat)


n_cores <- parallel::detectCores()
use_cores <- max(1, n_cores - 2)
cl <- makeCluster(use_cores)
registerDoParallel(cl)
nmf.options(parallel = use_cores)
parallel::clusterSetRNGStream(cl, iseed = 2025)

cat("Using", use_cores, "cores for parallel NMF runs\n")

# ------------------------------
#Parallelized NMF subsampling
# ------------------------------
results <- foreach(r = seq_len(n_runs), .packages = "NMF") %dopar% {
  
  cat("Subsample run:", r, "/", n_runs, "\n")
  
  # subsample SEs and samples
  sel_SE <- sample(SE_ids, size = round(frac_SE * S), replace = FALSE)
  sel_samples <- sample(colnames(rpk_log_mat), size = round(frac_smp * ncol(rpk_log_mat)), replace = FALSE)
  submat <- rpk_log_mat[sel_SE, sel_samples]

  # ensure non-negativity
  minv <- min(submat)
  if (minv < 0) submat <- submat - minv + 1e-6

  # run NMF (brunet method)
  nmf_sub <- nmf(submat, rank = k_best, method = "brunet", nrun = nrun_sub, .opt = "v")
  W_sub <- basis(nmf_sub)  # rows = SEs, cols = components
  H_sub <- coef(nmf_sub)   # rows = components, cols = samples


  W_full <- matrix(0, nrow = nrow(rpk_log_mat), ncol = ncol(W_sub),
                   dimnames = list(rownames(rpk_log_mat),
                                   paste0("run", r, ".c", seq_len(ncol(W_sub)))))
  W_full[sel_SE, ] <- W_sub

 
  H_full <- matrix(0, nrow = nrow(H_sub), ncol = ncol(rpk_log_mat),
                   dimnames = list(paste0("run", r, ".c", seq_len(nrow(H_sub))),
                                   colnames(rpk_log_mat)))
  H_full[, sel_samples] <- H_sub
  
  list(W = W_full, H = H_full)
}


W_list <- lapply(results, `[[`, "W")
H_list <- lapply(results, `[[`, "H")


stopCluster(cl)
registerDoSEQ()

saveRDS(W_list, file = "/home/user/Documents/Milan/Epigenomics/NMF_Wlist_rank10_parallel.rds")
saveRDS(H_list, file = "/home/user/Documents/Milan/Epigenomics/NMF_Hlist_rank10_parallel.rds")

##the NMF decomposition was done for the ranks: 6,7,8,9,10; only codes for rank 10 is shown here.

##building robust NMF consensus meta programs for ~20 MP
###############################################################################

library(NMF)
library(Matrix)
library(stats)
library(cluster)
library(igraph)
library(dplyr)
library(ggplot2)

rpk_log <- read.table(
  "/home/user/Documents/Milan/Epigenomics/SE_normalized.txt",
  header = TRUE,       
  row.names = 1,       
  sep = "\t",          
  check.names = FALSE  
)

sum(is.na(rpk_log))

rpk_log_mat <- as.matrix(rpk_log)
storage.mode(rpk_log_mat) <- "numeric"

rpk_log_mat <- rpk_log_mat[rowSums(rpk_log_mat, na.rm = TRUE) > 0, ]
rpk_log_mat <- rpk_log_mat[complete.cases(rpk_log_mat), ]


# ======================================================
# Hybrid NNLS × Jaccard meta-program pipeline
# (0.6 * Jaccard + 0.4 * PearsonCorr)
# ======================================================

library(dplyr)
library(Matrix)
library(igraph)
library(leidenbase)
library(reticulate)
library(matrixStats)
library(pheatmap)
library(viridis)

# --- Python NNLS setup ---
np <- import("numpy")
sp <- import("scipy.optimize")

base_dir <- "/home/user/Documents/Milan/Epigenomics"
W_files <- list(
  rank6  = file.path(base_dir, "NMF_Wlist_rank6_parallel.rds"),
  rank7  = file.path(base_dir, "NMF_Wlist_rank7_parallel.rds"),
  rank8  = file.path(base_dir, "NMF_Wlist_rank8_parallel.rds"),
  rank9  = file.path(base_dir, "NMF_Wlist_rank9_parallel.rds"),
  rank10 = file.path(base_dir, "NMF_Wlist_rank10_parallel.rds")
)
H_files <- list(
  rank6  = file.path(base_dir, "NMF_Hlist_rank6_parallel.rds"),
  rank7  = file.path(base_dir, "NMF_Hlist_rank7_parallel.rds"),
  rank8  = file.path(base_dir, "NMF_Hlist_rank8_parallel.rds"),
  rank9  = file.path(base_dir, "NMF_Hlist_rank9_parallel.rds"),
  rank10 = file.path(base_dir, "NMF_Hlist_rank10_parallel.rds")
)

# ======================================================
# Load NMF W/H results
# ======================================================
all_Wlist <- list()
all_Hlist <- list()

for (rn in names(W_files)) {
  wf <- W_files[[rn]]
  hf <- H_files[[rn]]
  if (file.exists(wf) && file.exists(hf)) {
    message("Loading ", rn, " ...")
    all_Wlist[[rn]] <- readRDS(wf)
    all_Hlist[[rn]] <- readRDS(hf)
  } else {
    message("Missing: ", rn, " -> skipping.")
  }
}
if (length(all_Wlist) == 0) stop("No W/H lists found!")

# ======================================================
# Collect programs and NNLS activity
# ======================================================

top_n_SE <- 1500
topn <- se_names[W[, comp] > quantile(W[, comp], 0.9)]
programs <- list()
program_activity <- list()

for (rname in names(all_Wlist)) {
  Wlist <- all_Wlist[[rname]]
  Hlist <- all_Hlist[[rname]]
  for (i in seq_along(Wlist)) {
    W <- Wlist[[i]]
    H <- Hlist[[i]]
    if (is.null(dim(W)) || ncol(W) < 1) next
    for (comp in seq_len(ncol(W))) {
      prog_id <- paste0(rname, "_run", i, "_c", comp)
      se_names <- rownames(W)
      ord <- order(W[, comp], decreasing = TRUE, na.last = TRUE)
      topn <- head(se_names[ord], min(top_n_SE, length(se_names)))
      programs[[prog_id]] <- topn
      if (!is.null(dim(H))) {
        program_activity[[prog_id]] <- as.numeric(H[comp, ])
      } else {
        program_activity[[prog_id]] <- as.numeric(H)
      }
    }
  }
}
prog_names <- names(programs)
n_prog <- length(programs)
message("Total programs collected: ", n_prog)
if (n_prog < 3) stop("Too few programs found.")

# ======================================================
# Compute Jaccard similarity
# ======================================================

library(future)
library(future.apply)
plan(multisession, workers = 32)


prog_names <- names(programs)
n_prog <- length(programs)
message("Computing Jaccard similarity for ", n_prog, " programs...")


jaccard_mat <- matrix(0, nrow = n_prog, ncol = n_prog,
                      dimnames = list(prog_names, prog_names))


compute_jaccard_row <- function(i) {
  ps_i <- programs[[i]]
  res <- numeric(n_prog - i)
  for (j in seq((i + 1), n_prog)) {
    ps_j <- programs[[j]]
    inter <- length(intersect(ps_i, ps_j))
    uni <- length(union(ps_i, ps_j))
    res[j - i] <- ifelse(uni > 0, inter / uni, 0)
  }
  res
}


res_list <- future_lapply(seq_len(n_prog - 1), compute_jaccard_row,
                          future.seed = TRUE)


for (i in seq_len(n_prog - 1)) {
  vals <- res_list[[i]]
  j_idx <- seq((i + 1), n_prog)
  jaccard_mat[i, j_idx] <- vals
  jaccard_mat[j_idx, i] <- vals
}
diag(jaccard_mat) <- 1

message("Jaccard similarity matrix computed successfully.")

# ------------------------------
#Compute program activity correlation (Pearson, vectorized)
# ------------------------------

library(future)
library(future.apply)

plan(multisession, workers = 32)

activity_matrix <- do.call(rbind, program_activity)
rownames(activity_matrix) <- prog_names


activity_matrix <- as.matrix(activity_matrix)
activity_matrix[is.na(activity_matrix)] <- 0

message("Computing Pearson correlation matrix across ", n_prog, " programs...")


corr_mat <- cor(t(activity_matrix), method = "pearson", use = "pairwise.complete.obs")


corr_mat[is.na(corr_mat)] <- 0
diag(corr_mat) <- 1

# Scale correlation values to [0, 1] range for combination with Jaccard
corr_mat_scaled <- (corr_mat + 1) / 2


# ------------------------------
#Combine structural + functional similarity
# ------------------------------

# Both matrices are [0, 1] range: Jaccard (structural) + scaled correlation (functional)
# Adjust weights if needed: 0.6 = structure, 0.4 = function
combined_mat <- 0.6 * jaccard_mat + 0.4 * corr_mat_scaled

# Set self-similarity to 0 before network building
diag(combined_mat) <- 0

message("Combined structural (Jaccard) and functional (correlation) similarity matrices computed.")

# ==========================================================
#Build network and adaptive Leiden clustering
# ==========================================================
library(igraph)

thresh <- quantile(combined_mat[combined_mat > 0], 0.7)
adj_mat <- combined_mat
adj_mat[adj_mat < thresh] <- 0

g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE)
E(g)$weight <- E(g)$weight

cat("Network constructed with", gorder(g), "nodes and", gsize(g), "edges\n")

set.seed(2025)
target_clusters <- 20
resolution <- 1.2

for (iter in 1:20) {
  cl <- cluster_leiden(
    g,
    objective_function = "CPM",
    resolution = resolution,
    weights = E(g)$weight
  )

  n_clust <- length(unique(membership(cl)))
  cat("Iteration", iter, ": resolution =", round(resolution, 6), "->", n_clust, "clusters\n")
  if (n_clust <= target_clusters) break
  resolution <- resolution * 0.8
}

membership_vec <- membership(cl)
names(membership_vec) <- names(programs)
cat("Final number of consensus meta-programs:", length(unique(membership_vec)), "\n")

# ==========================================================
#Build consensus meta-programs (aggregate SEs)
# ==========================================================

library(reticulate)
library(Matrix)

meta_programs <- list()

cat("Aggregating SEs for consensus meta-programs...\n")

for (cid in unique(membership_vec)) {
  members <- names(membership_vec[membership_vec == cid])
  merged_SEs <- unlist(programs[members], use.names = FALSE)
  freq_table <- sort(table(merged_SEs), decreasing = TRUE)
  
  
  top_consensus_SEs <- names(freq_table)[1:min(1500, length(freq_table))]
  
  meta_programs[[paste0("meta_", cid)]] <- top_consensus_SEs
}

cat("Consensus meta-programs generated:", length(meta_programs), "\n")


all_SEs <- unique(unlist(meta_programs))
W_meta <- matrix(0, nrow = length(all_SEs), ncol = length(meta_programs),
                 dimnames = list(all_SEs, names(meta_programs)))

cat("Building binary SE × meta-program matrix...\n")

for (m in names(meta_programs)) {
  W_meta[meta_programs[[m]], m] <- 1
}


V <- rpk_log_mat  

common_SEs <- intersect(rownames(V), rownames(W_meta))
cat("Common SEs between meta-programs and expression matrix:", length(common_SEs), "\n")

V <- V[common_SEs, , drop = FALSE]
W_meta <- W_meta[common_SEs, , drop = FALSE]

# ==========================================================
#Run NNLS for meta-program activity
# ==========================================================

cat("Running NNLS to infer meta-program activity...\n")


np <- import("numpy", convert = TRUE)
sp <- import("scipy.optimize", convert = TRUE)

# Ensure numeric conversion
W_meta_np <- np$array(as.matrix(W_meta))  

H_meta <- matrix(0, nrow = ncol(W_meta), ncol = ncol(V),
                 dimnames = list(colnames(W_meta), colnames(V)))

pb <- txtProgressBar(min = 0, max = ncol(V), style = 3)

for (i in seq_len(ncol(V))) {
  y_vector <- np$array(as.numeric(V[, i]))
  res <- sp$nnls(W_meta_np, y_vector)
  H_meta[, i] <- as.numeric(res[[1]])
  setTxtProgressBar(pb, i)
}
close(pb)

write.csv(H_meta, "H_meta_program_activity.csv")
cat("\nSaved NNLS meta-program activity to 'H_meta_program_activity.csv'\n")

# ==========================================================
#Save all outputs
# ==========================================================

save_path <- "/home/user/Documents/Milan/Epigenomics/meta_programs_rank6_8_consensus_SE+activity_leiden_NNLS.rds"

saveRDS(
  list(
    jaccard_matrix = jaccard_mat,
    corr_matrix = corr_mat,
    combined_matrix = combined_mat,
    network = g,
    membership = membership_vec,
    meta_programs = meta_programs,
    H_meta = H_meta,
    hc_meta = if (exists("hc_meta")) hc_meta else NULL,
    resolution = resolution,
    threshold = thresh
  ),
  file = save_path
)

cat("Saved integrated meta-programs (structural + functional + NNLS activity):\n",
    save_path, "\n")


#######################################################################################

# Extract SEs for each meta-program and their NNLS activity


# ------------------------------
#Save SEs in each meta-program
# ------------------------------
meta_SEs_df <- data.frame(
  meta_program = rep(names(meta_programs), times = sapply(meta_programs, length)),
  SE = unlist(meta_programs, use.names = FALSE)
)

write.csv(meta_SEs_df, "meta_program_SEs.csv", row.names = FALSE)
cat("Saved SEs in each meta-program → meta_program_SEs.csv\n")

# ------------------------------
#Extract meta-program activity for first 25 samples
# ------------------------------

if (nrow(H_meta) < ncol(H_meta)) {
  H_25 <- H_meta[, seq_len(min(25, ncol(H_meta))), drop = FALSE]
} else {
  H_25 <- H_meta[seq_len(min(25, nrow(H_meta))), , drop = FALSE]
}

colnames(H_25) <- paste0("Sample_", seq_len(ncol(H_25)))


library(reshape2)
H_25_long <- melt(H_25, varnames = c("Meta_program", "Sample"), value.name = "Activity")

write.csv(H_25, "H_meta_25samples.csv")
write.csv(H_25_long, "H_meta_25samples_long.csv", row.names = FALSE)

cat("Saved NNLS meta-program activity for 25 samples:\n   → H_meta_25samples.csv\n   → H_meta_25samples_long.csv\n")

# ==========================================================
#Compute SE-level activity within each meta-program using NNLS
# ==========================================================
library(reticulate)

np <- import("numpy", convert = FALSE)
sp <- import("scipy.optimize", convert = FALSE)

SE_activity_all <- list()

cat("\n Running per-meta-program NNLS for SE-level activity...\n")

for (mp_name in names(meta_programs)) {
  cat("Processing meta-program:", mp_name, "\n")
  
  SEs <- meta_programs[[mp_name]]
  SEs <- intersect(SEs, rownames(rpk_log_mat))  
  
  if (length(SEs) < 2) {
    cat("     Skipping", mp_name, "(<2 SEs)\n")
    next
  }
  
  V_sub <- rpk_log_mat[SEs, , drop = FALSE]
  n_SE <- nrow(V_sub)
  n_samp <- ncol(V_sub)
  
  
  W_sub <- diag(n_SE)
  rownames(W_sub) <- SEs
  colnames(W_sub) <- SEs
  
  H_SE <- matrix(0, nrow = n_SE, ncol = n_samp,
                 dimnames = list(SEs, colnames(V_sub)))
  
  for (i in seq_len(n_samp)) {
    y_vector <- np$array(as.numeric(V_sub[, i]))
    
    
    res <- sp$nnls(np$array(W_sub), np$array(y_vector))
res_r <- py_to_r(res)  
H_SE[, i] <- as.numeric(res_r[[1]])
  }
  
  SE_activity_all[[mp_name]] <- H_SE
}

saveRDS(SE_activity_all, file = "SE_activity_per_meta_program_NNLS.rds")
cat("Saved NNLS-based SE-level activity per meta-program → SE_activity_per_meta_program_NNLS.rds\n")

######################################################################################################
### heatmap: SE vs samples 


library(dplyr)
library(matrixStats)
library(pheatmap)
library(RColorBrewer)
library(viridis)

# ------------------------------
#Load NNLS-based SE activity matrices
# ------------------------------
SE_activity_all <- readRDS("SE_activity_per_meta_program_NNLS.rds")

# Combine all SEs
all_SEs <- unique(unlist(lapply(SE_activity_all, rownames)))
sample_names <- colnames(SE_activity_all[[1]])

# ------------------------------
#Build combined SE × Sample matrix
# ------------------------------
SE_activity_combined <- matrix(0, nrow = length(all_SEs), ncol = length(sample_names),
                               dimnames = list(all_SEs, sample_names))

for (mp_name in names(SE_activity_all)) {
  mat <- SE_activity_all[[mp_name]]
  overlap <- intersect(rownames(mat), all_SEs)
  SE_activity_combined[overlap, ] <- mat[overlap, ]
}

# ------------------------------
#Assign SEs to their meta-program
# ------------------------------
SE_to_meta <- sapply(rownames(SE_activity_combined), function(se) {
  mp_name <- names(SE_activity_all)[sapply(SE_activity_all, function(x) se %in% rownames(x))]
  mp_name[1]
})

annotation_row <- data.frame(Meta_program = factor(SE_to_meta))
rownames(annotation_row) <- rownames(SE_activity_combined)

# ------------------------------
#Normalize (Z-score per SE)
# ------------------------------
SE_activity_scaled <- t(scale(t(SE_activity_combined)))

# ------------------------------
#Order SEs by meta-program
# ------------------------------
ordered_SEs <- rownames(annotation_row)[order(annotation_row$Meta_program)]
annotation_row_ordered <- annotation_row[ordered_SEs, , drop = FALSE]
SE_activity_ordered <- SE_activity_scaled[ordered_SEs, , drop = FALSE]


num_MPs <- length(unique(SE_to_meta))


base_colors <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33",
                 "#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB")
meta_colors <- colorRampPalette(base_colors)(num_MPs)
names(meta_colors) <- unique(SE_to_meta)
ann_colors <- list(Meta_program = meta_colors)


activity_palette <- colorRampPalette(c("#FDFDFD", "#FFDAB9", "#CC5500", "#4B0082"))(100)


activity_breaks <- unique(round(seq(-2.5, 2.5, length.out = 101), 3))


gaps <- cumsum(table(annotation_row_ordered$Meta_program))
gaps <- gaps[-length(gaps)] 


png("SExSample_NNLS_activity_heatmap_light_background1.png",
    width = 5000, height = 4000, res = 300)

pheatmap(SE_activity_ordered,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         clustering_method = "ward.D2",
         gaps_row = gaps,
         color = activity_palette,
         breaks = activity_breaks,
         annotation_row = annotation_row_ordered,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         main = "Super-Enhancer Activity per Sample (Light Background + Dark Signal + Top Magenta)",
         fontsize_col = 8,
         border_color = NA)

dev.off()



# ===============================
# SE x SE Correlation Heatmap 
# ===============================
library(dplyr)
library(matrixStats)
library(pheatmap)
library(RColorBrewer)
library(viridis)

# ------------------------------
#Load NNLS-based SE activity matrices
# ------------------------------
SE_activity_all <- readRDS("SE_activity_per_meta_program_NNLS.rds")

# ------------------------------
#Union of all SEs across programs
# ------------------------------
all_SEs <- unique(unlist(lapply(SE_activity_all, rownames)))
prog_names <- names(SE_activity_all)
n_prog <- length(prog_names)

# ------------------------------
#Build SE × MP activity matrix (mean across samples)
# ------------------------------
SE_activity_mat <- matrix(0, nrow = length(all_SEs), ncol = n_prog,
                          dimnames = list(all_SEs, prog_names))

for (p in seq_along(SE_activity_all)) {
  mat <- SE_activity_all[[p]]
  overlap <- intersect(all_SEs, rownames(mat))
  SE_activity_mat[overlap, p] <- rowMeans(mat[overlap, , drop = FALSE])
}

# ------------------------------
#Z-score normalize per SE
# ------------------------------
SE_activity_mat <- t(scale(t(SE_activity_mat)))

# ------------------------------
#Assign SEs to their dominant meta-program
# ------------------------------
SE_to_meta <- apply(SE_activity_mat, 1, function(x) {
  names(which.max(x))
})
annotation_row <- data.frame(Meta_program = factor(SE_to_meta))
rownames(annotation_row) <- rownames(SE_activity_mat)

# ------------------------------
#Order SEs by Meta-Program
# ------------------------------
ordered_SEs <- rownames(annotation_row)[order(annotation_row$Meta_program)]
annotation_row_ordered <- annotation_row[ordered_SEs, , drop = FALSE]

# ------------------------------
#Compute SE × SE correlation
# ------------------------------
SE_corr <- cor(t(as.matrix(SE_activity_mat)), method = "pearson")
SE_corr_ordered <- SE_corr[ordered_SEs, ordered_SEs]


num_MPs <- length(unique(SE_to_meta))


base_colors <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33",
                 "#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB")
meta_colors <- colorRampPalette(base_colors)(num_MPs)
names(meta_colors) <- unique(SE_to_meta)
ann_colors <- list(Meta_program = meta_colors)


custom_very_dark_palette <- colorRampPalette(c("#FFFFFF", "#CC5500", "#4B0082"))(100)


breaks <- unique(round(c(
  seq(-1, -0.3, length.out = 10),
  seq(-0.3, 0.3, length.out = 20),
  seq(0.3, 1, length.out = 70)
), 3))


png("SExSE_correlation_heatmap_NNLS_MP_annotation_very_dark1.png",
    width = 4000, height = 4000, res = 300)

pheatmap(SE_corr_ordered,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = custom_very_dark_palette,
         breaks = breaks,
         annotation_row = annotation_row_ordered,
         annotation_col = annotation_row_ordered,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Super-Enhancer × Super-Enhancer Correlation (Bright Background + Very Dark Signal + Top Magenta)",
         border_color = NA)

dev.off()



