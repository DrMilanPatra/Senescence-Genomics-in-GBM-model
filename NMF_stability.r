########################################################################
# Computes NMF stability statistics
########################################################################


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(NMF)
  library(future.apply)
  library(parallel)
  library(doParallel)
  library(igraph)
  library(cluster)
  library(clue)        # solve_LSAP (Hungarian)
  library(proxy)       
  library(ggplot2)
  library(pheatmap)
  library(readr)
})


top_n <- 750              
ranks_to_use <- 6:10      
wlist_pattern <- "/home/user/Documents/Milan/Epigenomics/NMF_Wlist_rank%d_parallel.rds"
hlist_pattern <- "/home/user/Documents/Milan/Epigenomics/NMF_Hlist_rank%d_parallel.rds"
rpk_path <- "/home/user/Documents/Milan/Epigenomics/SE_normalized.txt"

out_dir <- "/home/user/Documents/Milan/Epigenomics/NMF_stability_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_boot <- 500             
n_split <- 50             
nmf_nrun_split <- 5       
nmf_method <- "brunet"
use_cores <- max(1, parallel::detectCores() - 2)
registerDoParallel(cores = use_cores)
message("Using cores: ", use_cores)

set.seed(2025)

# ------------------------------
# Load rpk matrix (for split-half)
# ------------------------------
rpk_log <- read.table(rpk_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
rpk_log_mat <- as.matrix(rpk_log)
storage.mode(rpk_log_mat) <- "numeric"
rpk_log_mat <- rpk_log_mat[rowSums(is.na(rpk_log_mat)) == 0, , drop = FALSE]  


# cosine similarity
cos_sim <- function(a, b) {

  proxy::simil(a, b, method = "cosine", by_rows = FALSE)
}

# jaccard between two sets
jaccard_set <- function(a, b) {
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(0)
  inter / uni
}


extract_top_sets <- function(Wmat, n = top_n) {
  sets <- lapply(seq_len(ncol(Wmat)), function(c) {
    rn <- rownames(Wmat)
    ord <- order(Wmat[, c], decreasing = TRUE, na.last = TRUE)
    top <- rn[ord][seq_len(min(n, length(rn)))]
    top
  })
  names(sets) <- paste0("comp", seq_len(ncol(Wmat)))
  sets
}



build_consensus_matrix <- function(all_sets, SE_universe) {

  SE_universe <- unique(SE_universe)
  nSE <- length(SE_universe)
  name_to_index <- setNames(seq_len(nSE), SE_universe)
  C <- matrix(0, nrow = nSE, ncol = nSE, dimnames = list(SE_universe, SE_universe))
  total_runs <- length(all_sets)
  for (run_sets in all_sets) {

    for (comp_set in run_sets) {
      idx <- name_to_index[intersect(comp_set, SE_universe)]
      if (length(idx) > 0) {
       
        C[idx, idx] <- C[idx, idx] + 1
      }
    }
  }
  C / total_runs
}

# compute cophenetic coefficient from consensus matrix
cophenetic_from_consensus <- function(Cmat) {

  dist_mat <- as.dist(1 - Cmat)
  hc <- hclust(dist_mat, method = "average")
  # cophenetic distances from hc
  cop_d <- cophenetic(hc)
  # compare cophenetic distances to original distance
  cor(as.vector(dist_mat), as.vector(cop_d), method = "pearson")
}


partition_from_consensus <- function(Cmat, k) {
  dist_mat <- as.dist(1 - Cmat)
  hc <- hclust(dist_mat, method = "average")
  cutree(hc, k = k)
}

# mean silhouette width given Cmat and partition
mean_sil_from_consensus <- function(Cmat, partition) {
  d <- as.dist(1 - Cmat)
  sil <- silhouette(partition, d)
  mean(sil[, "sil_width"], na.rm = TRUE)
}

# match components between two W matrices using cosine similarity and Hungarian algorithm
match_W_matrices <- function(W_A, W_B) {
  
  common <- intersect(rownames(W_A), rownames(W_B))
  A <- W_A[common, , drop = FALSE]
  B <- W_B[common, , drop = FALSE]
  # compute similarity matrix between columns (components)
  sim_mat <- proxy::simil(A, B, method = "cosine", by_rows = FALSE) # returns matrix components_A x components_B
  sim_mat <- as.matrix(sim_mat)
  # Hungarian maximum assignment -> convert to cost (1 - similarity)
  cost <- 1 - sim_mat
  # Solve LSAP requires square matrix; pad if needed
  n <- max(ncol(sim_mat), nrow(sim_mat))
  cost_p <- matrix(1, nrow = n, ncol = n)
  cost_p[1:nrow(cost), 1:ncol(cost)] <- cost
  assignment <- solve_LSAP(cost_p)
  assignment <- as.integer(assignment)

  matched <- rep(NA, ncol(sim_mat))
  for (i in seq_len(ncol(sim_mat))) {
    j <- assignment[i]
    if (j <= ncol(sim_mat)) matched[i] <- sim_mat[i, j]
  }
  list(assignment = assignment, matched = matched, sim_mat = sim_mat)
}

# ------------------------------
# Read W_lists for each rank and collate run-level component sets
# ------------------------------
all_rank_results <- list()

for (k in ranks_to_use) {
  wf <- sprintf(wlist_pattern, k)
  if (!file.exists(wf)) {
    warning("Missing Wlist for rank ", k, " -> skipping")
    next
  }
  Wlist <- readRDS(wf)  # list of matrices (W_full padded)
  message("Rank ", k, ": loaded ", length(Wlist), " runs")
  # For each run, we may have matrix with colnames representing components
  run_sets <- list()
  run_Ws <- list()
  for (i in seq_along(Wlist)) {
    Wfull <- Wlist[[i]]
    # extract top_n top SEs per component (if a column is all zeros - skip)
    nonzero_cols <- which(colSums(Wfull, na.rm = TRUE) > 0)
    if (length(nonzero_cols) == 0) next
    Wsub <- Wfull[, nonzero_cols, drop = FALSE]
    sets <- extract_top_sets(Wsub, n = top_n)
    names(sets) <- paste0("run", i, ".comp", seq_along(sets))
    run_sets[[paste0("run", i)]] <- sets
    run_Ws[[paste0("run", i)]] <- Wsub
  }
  all_rank_results[[paste0("rank", k)]] <- list(run_sets = run_sets, run_Ws = run_Ws)
}

# ------------------------------
# PER-RANK STABILITY: Jaccard across runs, consensus matrix, cophenetic, silhouette
# ------------------------------
rank_summaries <- list()

for (rk in names(all_rank_results)) {
  message("Processing ", rk)
  run_sets <- all_rank_results[[rk]]$run_sets
  run_Ws <- all_rank_results[[rk]]$run_Ws
  if (length(run_sets) < 2) {
    warning("Not enough runs for ", rk)
    next
  }
 
  SE_univ <- unique(unlist(lapply(run_sets, function(x) unlist(x))))
  
  runs_list <- lapply(run_sets, function(x) x)
  
  # 1) Pairwise Jaccard between component instances (across runs)
  
  comp_names <- unlist(lapply(seq_along(runs_list), function(i) {
    paste0(names(runs_list)[i], "_", names(runs_list[[i]]))
  }))
  comp_sets <- unlist(runs_list, recursive = FALSE)
  names(comp_sets) <- comp_names
  
  n_comp <- length(comp_sets)
  jacc_mat <- matrix(0, nrow = n_comp, ncol = n_comp, dimnames = list(comp_names, comp_names))
  for (i in seq_len(n_comp)) {
    for (j in i:n_comp) {
      v <- jaccard_set(comp_sets[[i]], comp_sets[[j]])
      jacc_mat[i, j] <- v
      jacc_mat[j, i] <- v
    }
  }
  
  # 2) Cluster components into k groups (k = rank number)
  k_num <- as.integer(sub("rank", "", rk))
  # use hierarchical clustering on 1-jaccard as distance
  hc_comp <- hclust(as.dist(1 - jacc_mat), method = "average")
  comp_cluster <- cutree(hc_comp, k = k_num)
  
  # For each consensus cluster, compute mean within-cluster Jaccard (stability)
  cluster_ids <- sort(unique(comp_cluster))
  cluster_stability <- sapply(cluster_ids, function(cid) {
    members <- which(comp_cluster == cid)
    if (length(members) <= 1) return(NA_real_)
    mean(jacc_mat[members, members][lower.tri(jacc_mat[members, members])], na.rm = TRUE)
  })
  names(cluster_stability) <- paste0("cons_comp_", cluster_ids)
  
  # 3) Consensus co-association matrix across runs (SE x SE)
  Cmat <- build_consensus_matrix(runs_list, SE_univ)
  
  # 4) Cophenetic coefficient from consensus
  cop_coef <- cophenetic_from_consensus(Cmat)
  
  # 5) Partition of SEs into k clusters from consensus and silhouette
  partition <- partition_from_consensus(Cmat, k = k_num)
  mean_sil <- mean_sil_from_consensus(Cmat, partition)
  
 
  rank_summaries[[rk]] <- list(
    jaccard_matrix = jacc_mat,
    comp_cluster = comp_cluster,
    cluster_stability = cluster_stability,
    consensus_matrix = Cmat,
    cophenetic = cop_coef,
    partition = partition,
    mean_silhouette = mean_sil,
    SE_universe = SE_univ,
    runs_list = runs_list,
    run_Ws = run_Ws
  )
  
 
  png(file.path(out_dir, paste0(rk, "_component_jaccard_heatmap.png")), width = 2000, height = 1800, res = 200)
  pheatmap(jacc_mat, main = paste0("Component Jaccard (", rk, ")"), cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE)
  dev.off()
  
  png(file.path(out_dir, paste0(rk, "_consensus_matrix_heatmap.png")), width = 2000, height = 2000, res = 200)
  pheatmap(Cmat, main = paste0("Consensus co-association (", rk, ")"), show_rownames = FALSE, show_colnames = FALSE)
  dev.off()
}

# ------------------------------
# BOOTSTRAP CI for cophenetic and silhouette per rank (resample runs)
# ------------------------------
message("Bootstrap CIs (", n_boot, " iterations) - this may take a while")

boot_stats <- list()
for (rk in names(rank_summaries)) {
  info <- rank_summaries[[rk]]
  runs_list <- info$runs_list
  n_runs_avail <- length(runs_list)
  if (n_runs_avail < 2) next
  SE_univ <- info$SE_universe
  k_num <- as.integer(sub("rank", "", rk))
  
  # function to compute stats for a bootstrap sample of runs (indices with replacement)
  compute_stats_from_indices <- function(idx_vec) {
    sampled_runs <- runs_list[idx_vec]
    Cb <- build_consensus_matrix(sampled_runs, SE_univ)
    cop_b <- cophenetic_from_consensus(Cb)
    partition_b <- partition_from_consensus(Cb, k = k_num)
    sil_b <- mean_sil_from_consensus(Cb, partition_b)
    # cluster stability: compute jaccard across sampled components and aggregate per cluster (approx)
    comp_sets_b <- unlist(sampled_runs, recursive = FALSE)
    if (length(comp_sets_b) < 2) {
      cs_mean <- NA
    } else {
      # quick pairwise jaccard mean
      comp_names_b <- seq_along(comp_sets_b)
      j_t <- numeric()
      for (i in seq_len(length(comp_sets_b) - 1)) {
        for (j in (i + 1):length(comp_sets_b)) j_t <- c(j_t, jaccard_set(comp_sets_b[[i]], comp_sets_b[[j]]))
      }
      cs_mean <- mean(j_t, na.rm = TRUE)
    }
    c(cophenetic = cop_b, silhouette = sil_b, mean_comp_jaccard = cs_mean)
  }
  
  # bootstrap run: sample indices
  boot_res <- future_lapply(seq_len(n_boot), function(b) {
    idx <- sample(seq_len(n_runs_avail), size = n_runs_avail, replace = TRUE)
    compute_stats_from_indices(idx)
  }, future.seed = TRUE)
  
  boot_mat <- do.call(rbind, boot_res)
  colnames(boot_mat) <- c("cophenetic", "silhouette", "mean_comp_jaccard")
  # compute medians and 95% CI
  stat_summary <- data.frame(
    metric = c("cophenetic", "silhouette", "mean_comp_jaccard"),
    observed = c(info$cophenetic, info$mean_silhouette, mean(info$cluster_stability, na.rm = TRUE)),
    boot_mean = apply(boot_mat, 2, mean, na.rm = TRUE),
    boot_sd = apply(boot_mat, 2, sd, na.rm = TRUE),
    CI_low = apply(boot_mat, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
    CI_high = apply(boot_mat, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
  )
  
  boot_stats[[rk]] <- list(boot_mat = boot_mat, summary = stat_summary)
  
  write.csv(stat_summary, file.path(out_dir, paste0(rk, "_bootstrap_stat_summary.csv")), row.names = FALSE)
  message("Rank ", rk, " boot summary saved.")
}

# ------------------------------
# Split-half reproducibility (will run NMF on split halves)
# ------------------------------

# ------------------------------
# Split-half reproducibility (NMF on random half splits)
# ------------------------------
message("Split-half reproducibility: running ", n_split, 
        " random splits with ", nmf_nrun_split, 
        " NMF repeats per split (method = ", nmf_method, ")")

split_stats <- list()

for (rep_i in seq_len(n_split)) {

  # Randomly split samples
  samples <- colnames(rpk_log_mat)
  half1 <- sample(samples, size = ceiling(length(samples)/2), replace = FALSE)
  half2 <- setdiff(samples, half1)

  
  if (length(half2) < 2) {
    half2 <- sample(samples, size = floor(length(samples)/2), replace = FALSE)
    half1 <- setdiff(samples, half2)
  }

  for (k in ranks_to_use) {

  
    mat1 <- rpk_log_mat[, half1, drop = FALSE]
    mat2 <- rpk_log_mat[, half2, drop = FALSE]

  
    if (min(mat1) < 0) mat1 <- mat1 - min(mat1) + 1e-6
    if (min(mat2) < 0) mat2 <- mat2 - min(mat2) + 1e-6

    
    mat1 <- mat1[rowSums(mat1) > 0 & apply(mat1, 1, var) > 0, , drop = FALSE]
    mat2 <- mat2[rowSums(mat2) > 0 & apply(mat2, 1, var) > 0, , drop = FALSE]

   
    common_genes <- intersect(rownames(mat1), rownames(mat2))
    mat1 <- mat1[common_genes, , drop = FALSE]
    mat2 <- mat2[common_genes, , drop = FALSE]

  
    if (nrow(mat1) < k + 1 || nrow(mat2) < k + 1) {
      message("Skipping split ", rep_i, " rank ", k, 
              " due to insufficient genes after filtering")
      next
    }

    # Run small NMFs
    nmf1 <- nmf(mat1, rank = k, method = nmf_method,
                nrun = nmf_nrun_split, .options = "v")
    nmf2 <- nmf(mat2, rank = k, method = nmf_method,
                nrun = nmf_nrun_split, .options = "v")

    W1 <- basis(nmf1)
    W2 <- basis(nmf2)

    # Match components by cosine similarity
    match_res <- match_W_matrices(W1, W2)
    mean_match <- mean(match_res$matched, na.rm = TRUE)
    median_match <- median(match_res$matched, na.rm = TRUE)

    # Compute Pearson correlation on matched pairs
    ass <- match_res$assignment
    matched_pairs <- data.frame(
      a = seq_len(ncol(match_res$sim_mat)),
      b = ass[seq_len(ncol(match_res$sim_mat))]
    )
    good_pairs <- matched_pairs[matched_pairs$b <= ncol(W2), ]

    cors <- sapply(seq_len(nrow(good_pairs)), function(z) {
      a_idx <- good_pairs$a[z]
      b_idx <- good_pairs$b[z]
      common <- intersect(rownames(W1), rownames(W2))
      cor(W1[common, a_idx], W2[common, b_idx], method = "pearson")
    })

    mean_cor <- mean(cors, na.rm = TRUE)

    split_stats[[paste0("split", rep_i, "_rank", k)]] <- list(
      mean_cosine = mean_match,
      median_cosine = median_match,
      mean_cor = mean_cor
    )

  } # end rank loop

  if (rep_i %% 10 == 0) message("Completed split ", rep_i)
}

# Summarize split-half reproducibility
split_summary <- do.call(rbind, lapply(ranks_to_use, function(k) {
  keys <- grep(paste0("_rank", k, "$"), names(split_stats), value = TRUE)
  vals <- do.call(rbind, lapply(keys, function(x) unlist(split_stats[[x]])))
  data.frame(
    rank = k,
    mean_cosine_mean = mean(vals[, "mean_cosine"], na.rm = TRUE),
    mean_cosine_sd   = sd(vals[, "mean_cosine"], na.rm = TRUE),
    mean_cor_mean    = mean(vals[, "mean_cor"], na.rm = TRUE),
    mean_cor_sd      = sd(vals[, "mean_cor"], na.rm = TRUE)
  )
}))

write.csv(split_summary,
          file = file.path(out_dir, "split_half_summary.csv"),
          row.names = FALSE)

message("Split-half summary saved.")


# ------------------------------
# Factor Match Score (FMS) approximate
# ------------------------------
fms_summary <- list()
for (rk in names(all_rank_results)) {
  info <- all_rank_results[[rk]]
  run_Ws <- info$run_Ws
  run_names <- names(run_Ws)
  if (length(run_names) < 2) next
  # sample up to 200 random pairs (or all if fewer)
  pairs <- combn(run_names, 2)
  npairs <- ncol(pairs)
  samp_pairs <- if (npairs > 200) { sample(seq_len(npairs), 200) } else seq_len(npairs)
  pair_res <- sapply(samp_pairs, function(j) {
    a <- run_Ws[[ pairs[1, j] ]]
    b <- run_Ws[[ pairs[2, j] ]]
    mr <- match_W_matrices(a, b)
    mean(mr$matched, na.rm = TRUE)
  })
  fms_summary[[rk]] <- data.frame(rank = as.integer(sub("rank", "", rk)),
                                  fms_mean = mean(pair_res, na.rm = TRUE),
                                  fms_sd = sd(pair_res, na.rm = TRUE),
                                  n_pairs = length(pair_res))
}
fms_df <- bind_rows(fms_summary)
write.csv(fms_df, file = file.path(out_dir, "FMS_summary.csv"), row.names = FALSE)
message("FMS summary saved.")

# ------------------------------
# Save rank-level summaries
# ------------------------------
rank_overview <- lapply(names(rank_summaries), function(rk) {
  info <- rank_summaries[[rk]]
  data.frame(
    rank = as.integer(sub("rank", "", rk)),
    n_runs = length(info$runs_list),
    cophenetic = info$cophenetic,
    mean_silhouette = info$mean_silhouette,
    mean_cluster_stability = mean(info$cluster_stability, na.rm = TRUE),
    n_SE_universe = length(info$SE_universe)
  )
})
rank_overview_df <- bind_rows(rank_overview)
write.csv(rank_overview_df, file = file.path(out_dir, "rank_overview_summary.csv"), row.names = FALSE)
message("Rank overview saved: ", file.path(out_dir, "rank_overview_summary.csv"))

# Save bootstrap summaries
for (rk in names(boot_stats)) {
  write.csv(boot_stats[[rk]]$summary, file = file.path(out_dir, paste0(rk, "_bootstrap_summary_table.csv")), row.names = FALSE)
}


saveRDS(rank_summaries, file = file.path(out_dir, "rank_summaries_all.rds"))
saveRDS(boot_stats, file = file.path(out_dir, "bootstrap_details.rds"))
saveRDS(split_stats, file = file.path(out_dir, "split_half_details.rds"))

message("All results saved to: ", out_dir)

# ------------------------------
#plotting
# ------------------------------
# 1) Cophenetic + CI per rank
cop_df <- bind_rows(lapply(names(boot_stats), function(rk) {
  s <- boot_stats[[rk]]$summary
  data.frame(rank = as.integer(sub("rank", "", rk)),
             metric = s$metric,
             observed = s$observed,
             CI_low = s$CI_low,
             CI_high = s$CI_high,
             boot_mean = s$boot_mean)
}))
cop_plot_df <- cop_df %>% filter(metric == "cophenetic")
ggplot(cop_plot_df, aes(x = factor(rank), y = observed)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  theme_bw() + xlab("NMF rank (k)") + ylab("Cophenetic coefficient") +
  ggtitle("Cophenetic coefficient per rank (95% bootstrap CI)")
ggsave(file.path(out_dir, "cophenetic_per_rank.png"), width = 7, height = 5, dpi = 300)

# 2) Silhouette
sil_plot_df <- cop_df %>% filter(metric == "silhouette")
ggplot(sil_plot_df, aes(x = factor(rank), y = observed)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  theme_bw() + xlab("NMF rank (k)") + ylab("Mean silhouette width") +
  ggtitle("Mean silhouette per rank (95% bootstrap CI)")
ggsave(file.path(out_dir, "silhouette_per_rank.png"), width = 7, height = 5, dpi = 300)

# 3) FMS summary
if (nrow(fms_df) > 0) {
  ggplot(fms_df, aes(x = factor(rank), y = fms_mean)) +
    geom_point() + geom_errorbar(aes(ymin = fms_mean - fms_sd, ymax = fms_mean + fms_sd), width = 0.2) +
    theme_bw() + xlab("NMF rank (k)") + ylab("FMS (mean matched cosine)") +
    ggtitle("Factor-match (pairwise component reproducibility)")
  ggsave(file.path(out_dir, "FMS_per_rank.png"), width = 7, height = 5, dpi = 300)
}

message("Plots saved to ", out_dir)
message("NMF stability analysis complete. Check: ", out_dir)

