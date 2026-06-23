
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

cat("\n[1] Loading GBMx_all...\n")
obj <- readRDS(input_file)

DefaultAssay(obj) <- "peaks"

cat("Total cells:", ncol(obj), "\n")


library(lme4)
library(dplyr)
library(ggplot2)
library(performance)

meta <- obj@meta.data


meta$patient <- meta$sample
meta$cluster <- meta$seurat_clusters

signatures <- c(
  "state1_2fc_GA_score",
  "state2_2fc_GA_score",
  "state5_2fc_GA_score",
  "p16_top150_GA_score",
  "senmayo_GA_score"
)

library(ggplot2)

for (sig in signatures) {

  p <- ggplot(meta, aes(x = .data[[sig]])) +
  geom_density(fill = "grey70", color = "black", linewidth = 0.3) +
  theme_classic(base_size = 12) +
  labs(x = sig, y = "Density") +
  theme(
    plot.title = element_text(face = "bold"),
    axis.line = element_line()
  )

  ggsave(
    filename = paste0("Fig1_", sig, ".tiff"),
    plot = p,
    device = "tiff",
    width = 5,
    height = 4,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
}

run_lmm <- function(sig) {

  sub <- meta[!is.na(meta[[sig]]), ]

  m1 <- lmer(
    as.formula(paste(sig, "~ cluster + (1|patient)")),
    data = sub, REML = FALSE
  )

  m0 <- lmer(
    as.formula(paste(sig, "~ 1 + (1|patient)")),
    data = sub, REML = FALSE
  )

  lrt <- anova(m0, m1)
  r2_vals <- performance::r2(m1)

  data.frame(
    signature = sig,
    p_value = lrt$`Pr(>Chisq)`[2],
    R2_marginal = r2_vals$R2_marginal,
    R2_conditional = r2_vals$R2_conditional
  )
}
fig2_stats <- bind_rows(lapply(signatures, run_lmm))

write.csv(fig2_stats, "Fig2_mixed_model.csv", row.names = FALSE)


library(tidyr)


fig2_stats_long <- fig2_stats %>%
  pivot_longer(
    cols = c(R2_marginal, R2_conditional),
    names_to = "type",
    values_to = "R2"
  )

p <- ggplot(fig2_stats_long, aes(x = signature, y = R2, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = expression(R^2), fill = NULL) +
  theme(
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank()
  )

ggsave(
  "Fig2_R2_barplot.tiff",
  plot = p,
  device = "tiff",
  width = 5,
  height = 4,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# ==============================
# CONTINUOUS TRANSITION MODEL
# ==============================


meta_sub <- meta %>%
  dplyr::filter(
    !is.na(state1_2fc_GA_score),
    !is.na(state2_2fc_GA_score),
    !is.na(state5_2fc_GA_score)
  )

# FULL MODEL (mechanistic)
m_full <- lmer(
  state5_2fc_GA_score ~ 
    state2_2fc_GA_score + 
    state1_2fc_GA_score + 
    (1|patient),
  data = meta_sub,
  REML = FALSE
)

# REDUCED MODELS (for LRT)
m_no_state2 <- lmer(
  state5_2fc_GA_score ~ 
    state1_2fc_GA_score + 
    (1|patient),
  data = meta_sub,
  REML = FALSE
)

m_no_state1 <- lmer(
  state5_2fc_GA_score ~ 
    state2_2fc_GA_score + 
    (1|patient),
  data = meta_sub,
  REML = FALSE
)

# ==============================
# MODEL COMPARISON
# ==============================

lrt_state2 <- anova(m_no_state2, m_full)
lrt_state1 <- anova(m_no_state1, m_full)

# ==============================
# EFFECT SIZES
# ==============================

coef_table <- summary(m_full)$coefficients

beta_state2 <- coef_table["state2_2fc_GA_score", "Estimate"]
beta_state1 <- coef_table["state1_2fc_GA_score", "Estimate"]

# ==============================
# VARIANCE EXPLAINED
# ==============================

r2_vals <- performance::r2(m_full)



transition_stats <- data.frame(
  beta_state2 = beta_state2,
  beta_state1 = beta_state1,
  p_state2 = lrt_state2$`Pr(>Chisq)`[2],
  p_state1 = lrt_state1$`Pr(>Chisq)`[2],
  R2_marginal = r2_vals$R2_marginal,
  R2_conditional = r2_vals$R2_conditional
)

write.csv(transition_stats,
          "State_transition_model.csv",
          row.names = FALSE)

print(transition_stats)
###################################################

library(Seurat)
library(Signac)
library(dplyr)
library(lme4)
library(Matrix)
library(GenomicRanges)
library(rtracklayer)
library(future)
library(dplyr)
library(ggplot2)

setwd("~/Documents/Milan/Epigenomics/data/scatac")

plan("multicore", workers = 8)   
options(future.globals.maxSize = 50 * 1024^3)


obj <- readRDS("Cancer_scATACseq_data/GBMx_all.rds")

bulk_senescence_peaks_file <- "Cancer_scATACseq_data/atac_up.bed"

cat("\n[7] ATAC_up in PIS integration...\n")

bulk_peaks <- import(bulk_senescence_peaks_file)
sc_peaks <- StringToGRanges(rownames(obj))

hits <- findOverlaps(sc_peaks, bulk_peaks)
length(hits)

idx <- unique(queryHits(hits))


peak_names <- rownames(obj[["peaks"]])


if (length(idx) == 0) stop("No overlapping peaks found!")
if (max(idx) > length(peak_names)) stop("Index exceeds peak names!")


# Add module score
obj <- AddModuleScore(
  object = obj,
  features = list(ATAC_up = peak_names[idx]),
  assay = "peaks",
  name = "ATAC_up"
)


obj$ATAC_up_score <- obj$ATAC_up1

obj$ATAC_up1 <- NULL

meta <- obj@meta.data


meta$patient <- meta$sample

states <- c(
  "state2_2fc_GA_score",   
  "state5_2fc_GA_score"    
)

signatures <- c(
  "state1_2fc_GA_score",   
  "p16_top150_GA_score",
  "senmayo_GA_score",
  "IFN_I_response_GA_score",
  "MHC1_antigen_presentation_GO_0002428_GA_score",
  "SASP_GA_score",
  "ISG_up_AVIV_GA_score",
  "ISG_nature_GA_score",
  "cytoDNA_sensor_hsa_04623_GA_score",
  "Hallmark_IFNg_resp_GA_score",
  "Hallmark_IFNa_resp_GA_score",
  "Hallmark_Inflammatory_resp_GA_score",
  "ATAC_up_score",
  "E2F_TARGETS_HALLMARK_GA_score",
  "CELL_CYCLE_REGEV_GA_score"
)

# ==============================
# FUNCTION: LMM + LRT
# ==============================
run_model <- function(state, sig) {

  df <- meta %>%
    dplyr::select(all_of(c(state, sig, "patient"))) %>%
    na.omit()


  if (nrow(df) < 100) return(NULL)

 
  m0 <- lmer(
    as.formula(paste(state, "~ 1 + (1|patient)")),
    data = df,
    REML = FALSE
  )

  m1 <- lmer(
    as.formula(paste(state, "~", sig, "+ (1|patient)")),
    data = df,
    REML = FALSE
  )

  # LRT
  lrt <- anova(m0, m1)

 
  beta <- fixef(m1)[sig]

  data.frame(
    state = state,
    signature = sig,
    beta = beta,
    chisq = lrt$Chisq[2],
    p_value = lrt$`Pr(>Chisq)`[2],
    n_cells = nrow(df)
  )
}

# ==============================
# RUN ALL MODELS
# ==============================
cat("\n[1] Running LMM association models...\n")

results_list <- list()

for (st in states) {
  for (sig in signatures) {

    cat("Processing:", st, "vs", sig, "\n")

    res <- run_model(st, sig)

    if (!is.null(res)) {
      results_list[[paste(st, sig, sep = "_")]] <- res
    }
  }
}

results_table <- bind_rows(results_list)

# ==============================
# MULTIPLE TEST CORRECTION
# ==============================
results_table <- results_table %>%
  mutate(
    FDR = p.adjust(p_value, method = "fdr"),
    direction = ifelse(beta > 0, "positive", "negative")
  )


write.csv(results_table,
          "State_Signature_LMM_final.csv",
          row.names = FALSE)



# ==============================
# SUMMARY
# ==============================
summary_table <- results_table %>%
  group_by(state) %>%
  arrange(p_value) %>%
  slice(1:5)

write.csv(summary_table,
          "Top_associations_per_state.csv",
          row.names = FALSE)



# ==============================
# ATAC_up WILCOXON 
# ==============================

cat("\n[2] Running ATAC_up enrichment...\n")

meta <- obj@meta.data

run_wilcox <- function(state_var) {

  # define High / Low (top vs bottom 25%)
  q30 <- quantile(meta[[state_var]], 0.25, na.rm = TRUE)
  q70 <- quantile(meta[[state_var]], 0.75, na.rm = TRUE)

  group <- ifelse(
    meta[[state_var]] > q70, "High",
    ifelse(meta[[state_var]] < q30, "Low", NA)
  )

  df <- data.frame(
    group = group,
    ATAC_up_score = meta$ATAC_up_score
  ) %>% na.omit()

  df$group <- factor(df$group, levels = c("Low", "High"))

  # Wilcoxon test
  w <- wilcox.test(ATAC_up_score ~ group, data = df)

  # effect size
  med_high <- mean(df$ATAC_up_score[df$group == "High"])
  med_low  <- mean(df$ATAC_up_score[df$group == "Low"])

  data.frame(
    state = state_var,
    p_value = w$p.value,
    mean_high = med_high,
    mean_low = med_low,
    delta_mean = med_high - med_low,
    n_cells = nrow(df)
  )
}


wilcox_results <- bind_rows(
  lapply(states, run_wilcox)
)

write.csv(wilcox_results,
          "ATAC_up_enrichment_wilcox.csv",
          row.names = FALSE)

print(wilcox_results)



##################################################CRISPR CRE comparisons ###########


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

############################################################
# CRISPR TARGET REGIONS 
############################################################

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

############################################################
#ATAC_up MODULE SCORE 
############################################################

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

all_peaks <- rownames(obj[["peaks"]])

all_gr <- StringToGRanges(all_peaks, sep = c("-", "-"))

hits <- findOverlaps(crispr_gr, ctcf_gr)
obs <- length(unique(queryHits(hits)))

cat("\nObserved overlaps:", obs, "\n")

# ==============================
# PERMUTATION TEST 10000
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
# STATE5 ENRICHMENT vs CRISPR LOCI
############################################################

library(GenomicRanges)
library(Signac)

all_peaks <- rownames(obj[["peaks"]])

all_gr <- StringToGRanges(
  all_peaks,
  sep = c("-", "-")
)


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
###################################################state2 permutation test
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


hits <- findOverlaps(crispr_gr, state2_gr)
obs <- length(unique(queryHits(hits)))

cat("\nObserved state2 overlaps:", obs, "\n")

# ------------------------------
# PERMUTATION TEST (10,000)
# ------------------------------
set.seed(1)
n_perm <- 10000

perm_state2 <- replicate(n_perm, {


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





