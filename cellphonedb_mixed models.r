

library(data.table)
library(Seurat)
library(SingleCellExperiment)
library(S4Vectors)
library(data.table)
library(ktplots)


base_path <- "/home/user/CellPhoneDB_cellstate_outputs/results"

patients <- list.dirs(base_path, full.names = FALSE, recursive = FALSE)

pval_thresh <- 0.05


for (patient in patients) {

  cat("\n=============================\n")
  cat("Processing:", patient, "\n")
  cat("=============================\n")

  #########################################
  # Define patient paths (INSIDE loop)
  #########################################
  out_path <- file.path(base_path, patient)

  means_file <- file.path(out_path, "statistical_analysis_means_cellstate.txt")
  pvals_file <- file.path(out_path, "statistical_analysis_pvalues_cellstate.txt")

  if (!file.exists(means_file) | !file.exists(pvals_file)) {
    cat("Skipping (missing files)\n")
    next
  }


  means_stat <- fread(means_file, data.table = FALSE)
  pvals_stat <- fread(pvals_file, data.table = FALSE)


  state1_cols <- grep("State1_Glioma", colnames(pvals_stat), value = TRUE)
  state2_cols <- grep("State2_Glioma", colnames(pvals_stat), value = TRUE)
  state5_cols <- grep("State5_Glioma", colnames(pvals_stat), value = TRUE)

  sig_cols <- c(state1_cols, state2_cols, state5_cols)

  sig_idx <- apply(
    pvals_stat[, sig_cols],
    1,
    function(x) any(x < pval_thresh, na.rm = TRUE)
  )

  cat("Total interactions:", nrow(pvals_stat), "\n")
  cat("Significant interactions:", sum(sig_idx), "\n")

  means_sig <- means_stat[sig_idx, ]
  pvals_sig <- pvals_stat[sig_idx, ]

  saveRDS(means_sig,
          file = file.path(out_path, paste0("means_SIG_State1_2_5_", patient, ".rds")))

  saveRDS(pvals_sig,
          file = file.path(out_path, paste0("pvals_SIG_State1_2_5_", patient, ".rds")))

  fwrite(means_sig,
         file = file.path(out_path, paste0("means_SIG_State1_2_5_", patient, ".txt")),
         sep = "\t")

  fwrite(pvals_sig,
         file = file.path(out_path, paste0("pvals_SIG_State1_2_5_", patient, ".txt")),
         sep = "\t")

  cat("Saved results for:", patient, "\n")
}

cat("\nALL PATIENTS COMPLETED\n")



####################################################################################################################################


library(data.table)
library(lme4)
library(dplyr)
library(reshape2)

base_path <- "/home/user/CellPhoneDB_cellstate_outputs/results"
patients <- list.dirs(base_path, full.names = FALSE, recursive = FALSE)

all_df <- list()

for (patient in patients) {

  cat("\nProcessing:", patient, "\n")

  out_path <- file.path(base_path, patient)

  means_file <- file.path(out_path,
                          paste0("means_SIG_State1_2_5_", patient, ".rds"))

  if (!file.exists(means_file)) next

  means_sig <- readRDS(means_file)

  if (is.null(means_sig) || nrow(means_sig) == 0) {
    cat("No signal for:", patient, "\n")
    next
  }


  df <- means_sig

  state_cols <- grep("State1_Glioma|State2_Glioma|State5_Glioma",
                     colnames(df), value = TRUE)

  df_long <- melt(df,
                  id.vars = c("id_cp_interaction", "interacting_pair"),
                  measure.vars = state_cols,
                  variable.name = "cell_pair",
                  value.name = "interaction_strength")

  setDT(df_long)

  df_long[, c("cell_type", "state") := tstrsplit(cell_pair, "\\|")]



  df_long[, state := gsub("_Glioma", "", state)]



  df_long[, patient := patient]



  df_long <- df_long[interaction_strength > 0]

  all_df[[patient]] <- df_long
}


df_all <- rbindlist(all_df, fill = TRUE)


keep_states <- c("State1", "State2", "State5")
keep_cells  <- c("TCells", "Myeloid", "BCells")

df_all <- df_all[
  state %in% keep_states &
  cell_type %in% keep_cells
]



df_all[, state := factor(state)]
df_all[, patient := factor(patient)]
df_all[, cell_type := factor(cell_type)]
df_all[, interacting_pair := factor(interacting_pair)]

#########################################
#Linear Mixed Model
#########################################

library(emmeans)
library(lme4)
library(data.table)

df_all[, interaction_strength := log1p(interaction_strength)]



m2 <- lmer(
  interaction_strength ~ state * cell_type +
    (1 | patient) +
    (1 | interacting_pair),
  data = df_all
)

emm <- emmeans(m2, ~ state | cell_type)


pairs_res <- pairs(emm)

emm_df <- as.data.table(emm)
pairs_df <- as.data.table(pairs_res)

#summary(emm)

fwrite(emm_df,
       "CPDB_LMM_emmeans_state_by_LR_celltype.csv")
       
       fwrite(pairs_df,
       "CPDB_LMM_state_comparisons_by_LR_celltype.csv")
       
       saveRDS(m2, "CPDB_LMM_model.rds")

write.csv(df_all, "CPDB_LMM_long_format.csv", row.names = FALSE)

cat("\nLMM DATASET CREATED\n")


#######plotiing

library(ggplot2)

emm_df$state <- factor(
  emm_df$state,
  levels = c("State1","State2","State5")
)

tiff("cellphonedb_global.tiff",
     units = "in", width = 5, height = 5,
     res = 1200, compression = "lzw")

p <- ggplot(
  emm_df,
  aes(
    x = state,
    y = emmean,
    group = 1
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    width = 0.08
  ) +
 facet_wrap(
  ~cell_type,
  nrow = 1,
  scales = "fixed"
) +
  ylab("Adjusted interaction strength") +
  xlab("") +
  theme_classic(base_size = 14) +

  theme(
axis.text = element_text(color = "black"),
axis.line = element_line(linewidth = 0.8),
axis.ticks = element_line(linewidth = 0.8),
	axis.ticks.length = unit(0.45, "cm"),
    legend.position = "top",
    legend.text = element_text(size = 11),
  axis.title.y = element_text(size = 12),
  axis.text.y  = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 6),
  axis.text.x  = element_text(size = 9, color = "black")
)
  
p

dev.off()

########################################################################


library(Seurat)
library(SingleCellExperiment)
library(S4Vectors)
library(data.table)
library(ktplots)



path <- "/home/user/CellPhoneDB_cellstate_subtype_outputs/results/ndGBM-01"

means <- fread(file.path(path, "statistical_analysis_means_cellstate.txt"))
pvals <- fread(file.path(path, "statistical_analysis_pvalues_cellstate.txt"))

decon <- fread(file.path(path, "statistical_analysis_deconvoluted_cellstate.txt"))

means <- as.data.frame(means)
pvals <- as.data.frame(pvals)
decon <- as.data.frame(decon)

load("/home/user/Documents/Milan/Epigenomics/data/scRNA_seq_GBM4set_mice2set/GSE182109/gbm_GSE182109_cellstate_subtype_cellphonedb.RData")

sce_sce <- SingleCellExperiment(
  assays = list(counts = sce@assays$RNA@counts)
)

colData(sce_sce) <- S4Vectors::DataFrame(sce@meta.data)

colData(sce_sce)$celltype <- sce@meta.data$celltype_state


lr <- means[, c("interacting_pair", "gene_a", "gene_b")]



state5_up <- read.delim(
    "/home/user/Documents/Milan/scripts/state5_up.txt",
    header = TRUE
)

genes_up <- unique(state5_up$gene)


lr_keep <- means$gene_a %in% genes_up

means_filt <- means[lr_keep, ]
pvals_filt <- pvals[lr_keep, ]


pattern_lr <- paste(
  c(
    # Cytokine / SASP
    "TGFB|IL1|CSF1|HGF|LGALS|KITLG",
    
    # Chemokine
    "CCL|CXCL|CX3CL|CCR|CXCR",
    
    # Immune checkpoint
    "PVR|TIGIT|CD96|CD226|HLA|VSIR|LILR|NKG2|MICA|CD58|CD99|PTPRC|BST2",
    
    # Phagocytosis
    "APOE|GAS6|AXL|MERTK|TYRO3|PLAU|PLAUR|TREM2|PSAP|LRP5"
  ),
  collapse = "|"
)

lr_keep_pathway <- grepl(pattern_lr, means_filt$interacting_pair)

means_filt2 <- means_filt[lr_keep_pathway, ]
pvals_filt2 <- pvals_filt[lr_keep_pathway, ]


library(ggplot2)
library(ggrepel)



tiff("ktplot_cd8_state5_upregulated_pathwayFiltered.tiff",
     units = "in", width = 18, height = 12, res = 1200, compression = "lzw")

p<- plot_cpdb2(
    scdata = sce_sce,
    cell_type1 = "CD8_TCells",
    cell_type2 = ".",
    celltype_key = "celltype",
    means = means_filt2,
    pvals = pvals_filt2,
    deconvoluted = decon,
    desiredInteractions = list(
        c("CD8_TCells", "State2_Glioma"),
        c("CD8_TCells", "State5_Glioma")
    )
)

p$layers[[3]] <- NULL

p <- p +
  ggrepel::geom_text_repel(
    aes(x = x, y = y, label = label), force=2,
    size = 6
  ) +
  scale_size_continuous(range = c(3, 10)) +
  scale_colour_manual(
    values = c(
      "CD8_TCells" = "#1B9E77",
      "State2_Glioma" = "#377EB8",
      "State5_Glioma" = "#E08B00"
    )
  ) +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),

    legend.position = c(1.08, 0.5),
    legend.justification = c(0, 0.5)
  )
  
  print(p)
  
dev.off()


tiff("ktplot_mdsc_state5_upregulated_pathwayFiltered1.tiff",
     units = "in", width = 18, height = 12, res = 1200, compression = "lzw")
p <- plot_cpdb2(
    scdata = sce_sce,
    cell_type1 = "MDSC",
    cell_type2 = ".",
    celltype_key = "celltype",
    means = means_filt2,
    pvals = pvals_filt2,
    deconvoluted = decon,
    desiredInteractions = list(
        c("MDSC", "State2_Glioma"),
        c("MDSC", "State5_Glioma")
    )
)

p$layers[[3]] <- NULL

p <- p +
  ggrepel::geom_text_repel(
    aes(x = x, y = y, label = label), force=2,
    size = 6
  ) +
  scale_size_continuous(range = c(3, 10)) +
  scale_colour_manual(
    values = c(
      "MDSC" = "#1B9E77",
      "State2_Glioma" = "#377EB8",
      "State5_Glioma" = "#E08B00"
    )
  ) +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),

    legend.position = c(1.08, 0.5),
    legend.justification = c(0, 0.5)
  )

print(p)
dev.off()


############################################################################################################################




library(data.table)
library(lme4)
library(dplyr)
library(reshape2)
library(limma)
library(statmod)

#base_path <- "/home/user/CellPhoneDB_cellstate_outputs/results"
base_path <- "/home/user/CellPhoneDB_cellstate_subtype_outputs/results"
patients <- list.dirs(base_path, full.names = FALSE, recursive = FALSE)

all_df <- list()


for (patient in patients) {

  cat("\nProcessing:", patient, "\n")

  out_path <- file.path(base_path, patient)

  means_file <- file.path(out_path,
                          paste0("means_SIG_State1_2_5_", patient, ".rds"))

  if (!file.exists(means_file)) next

  means_sig <- readRDS(means_file)

  if (is.null(means_sig) || nrow(means_sig) == 0) {
    cat("No signal for:", patient, "\n")
    next
  }



  df <- means_sig

  state_cols <- grep("State1_Glioma|State2_Glioma|State5_Glioma",
                     colnames(df), value = TRUE)

  df_long <- melt(df,
                  id.vars = c("id_cp_interaction", "interacting_pair"),
                  measure.vars = state_cols,
                  variable.name = "cell_pair",
                  value.name = "interaction_strength")

  setDT(df_long)

  df_long[, c("cell_type", "state") := tstrsplit(cell_pair, "\\|")]



  df_long[, state := gsub("_Glioma", "", state)]



  df_long[, patient := patient]



  df_long <- df_long[interaction_strength > 0]

  all_df[[patient]] <- df_long
}


df_all <- rbindlist(all_df, fill = TRUE)

#saveRDS(df_all, file = "df_all_CellPhoneDB_Statewise.rds")

###Check for particular L_R interactions
df_all <- readRDS("df_all_CellPhoneDB_Statewise.rds")

pattern_lr <- paste(
  c(
    # Cytokine / SASP
    "TGFB|IL1|CSF1|HGF|LGALS|KITLG",
    
    # Chemokine
    "CCL|CXCL|CX3CL|CCR|CXCR",
    
    # Immune checkpoint
    "PVR|TIGIT|CD96|CD226|HLA|VSIR|LILR|NKG2|MICA|CD58|CD99|PTPRC|BST2",
    
    # Phagocytosis
    "APOE|GAS6|AXL|MERTK|TYRO3|PLAU|PLAUR|TREM2|PSAP|LRP5"
  ),
  collapse = "|"
)



df_lr <- df_all[grepl(pattern_lr, interacting_pair)]

#df_lr <- copy(df_all)

library(lmerTest)

df_cd8 <- df_lr[
  cell_type == "CD8_TCells" &
  state %in% c("State2", "State5")
]

df_cd8$state <- relevel(factor(df_cd8$state), ref = "State2")




###Model for subset of immune cell interaction in state2 (ref) vs state5
model_cd8 <- lmer(
  log1p(interaction_strength) ~ state +
    (1|patient) +
    (1|interacting_pair),
  data=df_cd8
)

summary(model_cd8)


df_NKCells <- df_lr[
  cell_type == "NKCells" &
  state %in% c("State2", "State5")
]

df_NKCells$state <- relevel(factor(df_NKCells$state), ref = "State2")

model_NKCells <- lmer(
  log1p(interaction_strength) ~ state +
    (1|patient) +
    (1|interacting_pair),
  data=df_NKCells
)

summary(model_NKCells)


df_MDSC <- df_lr[
  cell_type == "MDSC" &
  state %in% c("State2", "State5")
]

df_MDSC$state <- relevel(factor(df_MDSC$state), ref = "State2")

model_MDSC <- lmer(
  log1p(interaction_strength) ~ state +
    (1|patient) +
    (1|interacting_pair),
  data=df_MDSC
)

summary(model_MDSC)


df_Tregs <- df_lr[
  cell_type == "Tregs" &
  state %in% c("State2", "State5")
]

df_Tregs$state <- relevel(factor(df_Tregs$state), ref = "State2")


model_Tregs <- lmer(
  log1p(interaction_strength) ~ state +
    (1|patient) +
    (1|interacting_pair),
  data=df_Tregs
)

summary(model_Tregs)


df_s_mac_1 <- df_lr[
  cell_type == "s_mac_1" &
  state %in% c("State2", "State5")
]

df_s_mac_1$state <- relevel(factor(df_s_mac_1$state), ref = "State2")

model_s_mac_1 <- lmer(
  log1p(interaction_strength) ~ state +
    (1|patient) +
    (1|interacting_pair),
  data=df_s_mac_1
)

summary(model_s_mac_1)


df_s_mac_2 <- df_lr[
  cell_type == "s_mac_2" &
  state %in% c("State2", "State5")
]

df_s_mac_2$state <- relevel(factor(df_s_mac_2$state), ref = "State2")



model_s_mac_2 <- lmer(
  log1p(interaction_strength) ~ state +
    (1|patient) +
    (1|interacting_pair),
  data=df_s_mac_2
)

summary(model_s_mac_2)


###################################################################
###mixed model for each interacting L-R pair

df_cd8 <- df_lr[
  cell_type == "CD8_TCells" &
  state %in% c("State2", "State5")
]

df_cd8$state <- relevel(factor(df_cd8$state), ref = "State2")


model_cd8 <- lmer(
  log1p(interaction_strength) ~  state + interacting_pair  + (1|patient),
  data=df_cd8
)

summary(model_cd8)


df_NKCells <- df_lr[
  cell_type == "NKCells" &
  state %in% c("State2", "State5")
]

df_NKCells$state <- relevel(factor(df_NKCells$state), ref = "State2")



model_NKCells <- lmer(
  log1p(interaction_strength) ~ state + interacting_pair + (1 | patient),
  data = df_NKCells
)

summary(model_NKCells)



####################################

df_MDSC <- df_lr[
  cell_type == "MDSC" &
  state %in% c("State2", "State5")
]

df_MDSC$state <- relevel(factor(df_MDSC$state), ref = "State2")



model_MDSC <- lmer(
  log1p(interaction_strength) ~ state + interacting_pair + (1 | patient),
  data = df_MDSC
)

summary(model_MDSC)


####################################

df_Tregs <- df_lr[
  cell_type == "Tregs" &
  state %in% c("State2", "State5")
]

df_Tregs$state <- relevel(factor(df_Tregs$state), ref = "State2")



model_Tregs <- lmer(
  log1p(interaction_strength) ~ state + interacting_pair + (1 | patient),
  data = df_Tregs
)

summary(model_Tregs)


####################################

df_s_mac_1 <- df_lr[
  cell_type == "s_mac_1" &
  state %in% c("State2", "State5")
]

df_s_mac_1$state <- relevel(factor(df_s_mac_1$state), ref = "State2")



model_s_mac_1 <- lmer(
  log1p(interaction_strength) ~ state + interacting_pair + (1 | patient),
  data = df_s_mac_1
)

summary(model_s_mac_1)

####################################

df_s_mac_2 <- df_lr[
  cell_type == "s_mac_2" &
  state %in% c("State2", "State5")
]

df_s_mac_2$state <- relevel(factor(df_s_mac_2$state), ref = "State2")


model_s_mac_2 <- lmer(
  log1p(interaction_strength) ~ state + interacting_pair + (1 | patient),
  data = df_s_mac_2
)

summary(model_s_mac_2)


####################################################Plotting
library(dplyr)
library(broom.mixed)
library(ggplot2)
library(grid)

# Models
models <- list(
  s_mac_1 = model_s_mac_1,
  s_mac_2 = model_s_mac_2,
  Treg    = model_Tregs,
  MDSC    = model_MDSC,
  NK      = model_NKCells,
  CD8     = model_cd8
)

for (cell in names(models)) {

  # -----------------------------
  # 1. Extract ONLY LR × state interaction terms
  # -----------------------------
res_full <- broom.mixed::tidy(models[[cell]], effects = "fixed") %>%
  filter(
    term != "(Intercept)",
    term != "stateState5"
  ) %>%
  mutate(
    term = gsub("^interacting_pair", "", term),
    CI_low  = estimate - 1.96 * std.error,
    CI_high = estimate + 1.96 * std.error,
    FDR = p.adjust(p.value, method = "BH")
  )

  # Skip if no interaction terms found
  if (nrow(res_full) == 0) next


  res_adj <- res_full %>%
  filter(FDR < 0.05)
  
  # Save full table
  write.csv(
    res_adj,
    paste0(cell, "_LR_state5_vs_state2_all_results.csv"),
    row.names = FALSE
  )

  # -----------------------------
  # 3. Select significant LR rewiring events
  # -----------------------------

res <- res_full %>%
  filter(FDR < 0.05) %>%
  arrange(desc(abs(estimate))) %>%
  slice_head(n = 20)
  
  
  if (nrow(res) == 0) next

  res$term <- factor(res$term, levels = rev(res$term))

  # -----------------------------
  # 4. Plot (TRUE LR rewiring forest plot)
  # -----------------------------
  p <- ggplot(res, aes(x = estimate, y = term)) +

    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +

    geom_errorbarh(
      aes(xmin = CI_low, xmax = CI_high),
      height = 0.25,
      linewidth = 0.7
    ) +

    geom_point(size = 2.5) +

    labs(
      title = paste0(cell, ": State5 vs State2 LR rewiring"),
      x = expression(beta~"(State interaction effect)"),
      y = NULL
    ) +

    theme_classic(base_size = 22) +

    theme(
      plot.title = element_text(
        size = 10,
        face = "bold",
        hjust = 0.5
      ),
      axis.title.x = element_text(size = 9),
      axis.text.x = element_text(size = 18, color = "black"),
      axis.text.y = element_text(size = 18, color = "black"),
      axis.line = element_line(size = 1.05),
      axis.ticks = element_line(size = 1.05),
      axis.ticks.length = unit(0.5, "cm")
    )

  # -----------------------------
  # 5. Save figure
  # -----------------------------
  tiff(
    paste0(cell, "_LR_rewiring_top20.tiff"),
    units = "in",
    width = 6.5,
    height = 6,
    res = 1200,
    compression = "lzw"
  )

  print(p)
  dev.off()
}