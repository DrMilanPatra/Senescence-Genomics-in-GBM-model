library(Seurat)
library(nebula)

load(file="/home/user/Documents/Milan/Epigenomics/data/scRNA_seq_GBM4set_mice2set/GSE182109/gbm_GSE182109_filtered_scores.RData")

dir <- "~/nebula_gbm/"
if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

t1 <- Sys.time()
message("Running NEBULA for CDKN2A p16+ vs p16-")


cdkn2a_expr <- FetchData(se, vars = "CDKN2A")[,1]

se@meta.data$p16_class <- ifelse(
  cdkn2a_expr > 0.25, "p16_pos", "p16_neg"
)

se@meta.data$p16_class <- factor(
  se@meta.data$p16_class,
  levels = c("p16_neg", "p16_pos")
)

se <- CellCycleScoring(
  se,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = FALSE
)


subset_cells <- !is.na(se@meta.data$p16_class)
se_sub <- subset(se, cells = colnames(se)[subset_cells])


data_neb <- scToNeb(
  obj   = se_sub,
  assay = "RNA",
  id    = "Patient",
  pred  = "p16_class"
)

table(data_neb$pred)
str(data_neb$pred)



pred_df <- data_neb$pred

pred_df$p16_class <- factor(pred_df$p16_class,
                             levels=c("p16_neg","p16_pos"))

pred_df$S.Score   <- se_sub$S.Score
pred_df$G2M.Score <- se_sub$G2M.Score

pred_matrix <- model.matrix(~ p16_class + S.Score + G2M.Score,
                             data=pred_df)



dim(pred_matrix)



count_mat <- data_neb$count
id_vec <- as.factor(data_neb$id)

offset_vec <- log1p(colSums(count_mat))
names(offset_vec) <- NULL



##########################################
#Run NEBULA
##########################################
neb_res <- nebula(
  count  = count_mat,
  id     = id_vec,
  pred   = pred_matrix,
  offset = offset_vec,
  verbose = TRUE
)



neb_summary <- neb_res$summary

neb_summary$padj_p16 <- p.adjust(
  neb_summary$p_p16_classp16_pos,
  method = "BH"
)

write.table(
  neb_summary,
  file = paste0(dir, "NEBULA_gbm_CDKN2A_summary.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



sig_genes <- subset(
  neb_summary,
  padj_p16 < 0.05 &
  abs(logFC_p16_classp16_pos) > 0.25
)

sig_genes <- sig_genes[
  order(-sig_genes$logFC_p16_classp16_pos),
]

write.table(
  sig_genes,
  file = paste0(dir, "NEBULA_gbm_CDKN2A_sig.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


message("Finished in ", round(Sys.time() - t1, 2), " seconds.")