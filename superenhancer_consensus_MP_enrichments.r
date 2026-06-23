# ==========================================================
# META-PROGRAM → GENE ANNOTATION + GSEA + OVERLAP

# ==========================================================

library(dplyr)
library(tidyr)
library(tibble)
library(matrixStats)

base_dir <- "/home/user/Documents/Milan/Epigenomics"


hichip <- read.delim(
  file.path(base_dir, "SE_hichip_annotation_coordinates.txt"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)



gmt_path <- file.path(base_dir, "curated_list_gbm.gmt")

read_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file)
  
  gmt_list <- lapply(lines, function(x) {
    parts <- strsplit(x, "\t")[[1]]
    set_name <- parts[1]
    genes <- parts[-c(1,2)]  # remove name + description
    return(genes)
  })
  
  names(gmt_list) <- sapply(lines, function(x) strsplit(x, "\t")[[1]][1])
  return(gmt_list)
}

gmt_list <- read_gmt(gmt_path)

cat("Loaded gene sets:", length(gmt_list), "\n")


# ==========================================================
#MATCH SEs with meta-programs
# ==========================================================

common_SEs <- intersect(rownames(rpk_log_mat), unique(unlist(meta_programs)))
rpk_log_mat <- rpk_log_mat[common_SEs, ]

cat("Matched SEs:", length(common_SEs), "\n")

# ==========================================================
# CORE SE SELECTION (~top 40%)
# ==========================================================

meta_programs_core <- list()

for (mp in names(meta_programs)) {
  
  SEs <- intersect(meta_programs[[mp]], rownames(rpk_log_mat))
  
  if (length(SEs) < 10) next
  
  se_signal <- rowMeans(rpk_log_mat[SEs, , drop = FALSE])
  
  thresh <- quantile(se_signal, 0.6, na.rm = TRUE)
  core_SEs <- SEs[se_signal >= thresh]
  
  if (length(core_SEs) < 20) {
    core_SEs <- SEs[order(se_signal, decreasing = TRUE)][1:min(20, length(SEs))]
  }
  
  meta_programs_core[[mp]] <- core_SEs
}

cat("Core SEs:", length(meta_programs_core), "\n")

# ==========================================================
# SE → MP mapping
# ==========================================================

SE_to_MP <- data.frame(
  SE_id = unlist(meta_programs_core),
  meta_program = rep(names(meta_programs_core),
                     lengths(meta_programs_core))
)


library(dplyr)
SE_coords <- hichip %>%
  select(SE_id, chr_se, start_se, end_se) %>%
  distinct()

SE_MP_coords <- SE_to_MP %>%
  inner_join(SE_coords, by = "SE_id")

write.table(
  SE_MP_coords,
  file.path(base_dir, "MP_SE_coordinates.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("SE coordinates saved\n")



SE_gene_MP <- hichip %>%
  inner_join(SE_to_MP, by = "SE_id")



SE_signal_df <- rpk_log_mat %>%
  as.data.frame() %>%
  rownames_to_column("SE_id") %>%
  mutate(SE_signal = rowMeans(across(-SE_id)))

SE_gene_MP <- SE_gene_MP %>%
  inner_join(SE_signal_df[, c("SE_id", "SE_signal")], by = "SE_id")

# ==========================================================
# GENE PENALTY
# ==========================================================

gene_freq <- table(SE_gene_MP$gene)

SE_gene_MP <- SE_gene_MP %>%
  mutate(gene_penalty = log1p(gene_freq[gene]))

# ==========================================================
# BUILD GENE SCORES
# ==========================================================

SE_gene_weighted <- SE_gene_MP %>%
  group_by(meta_program, gene) %>%
  summarise(
    total_loop = sum(loop_score),
    max_SE_signal = max(SE_signal),
    penalty = mean(gene_penalty),
    .groups = "drop"
  )

# ==========================================================
# FINAL SCORE
# ==========================================================

SE_gene_ranked <- SE_gene_weighted %>%
  group_by(meta_program) %>%
  mutate(
    loop_z = scale(total_loop)[,1],
    se_z   = scale(max_SE_signal)[,1],
    pen_z  = scale(penalty)[,1],
    
    final_score = loop_z + se_z - pen_z,
    final_score = final_score + runif(n(), 0, 1e-6),
    rank = rank(-final_score, ties.method = "first")
  ) %>%
  arrange(meta_program, desc(final_score)) %>%
  ungroup()

cat("Gene scoring done\n")

# ==========================================================
# TOP GENES
# ==========================================================

top_genes_per_MP <- SE_gene_ranked %>%
  group_by(meta_program) %>%
  slice_max(order_by = final_score, n = 300) %>%
  ungroup()

write.table(
  SE_gene_ranked,
  file.path(base_dir, "MP_gene_ranked_FULL.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  top_genes_per_MP,
  file.path(base_dir, "MP_gene_top100.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ==========================================================
# GSEA INPUT 
# ==========================================================

ranked_lists <- split(SE_gene_ranked, SE_gene_ranked$meta_program)

ranked_lists <- lapply(ranked_lists, function(df) {
  v <- df$final_score
  names(v) <- df$gene
  sort(v, decreasing = TRUE)
})

saveRDS(
  ranked_lists,
  file.path(base_dir, "MP_ranked_lists_for_GSEA.rds")
)

cat("GSEA lists ready\n")

# ==========================================================
# CURATED OVERLAP (HYPERGEOMETRIC) 
# ==========================================================

all_genes <- unique(SE_gene_ranked$gene)

hyper_results <- list()

for(mp in names(ranked_lists)) {

  genes_mp <- names(ranked_lists[[mp]])[1:300]

  for(gs in names(gmt_list)) {

    geneset <- intersect(gmt_list[[gs]], all_genes)

    k <- length(intersect(genes_mp, geneset))
    M <- length(all_genes)
    n <- length(geneset)
    N <- length(genes_mp)

    if(n < 10) next

    pval <- phyper(k-1, n, M-n, N, lower.tail = FALSE)

    hyper_results[[paste(mp, gs, sep="_")]] <- data.frame(
      meta_program = mp,
      gene_set = gs,
      overlap = k,
      geneset_size = n,
      p_value = pval
    )
  }
}

hyper_df <- bind_rows(hyper_results)

# adjust p-values
hyper_df$FDR <- p.adjust(hyper_df$p_value, method = "BH")

write.table(
  hyper_df,
  file.path(base_dir, "MP_GMT_hypergeometric.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)


#######################GSEA RUN#####################


library(dplyr)
library(readr)
library(reticulate)

#use_condaenv("seurat411_env", required = TRUE)
gp <- import("gseapy")



output_dir <- "/home/user/Documents/Milan/Epigenomics/GSEA_results1/"
#output_dir <- "/home/user/Documents/Milan/Epigenomics/GSEA_results_curated1/"
dir.create(output_dir, showWarnings = FALSE)
gmt_path <- "/home/user/Documents/Milan/Epigenomics/h.all.v2025.1.Hs.symbols.gmt"

#gmt_path <- "/home/user/Documents/Milan/Epigenomics/curated_list_gbm.gmt"
all_results <- list()


for(mp in names(ranked_lists)) {

  cat("\nRunning GSEA for:", mp, "\n")

  vec <- ranked_lists[[mp]]
df_r <- data.frame(
  gene = names(vec),
  score = as.numeric(vec),
  stringsAsFactors = FALSE
)


df_r <- df_r[!is.na(df_r$score), ]
df_r <- df_r[!duplicated(df_r$gene), ]


df_r <- df_r[order(df_r$score, decreasing = TRUE), ]

# convert to pandas with correct structure
df_pd <- reticulate::r_to_py(df_r)

  mp_outdir <- file.path(output_dir, mp)
  dir.create(mp_outdir, showWarnings = FALSE)

  gp$prerank(
    rnk = df_pd,
    gene_sets = gmt_path,
    min_size = as.integer(10),   
    max_size = as.integer(1000),
    permutation_num = as.integer(2000),  
    outdir = mp_outdir,
    format = "png",
    seed = as.integer(42),
    verbose = TRUE
  )
}

##########################################
##extraction of gsea results
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

#base_dir <- "/home/user/Documents/Milan/Epigenomics/GSEA_results1"

base_dir <- "/home/user/Documents/Milan/Epigenomics/GSEA_results_curated1"

mp_dirs <- list.dirs(base_dir, recursive = FALSE)

all_gsea <- list()

for (mp_dir in mp_dirs) {
  
  mp_name <- basename(mp_dir)
  files <- list.files(mp_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(files) == 0) next
  
  df <- read_csv(files[1])
  df$meta_program <- mp_name
  
  all_gsea[[mp_name]] <- df
}

gsea_df <- bind_rows(all_gsea)

cat("Loaded GSEA results:", nrow(gsea_df), "\n")

library(dplyr)
library(stringr)

gsea_df <- gsea_df %>%
  mutate(
    Lead_genes = str_split(Lead_genes, ";")
  )
  
  mp_genes <- readRDS("/home/user/Documents/Milan/Epigenomics/MP_ranked_lists_for_GSEA.rds")

top_genes <- lapply(mp_genes, function(x) {
  names(x)[1:min(300, length(x))]
})

results <- list()

for (mp in unique(gsea_df$meta_program)) {
  
  mp_set <- top_genes[[mp]]
  if (is.null(mp_set)) next
  
  df_mp <- gsea_df %>% filter(meta_program == mp)
  
  res <- lapply(seq_len(nrow(df_mp)), function(i) {
    
    pathway_genes <- df_mp$Lead_genes[[i]]
    
    inter <- length(intersect(mp_set, pathway_genes))
    
    overlap_frac <- inter / length(mp_set)
    pathway_cov   <- inter / length(pathway_genes)
    
    nes <- df_mp$NES[i]
    
    # FINAL ROBUST SCORE **VVIP
    final_score <- (0.6 * nes) + (0.4 * overlap_frac)
    
    data.frame(
      meta_program = mp,
      pathway = df_mp$Term[i],
      NES = nes,
      overlap = inter,
      overlap_frac = overlap_frac,
      pathway_coverage = pathway_cov,
      FDR = df_mp$`FDR q-val`[i],
      score = final_score
    )
  })
  
  results[[mp]] <- bind_rows(res)
}

final_gsea_matrix <- bind_rows(results)

library(tidyr)

heatmap_mat <- final_gsea_matrix %>%
  select(meta_program, pathway, score) %>%
  pivot_wider(names_from = pathway, values_from = score, values_fill = 0) %>%
  as.data.frame()

rownames(heatmap_mat) <- heatmap_mat$meta_program
heatmap_mat <- as.matrix(heatmap_mat[,-1])

library(pheatmap)
library(stringr)

heatmap_mat <- heatmap_mat[order(as.numeric(str_extract(rownames(heatmap_mat), "\\d+"))), ]




png("MP_GSEA_consensus_heatmap.png", width = 4500, height = 3500, res = 300)
#png("MP_GSEA_consensus_heatmap_curated1.png", width = 4500, height = 3500, res = 300)

pheatmap(
  heatmap_mat,
  color = colorRampPalette(c("#FFFFFF", "#FF8C00", "#4B0082"))(100),
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  main = "Meta-program × Pathway Consensus Enrichment (NES + Overlap)"
)

dev.off()

