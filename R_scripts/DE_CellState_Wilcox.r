############################################################
# FindMarkers for CellState (one cell state -vs-rest) 
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)

load("C:/Users/mlnpa/Documents/glioblastoma_scrnaseq/gbm_GSE182109_filtered_scores_heterogeneity.RData")


if(!"CellState" %in% colnames(se@meta.data)){
  stop("CellState metadata missing in Seurat object.")
}

print(table(se$CellState))


Idents(se) <- se$CellState

# -------------------------------
#Differential expression loop
# -------------------------------
states <- levels(se$CellState)

for(st in states){
  
  message("Running DE for: ", st)
  
  # -------------------------------
  # Run FindMarkers (State vs rest)
  # -------------------------------
  res <- FindMarkers(
    se,
    ident.1 = st,
    ident.2 = NULL,
    assay = "RNA",
    slot = "data",
    test.use = "wilcox",
    logfc.threshold = 0.25
  )
  

  res$gene <- rownames(res)
  res$state <- st
  
  # -------------------------------
  # Multiple hypothesis correction
  # -------------------------------
  res$padj <- p.adjust(res$p_val, method = "BH")
  

  write.table(
    res,
    file = paste0("CellState_DEG_", st, "_FULL.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # -------------------------------
  # Filter significant genes (padj < 0.1)
  # -------------------------------
  res_sig <- res %>%
    filter(padj < 0.1)

  write.table(
    res_sig,
    file = paste0("CellState_DEG_", st, "_padj0.1.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}



################################################################# 
library(Seurat)
library(dplyr)
library(ggplot2)

load("C:/Users/mlnpa/Documents/glioblastoma_scrnaseq/gbm_GSE182109_filtered_scores_heterogeneity.RData")


states <- c(
  "State1",
  "State2",
  "State3",
  "State5"
)



all_deg <- bind_rows(
  lapply(states, function(st){

    res <- read.delim(
      paste0("CellState_DEG_", st, "_padj0.1.txt"),
      stringsAsFactors = FALSE
    )

    res$state <- st
    res
  })
)



top_genes <- all_deg %>%
  filter(avg_log2FC > 0) %>%
  group_by(state) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

genes_use <- unique(top_genes$gene)



genes_exclude <- c(
  "SLPI"
)

genes_use <- setdiff(
  genes_use,
  genes_exclude
)


genes_add <- c(
  "COL8A1",
  "IL1B",
  "IFI27",
  "CCL2",
  "CCL3",
  "COL1A2",
  "IGFBP6",
  "FOSB",
  "FOS",
  "FOSL2",
  "JUN",
  "JUNB",
  "BCL2L1"
)

genes_use <- unique(c(
  genes_use,
  genes_add
))



genes_use <- intersect(
  genes_use,
  rownames(se)
)

cat("Selected", length(genes_use), "genes\n")

print(genes_use)



Idents(se) <- "CellState"

se_main <- subset(
  se,
  idents = states
)

se_main$CellState <- factor(
  se_main$CellState,
  levels = states
)

Idents(se_main) <- "CellState"



se_main <- ScaleData(
  se_main,
  features = genes_use,
  verbose = FALSE
)



p <- DoHeatmap(
  se_main,
  features = genes_use,
  group.by = "CellState"
) +
  NoLegend() +
  theme(
    axis.text.y = element_text(
      size = 8,
      face = "italic"
    ),
    axis.text.x = element_text(
      size = 12,
      face = "bold"
    )
  )

p

ggsave(
  "CellState_Top10Markers_Custom_Heatmap_MainFigure1.tiff",
  p,
  device = "tiff",
  width = 12,
  height = 14,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


