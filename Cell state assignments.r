##Senescence based cell state heterogeneity


library(Seurat)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(mclust)

load("C:/Users/mlnpa/Documents/glioblastoma_scrnaseq/gbm_GSE182109_filtered_scores.RData")

df <- se@meta.data[, c(
  "Patient",
  "E2F_TARGETS_HALLMARK.score1",
  "CELL_CYCLE_REGEV.score1",
  "IFN_I_response.score1",
  "SENESCENCE_GOBP.score1",
  "P53_HALLMARK.score1",
  "CELLAGE_SEN.score1",
  "repli_top150.score1",
  "ras_top150.score1",
  "IR_top150.score1",
  "SASP.score1",
  "ISG_nature.score1",
  "Hallmark_IFNg_resp.score1",
  "Hallmark_IFNa_resp.score1",
  "Hallmark_Inflammatory_resp.score1",
  "p16_top150.score1",
  "p15_top150.score1",
  "senmayo.score1"
)]

df <- na.omit(df)

# -------------------------------
#Prepare matrix
# -------------------------------
scores <- df %>% select(-Patient)
scores <- as.data.frame(lapply(scores, as.numeric))

scores_scaled <- scale(scores)
rownames(scores_scaled) <- rownames(df)


set.seed(123)

scores_pca <- prcomp(scores_scaled)$x[,1:10]

kmeans_res <- kmeans(scores_pca, centers = 6, nstart = 50, iter.max = 100)

cell_states_all <- kmeans_res$cluster

# -------------------------------
#State composition
# -------------------------------
annotation_all <- data.frame(
  Patient = df$Patient,
  State = factor(cell_states_all)
)

state_composition <- prop.table(
  table(annotation_all$Patient, annotation_all$State),
  1
)

print(state_composition)

# -------------------------------
#State means (biology)
# -------------------------------
scores_df <- as.data.frame(scores_scaled)
scores_df$State <- cell_states_all

state_means <- scores_df %>%
  group_by(State) %>%
  summarise(across(everything(), mean))

# =========================================
#HEATMAP 
# =========================================

set.seed(123)

# stratified sampling (better than random)
cells_use <- unlist(
  lapply(split(rownames(scores_scaled), cell_states_all), function(x){
    sample(x, min(length(x), 250))
  })
)

mat <- scores_scaled[cells_use, , drop = FALSE]

# similarity
cell_cor <- cor(t(mat), method = "spearman")

# enhance contrast
cell_cor[cell_cor > 0.7] <- 0.7
cell_cor[cell_cor < -0.3] <- -0.3
cell_cor_adj <- sign(cell_cor) * (abs(cell_cor)^1.5)



# annotation (IMPORTANT: use kmeans states)
annotation_col <- data.frame(
  Patient = factor(df[cells_use, "Patient"]),
  State = factor(cell_states_all[cells_use])
)

rownames(annotation_col) <- cells_use

# colors
ann_colors <- list(
  Patient = setNames(
    colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(annotation_col$Patient))),
    levels(annotation_col$Patient)
  ),
  State = setNames(
    brewer.pal(length(unique(annotation_col$State)), "Set2"),
    levels(annotation_col$State)
  )
)



# order within each state to improve visualization
ord <- unlist(lapply(split(1:nrow(annotation_col), annotation_col$State), function(idx){
  sub_cor <- cell_cor_adj[idx, idx]
  hc <- hclust(as.dist(1 - sub_cor))
  idx[hc$order]
}))

cell_cor_adj <- cell_cor_adj[ord, ord]
annotation_col <- annotation_col[ord, ]

# save heatmap


tiff("senescence_heterogeneity_FINAL.tiff",
     units = "in", width = 15, height = 15,
     res = 600, compression = "lzw")
	 
pheatmap(cell_cor_adj,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = annotation_col,
         annotation_row = annotation_col,
         annotation_colors = ann_colors,
         color = colorRampPalette(rev(brewer.pal(11, "BrBG")))(99),
         breaks = seq(-0.7, 0.7, length.out = 100),
         main = "Cell–cell similarity reveals senescence state heterogeneity")

dev.off()

######################################################################################


# Create labeled states
state_labels <- paste0("State", cell_states_all)

# Keep names for safe mapping
names(state_labels) <- rownames(scores_scaled)

# Add to Seurat metadata (safe alignment)
se$CellState <- state_labels[colnames(se)]

# Convert to ordered factor
se$CellState <- factor(se$CellState,
                       levels = paste0("State", 1:6))


print(table(se$CellState))


DimPlot(se, group.by = "CellState", label = TRUE)




y<- rownames(se[["RNA"]]@data)
se <- ScaleData(se, features = y, assay = "RNA", do.center=TRUE, do.scale=TRUE)


cytotoxic_genes <- c(
"GZMA","GZMB","GZMH","PRF1",
"NKG7","IFNG","TNF","CCL5",
"KLRD1","KLRK1","CTSW"
)

immunosuppressive_genes <- c(
"FOXP3","IL10","TGFB1","CTLA4",
"PDCD1","PDCD1LG2","HAVCR2",
"LAG3","CD163","MRC1",
"ARG1","IDO1","CCR4","CCR8"
)


common_genes <- rownames(se[["RNA"]]@scale.data)

cytotoxic_use <- cytotoxic_genes[cytotoxic_genes %in% common_genes]
immunosuppressive_use <- immunosuppressive_genes[immunosuppressive_genes %in% common_genes]


se <- AddModuleScore(se,
                     features = list(cytotoxic_use),
                     slot = "scale.data",
                     name = "Cytotoxic",
                     seed = 1)

se <- AddModuleScore(se,
                     features = list(immunosuppressive_use),
                     slot = "scale.data",
                     name = "Immunosuppressive",
                     seed = 1)


colnames(se@meta.data)[grep("^Cytotoxic1$", colnames(se@meta.data))] <- "Cytotoxic_score"
colnames(se@meta.data)[grep("^Immunosuppressive1$", colnames(se@meta.data))] <- "Immunosuppressive_score"


head(se@meta.data[, c("Cytotoxic_score", "Immunosuppressive_score")])


FeaturePlot(se, features = c("Cytotoxic_score", "Immunosuppressive_score"))



#save(se, file="C:/Users/mlnpa/Documents/glioblastoma_scrnaseq/gbm_GSE182109_filtered_scores_heterogeneity.RData")

###annotate CellState


state_means_mat <- as.data.frame(state_means)


print(state_means_mat)


library(pheatmap)

pheatmap(as.matrix(state_means_mat[,-1]),
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "State-wise signature enrichment")


##########################################

library(Seurat)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(gplots)
library("ggsci")
library("ggplot2")
library("gridExtra")
library(RColorBrewer)
library(ggbreak)
library(ggpubr)

load("C:/Users/mlnpa/Documents/glioblastoma_scrnaseq/gbm_GSE182109_filtered_scores_heterogeneity.RData")


se$expr_CDKN2A <- FetchData(se, vars = "CDKN2A")[,1]


df_plot <- data.frame(
  State = se$CellState,
  IR_Senescence = se$IR_top150.score1,
  p16_associated = se$p16_top150.score1,
  SenMayo = se$senmayo.score1,
  p16_expression = se$expr_CDKN2A
)

signatures <- c(
  "IR_Senescence",
  "p16_associated",
  "SenMayo",
  "p16_expression"
)

library(ggpubr)

df_plot$State <- factor(
  df_plot$State,
  levels = paste0("State", 1:6)
)

dir.create("signature_boxplots", showWarnings = FALSE)

state_cols <- c(
  "#bddf26",
  "#29af7f",
  "#2e6f8e",
  "#DCCB8F",
  "#C97B3C",
  "#A63D2E"
)



for(sig in signatures){

  p <- ggviolin(
    df_plot,
    x = "State",
    y = sig,
    fill = "State",
    palette = state_cols,
    trim = TRUE,
    scale = "width"
  ) +

    geom_boxplot(
      width = 0.12,
      fill = "white",
      color = "black",
      outlier.shape = NA,
      size = 0.7
    ) +

    rotate_x_text(angle = 30) +

    theme_classic() +

    theme(
      axis.text.x = element_text(
        face = "bold",
        size = 10
      ),
      axis.text.y = element_text(
        face = "bold",
        size = 10
      ),
      axis.title = element_text(
        face = "bold",
        size = 11
      ),
      axis.line = element_line(size = 1.1),
      axis.ticks = element_line(size = 1.1),
      axis.ticks.length = unit(0.45, "cm"),
      legend.position = "none"
    ) +

    ggtitle(sig) +
    labs(x = NULL, y = "Score")

  tiff(
    filename = paste0(
      "signature_boxplots/",
      sig,
      "_violin.tiff"
    ),
    units = "in",
    width = 4.2,
    height = 4,
    res = 300,
    compression = "lzw"
  )

  print(p)

  dev.off()
}
  

#########################################################################

library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(rstatix)
library(effectsize)
library(ggpubr)

load("C:/Users/mlnpa/Documents/glioblastoma_scrnaseq/gbm_GSE182109_filtered_scores_heterogeneity.RData")


se$expr_CDKN2A <- FetchData(se, vars = "CDKN2A")[,1]
se$CellState <- factor(se$CellState, levels = paste0("State", 1:6))

df <- se@meta.data


state_colors <- setNames(brewer.pal(6, "Set2"), levels(se$CellState))

features <- c("expr_CDKN2A",
  "E2F_TARGETS_HALLMARK.score1",
  "CELL_CYCLE_REGEV.score1",
  "IFN_I_response.score1",
  "SENESCENCE_GOBP.score1",
  "P53_HALLMARK.score1",
  "CELLAGE_SEN.score1",
  "repli_top150.score1",
  "ras_top150.score1",
  "IR_top150.score1",
  "SASP.score1",
  "ISG_nature.score1",
  "Hallmark_IFNg_resp.score1",
  "Hallmark_IFNa_resp.score1",
  "Hallmark_Inflammatory_resp.score1",
  "p16_top150.score1",
  "p15_top150.score1",
  "senmayo.score1"
)

# =========================================
#Statistical pipeline
# =========================================
stats_list <- list()

for (f in features) {

  form <- as.formula(paste(f, "~ CellState"))

  # Kruskal-Wallis (global test)
  kw <- rstatix::kruskal_test(df, formula = form)

  # Effect size (epsilon squared)
  model <- aov(df[[f]] ~ df$CellState)
  eff <- suppressMessages(effectsize::epsilon_squared(model))

  # Dunn posthoc test (BH correction)
  dunn <- rstatix::dunn_test(
    df,
    formula = form,
    p.adjust.method = "BH"
  )

  stats_list[[f]] <- list(
    kw = kw,
    effect = eff,
    dunn = dunn
  )
}


stats_list$expr_CDKN2A$kw
stats_list$expr_CDKN2A$effect

stats_list$p16_top150.score1$kw
stats_list$p16_top150.score1$effect

stats_list$senmayo.score1$kw
stats_list$senmayo.score1$effect

stats_list$IR_top150.score1$kw
stats_list$IR_top150.score1$effect


# =========================================
#GRADIENT ANALYSIS 
# =========================================

# -------------------------------
#PCA on signature scores
# -------------------------------
pca_full <- prcomp(scores_scaled)

pc1 <- pca_full$x[,1]

# -------------------------------
#Correlation of each signature with PC1
# -------------------------------
sig_cor <- apply(scores_scaled, 2, function(x){
  cor(x, pc1, method = "spearman")
})

sig_cor <- sort(sig_cor)

# -------------------------------
#Plot signature ordering along gradient
# -------------------------------
png("FigS3_signature_gradient.png", 900, 700)

par(mar = c(12, 4, 4, 2))  # 🔥 increase bottom margin

barplot(sig_cor,
        las = 2,              # vertical labels
        cex.names = 0.7,      # shrink text
        col = "steelblue",
        main = "Signature association with PC1",
        ylab = "Spearman correlation with PC1")

abline(h = 0, lty = 2)

dev.off()


df_plot <- data.frame(
  PC1 = pc1,
  State = factor(cell_states_all)
)

png("FigS3_PC1_distribution.png", 700, 600)

boxplot(PC1 ~ State,
        data = df_plot,
        col = "lightblue",
        main = "Cell states distributed along transcriptional gradient",
        ylab = "PC1 score")

dev.off()


key_sigs <- c(
  "E2F_TARGETS_HALLMARK.score1",   # proliferation
  "p16_top150.score1",
"senmayo.score1", 
  "ISG_nature.score1"            
)

df_long2 <- df_long %>% filter(Signature %in% key_sigs)



# =========================================
#HETEROGENEITY
# =========================================

calc_entropy <- function(p){
  p <- p[p > 0]
  -sum(p * log(p))
}

entropy_scores <- apply(state_composition, 1, calc_entropy)

group <- ifelse(grepl("LGG", names(entropy_scores)), "LGG", "GBM")

t.test(entropy_scores ~ group)

png("FigS3_entropy.png", 600, 600)
boxplot(entropy_scores ~ group, col = c("skyblue","salmon"))
stripchart(entropy_scores ~ group, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 16)
dev.off()



comp_df <- as.data.frame(state_composition)
comp_df$Patient <- comp_df$Var1


#png("FigS3_composition.png", 800, 600)
tiff('FigS3_composition.tiff', units="in", width=6, height=5, res=900, compression = 'lzw')
ggplot(comp_df, aes(x = Patient, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity") +


theme_classic(base_size = 14) +
theme(
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(size = 1.05),
    axis.ticks = element_line(size = 1.05),
    axis.ticks.length = unit(0.5, "cm"),
    legend.position = "right"
  ) + labs(fill = "Cell state")
  
  dev.off()
  
  

################################
##intra tumor heterogeneity

library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(ggplot2)


scores_df <- as.data.frame(scores_scaled)
scores_df$Patient <- df$Patient
scores_df$State <- factor(cell_states_all) 

# -------------------------------
#Intra-tumoral Shannon entropy per patient
# -------------------------------
calc_entropy <- function(p){
  p <- p[p > 0]  
  -sum(p * log(p))
}


state_comp <- prop.table(table(scores_df$Patient, scores_df$State), 1)

# entropy per patient
entropy_scores <- apply(state_comp, 1, calc_entropy)
entropy_scores

# Optional: add group info for GBM vs LGG
group <- ifelse(grepl("LGG", names(entropy_scores)), "LGG", "GBM")

# Test difference
t_test_res <- t.test(entropy_scores ~ group)
print(t_test_res)


tiff("violin_states_per_patient.tiff",
     units = "in", width = 16, height = 6,
     res = 600, compression = "lzw")

ggplot(scores_df, aes(x = Patient, y = as.numeric(State), fill = State)) +
  geom_violin(scale = "width", trim = FALSE, color = "black") +  # violin border
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Patient", y = "Transcriptional State", 
       title = "Distribution of transcriptional states across samples") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid.major = element_blank(),  # remove major grid
    panel.grid.minor = element_blank(),  # remove minor grid
    panel.border = element_rect(color = "black", fill = NA, size = 1.2)
  )

dev.off()

#intra tumor heterogeneity

# state composition per patient
state_composition <- prop.table(
  table(annotation_all$Patient, annotation_all$State), 1
)

# entropy function
calc_entropy <- function(p){
  p <- p[p > 0]  # avoid log(0)
  -sum(p * log(p))
}

# intra-tumor entropy
entropy_scores <- apply(state_composition, 1, calc_entropy)
entropy_scores


library(ggplot2)

# Prepare data for violin plot
scores_df$Patient <- df$Patient  # make sure patient info is there

# Calculate per-patient entropy (already have)
entropy_df <- data.frame(
  Patient = names(entropy_scores),
  Entropy = entropy_scores
)



# =========================================
#SUBSAMPLING ROBUSTNESS (50 RUNS)
# =========================================

get_clusters_subsample <- function(seed){
  set.seed(seed)
  
  # stratified sampling
  cells_use <- unlist(
    lapply(split(rownames(scores_scaled), cell_states_all), function(x){
      sample(x, min(length(x), 250))
    })
  )
  
  mat <- scores_scaled[cells_use, ]
  
  # PCA
  pca <- prcomp(mat)$x[,1:10]
  
  # kmeans
  km <- kmeans(pca, centers = 6, nstart = 30)
  
  # store cluster assignment
  out <- rep(NA, nrow(scores_scaled))
  names(out) <- rownames(scores_scaled)
  out[cells_use] <- km$cluster
  
  return(out)
}

#50 subsampling runs
cluster_list <- lapply(1:50, get_clusters_subsample)

# -------------------------------
#Pairwise ARI
# -------------------------------
aris <- c()

for(i in 1:49){
  for(j in (i+1):50){
    
    common <- intersect(
      names(cluster_list[[i]][!is.na(cluster_list[[i]])]),
      names(cluster_list[[j]][!is.na(cluster_list[[j]])])
    )
    
    cl1 <- cluster_list[[i]][common]
    cl2 <- cluster_list[[j]][common]
    
    aris <- c(aris, adjustedRandIndex(cl1, cl2))
  }
}

# -------------------------------
# Summary
# -------------------------------
cat("Mean ARI:", mean(aris), "\n")
cat("SD ARI:", sd(aris), "\n")

# -------------------------------
# Plot
# -------------------------------
png("FigS3_subsampling_ARI.png", 700, 600)

hist(aris,
     breaks = 20,
     col = "grey",
     border = "black",
     main = "Clustering stability across subsampling",
     xlab = "Adjusted Rand Index (ARI)")

abline(v = mean(aris), col = "red", lwd = 2)

dev.off()



# =========================================
#CONSENSUS HEATMAP
# =========================================

get_clusters_subsample <- function(seed){
  set.seed(seed)
  
  cells_sub <- unlist(
    lapply(split(rownames(scores_scaled), cell_states_all), function(x){
      sample(x, min(length(x), 250))
    })
  )
  
  mat_sub <- scores_scaled[cells_sub, ]
  pca <- prcomp(mat_sub)$x[,1:10]
  km <- kmeans(pca, centers = 6, nstart = 20)
  
  out <- rep(NA, nrow(scores_scaled))
  names(out) <- rownames(scores_scaled)
  out[cells_sub] <- km$cluster
  
  return(out)
}

# -------------------------------
# 50 subsampling runs
# -------------------------------
cluster_list <- lapply(1:50, get_clusters_subsample)

# -------------------------------
# SAME cells as main heatmap
# -------------------------------
n <- length(cells_use)

consensus_mat <- matrix(0, n, n)
counts_mat <- matrix(0, n, n)

rownames(consensus_mat) <- cells_use
colnames(consensus_mat) <- cells_use

# -------------------------------
# Build consensus matrix
# -------------------------------
for(run in 1:50){
  
  cl <- cluster_list[[run]][cells_use]
  valid <- !is.na(cl)
  
  idx <- which(valid)
  
  for(i in idx){
    for(j in idx){
      
      counts_mat[i,j] <- counts_mat[i,j] + 1
      
      if(cl[i] == cl[j]){
        consensus_mat[i,j] <- consensus_mat[i,j] + 1
      }
    }
  }
}

# normalize
consensus_mat <- consensus_mat / counts_mat
consensus_mat[is.na(consensus_mat)] <- 0

# -------------------------------
#MATCH SAME ORDER AS MAIN HEATMAP
# -------------------------------
consensus_mat <- consensus_mat[ord, ord]

consensus_scaled <- consensus_mat


tiff("senescence_heterogeneity_consensus_FINAL.tiff",
     units = "in", width = 15, height = 15,
     res = 600, compression = "lzw")
	 
	 pheatmap(consensus_scaled,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = annotation_col,
         annotation_row = annotation_col,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("white", "red"))(100),
         breaks = seq(0, 1, length.out = 100),
         main = "Consensus cell–cell similarity across 50 subsampling runs")
		 
		 dev.off()
		 



# =========================================
#K SELECTION (SUPPLEMENT)
# =========================================


wss <- c()

for (k in 1:10) {
  km <- kmeans(scores_scaled, centers = k, nstart = 20)
  wss[k] <- km$tot.withinss
}

png("FigS3_elbow.png", 600, 600)
plot(1:10, wss, type = "b",
     xlab = "Number of clusters (k)",
     ylab = "Total within-cluster sum of squares",
     main = "Elbow method for k selection")

dev.off()


sil_scores <- c()

for (k in 2:10) {
  
  km <- kmeans(pca, centers = k, nstart = 20)
  
  idx <- sample(1:nrow(pca), min(1000, nrow(pca)))
  
  d <- dist(pca[idx, ])
  sil <- silhouette(km$cluster[idx], d)
  
  sil_scores[k] <- mean(sil[, 3])
}



png("FigS3_silhouette.png", 600, 600)
plot(2:10, sil_scores[2:10], type = "b")
dev.off()

# =========================================
#SIGNATURE CORRELATION
# =========================================

png("FigS3_correlation.png", 600, 600)
heatmap(cor(scores_scaled, method = "spearman"))
dev.off()





