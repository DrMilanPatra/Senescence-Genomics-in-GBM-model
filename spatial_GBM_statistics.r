
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmnet)
library(FNN)
library(spdep)


data_dir <- "C:/Users/mlnpa/Downloads/GSE194329_RAW"

heatmap_dir <- file.path(data_dir,"VISIUM_HEATMAPS")
if(!dir.exists(heatmap_dir)) dir.create(heatmap_dir)



samples <- c(
"GBM1_spaceranger_out",
"GBM2_spaceranger_out",
"GBM3_spaceranger_out",
"GBM4_spaceranger_out",
"GBM5_1_spaceranger_out",
"GBM5_2_spaceranger_out"
)



beta <- read.table(
"genesets_milan/intg_all_ittai_gbm_spatial.txt",
sep="\t",
header=TRUE,
stringsAsFactors=FALSE
)

beta_gene <- list()

for(m in colnames(beta)){

genes <- beta[[m]]
genes <- genes[genes != ""]
genes <- unique(genes)

beta_gene[[m]] <- genes

}



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



results <- data.frame()

elastic_results <- data.frame()

all_moran_results <- data.frame()



for(s in samples){

cat("\n===========================\n")
cat("Processing sample:",s,"\n")
cat("===========================\n")



visium <- Load10X_Spatial(
data.dir=file.path(data_dir,s)
)



visium <- SCTransform(visium,assay="Spatial",verbose=FALSE)



for(gs in names(beta_gene)){

genes <- intersect(beta_gene[[gs]], rownames(visium))

if(length(genes) >= 3){

visium <- AddModuleScore(
visium,
features=list(genes),
name=gs
)

}

}



visium <- AddModuleScore(
visium,
features=list(intersect(cytotoxic_genes,rownames(visium))),
name="Cytotox"
)

visium <- AddModuleScore(
visium,
features=list(intersect(immunosuppressive_genes,rownames(visium))),
name="ImmuneSupp"
)



coords <- GetTissueCoordinates(visium)
meta <- visium@meta.data

df <- cbind(coords,meta)



score_cols <- grep("1$",colnames(df),value=TRUE)
df[score_cols] <- scale(df[score_cols])

############################
# Exact column detection
############################

IFN_IFNG <- grep("IFN_IFNG_gbm1",colnames(df),value=TRUE)
IFN_HALL <- grep("IFN_Hallmark_IFNg_resp1",colnames(df),value=TRUE)
ISG <- grep("ISG_nature1",colnames(df),value=TRUE)

SEN_IR <- grep("Sen_IR_top1501",colnames(df),value=TRUE)
SEN_MAYO <- grep("Sen_senmayo1",colnames(df),value=TRUE)
SEN_P16 <- grep("Sen_p16_top1501",colnames(df),value=TRUE)

SASP <- grep("^SASP1$",colnames(df),value=TRUE)

CYT <- "Cytotox1"
SUP <- "ImmuneSupp1"

############################
# Build interaction terms
############################

df$SASP_SENIR <- df[[SASP]] * df[[SEN_IR]]
df$SASP_SENMAYO <- df[[SASP]] * df[[SEN_MAYO]]
df$SASP_SENP16 <- df[[SASP]] * df[[SEN_P16]]

df$IFN_SENIR <- df[[IFN_IFNG]] * df[[SEN_IR]]
df$IFN_SENMAYO <- df[[IFN_IFNG]] * df[[SEN_MAYO]]
df$IFN_SENP16 <- df[[IFN_IFNG]] * df[[SEN_P16]]

############################
# Elastic net predictors
############################

predictors <- c(
SEN_IR,SEN_MAYO,SEN_P16,
IFN_IFNG,IFN_HALL,ISG,
SASP,
"SASP_SENIR","SASP_SENMAYO","SASP_SENP16",
"IFN_SENIR","IFN_SENMAYO","IFN_SENP16"
)

X <- as.matrix(df[,predictors])

############################
# Responses
############################

y_cyt <- df[[CYT]]
y_sup <- df[[SUP]]

############################
# Elastic net
############################

cv_cyt <- cv.glmnet(X,y_cyt,alpha=0.5)
cv_sup <- cv.glmnet(X,y_sup,alpha=0.5)

fit_cyt <- glmnet(X,y_cyt,alpha=0.5,lambda=cv_cyt$lambda.min)
fit_sup <- glmnet(X,y_sup,alpha=0.5,lambda=cv_sup$lambda.min)

############################
# Predictions
############################

pred_cyt <- predict(fit_cyt,X)
pred_sup <- predict(fit_sup,X)

df$Cytotox_pred <- as.numeric(pred_cyt)
df$Supp_pred <- as.numeric(pred_sup)

############################
# Add predictions to Seurat
############################

visium$Cytotox_pred <- df$Cytotox_pred
visium$Supp_pred <- df$Supp_pred

############################
# Model performance
############################

R2_cyt <- cor(pred_cyt,y_cyt)^2
R2_sup <- cor(pred_sup,y_sup)^2



# Features to plot
features_to_plot <- c(
  SEN_IR, IFN_IFNG, ISG,
  SEN_MAYO, SEN_P16, SASP,
  "Cytotox_pred", "Supp_pred"
)

# Generate individual plots
plots <- lapply(features_to_plot, function(f) {
  SpatialFeaturePlot(visium, features=f)
})

# Combine plots in a grid
# Here, 4 columns and 2 rows for 8 features
combined <- wrap_plots(plots, ncol=4)

# Save
ggsave(
  file.path(heatmap_dir, paste0(s,"_spatial_summary_expanded.png")),
  combined,
  width=16,
  height=8
)

############################
# Moran statistics
############################

coords_mat <- as.matrix(df[,c("imagecol","imagerow")])
neighbors <- knn2nb(knearneigh(coords_mat, k=6))
weights <- nb2listw(neighbors)

moran_features <- c(
  "Sen_IR_top1501",
  "Sen_p16_top1501",
  "Sen_senmayo1",
  "IFN_IFNG_gbm1",
  "IFN_Hallmark_IFNg_resp1",
  "ISG_nature1",
  "Cytotox1",
  "ImmuneSupp1"
)

moran_list <- lapply(moran_features, function(f){

  if(!f %in% colnames(df)) return(NULL)

  res <- moran.test(df[[f]], weights)

  data.frame(
    sample = s,
    feature = f,
    Moran_I = unname(res$estimate[1]),
    Moran_p = res$p.value
  )
})

moran_results <- do.call(rbind, moran_list)

all_moran_results <- rbind(
  all_moran_results,
  moran_results
)

write.csv(
  moran_results,
  file.path(data_dir, paste0(s, "_moran_results.csv")),
  row.names = FALSE
)

############################
# SUMMARY
############################
results <- rbind(results, data.frame(
  sample = s,
  Moran_Sen_mean = mean(moran_results$Moran_I[grepl("Sen", moran_results$feature)]),
  Moran_IFN_mean = mean(moran_results$Moran_I[grepl("IFN|ISG", moran_results$feature)]),
  Moran_overall_mean = mean(moran_results$Moran_I)
))



elastic_results <- rbind(elastic_results,data.frame(

sample=s,
R2_cytotoxic=R2_cyt,
R2_suppressive=R2_sup

))


write.csv(
df,
file.path(data_dir,paste0(s,"_spot_results.csv")),
row.names=FALSE
)


}


write.csv(
  results,
  file.path(data_dir,"SPATIAL_STATISTICS_RESULTS.csv"),
  row.names=FALSE
)

write.csv(
elastic_results,
file.path(data_dir,"ELASTIC_NET_RESULTS.csv"),
row.names=FALSE
)



write.csv(
  all_moran_results,
  file.path(data_dir,"ALL_MORAN_RESULTS.csv"),
  row.names = FALSE
)




cat("\nAnalysis complete\n")

############################
# Spatial Lag Regression
# All Samples
# Senescence × IFN model
############################

rm(list = ls())

library(spdep)
library(spatialreg)
library(dplyr)

############################################################
# Input
############################################################

data_dir <- "C:/Users/mlnpa/Downloads/GSE194329_RAW"

samples <- c(
  "GBM1_spaceranger_out",
  "GBM2_spaceranger_out",
  "GBM3_spaceranger_out",
  "GBM4_spaceranger_out",
  "GBM5_1_spaceranger_out",
  "GBM5_2_spaceranger_out"
)

outcomes <- c(
  "Cytotox1",
  "ImmuneSupp1"
)

senescence <- c(
  "Sen_CELLAGE_SEN1",
  "Sen_IR_top1501",
  "Sen_p16_top1501",
  "Sen_senmayo1",
  "Sen_repli_top1501",
  "Sen_ras_top1501"
)



spatial_list <- list()



for(s in samples){

  cat("\n====================================\n")
  cat("Processing:", s, "\n")
  cat("====================================\n")

 

  df <- read.csv(
    file.path(
      data_dir,
      paste0(s, "_spot_results.csv")
    )
  )

  ##########################################################
  # Spatial neighborhood graph
  ##########################################################

  coords_mat <- as.matrix(
    df[, c("imagecol", "imagerow")]
  )

  neighbors <- knn2nb(
    knearneigh(
      coords_mat,
      k = 6
    )
  )

  lw <- nb2listw(
    neighbors,
    style = "W",
    zero.policy = TRUE
  )

 

  sample_list <- list()

  ##########################################################
  # Run spatial lag models
  ##########################################################

  for(outcome in outcomes){

    for(sen in senescence){

      cat(
        "Running:",
        outcome,
        "|",
        sen,
        "\n"
      )

      formula_2way <- as.formula(
        paste0(
          outcome,
          " ~ ",
          sen,
          " * IFN_IFNG_gbm1"
        )
      )

      fit <- tryCatch(

        lagsarlm(
          formula_2way,
          data = df,
          listw = lw,
          zero.policy = TRUE
        ),

        error = function(e){

          cat(
            "FAILED:",
            outcome,
            sen,
            "\n"
          )

          return(NULL)
        }
      )

      if(is.null(fit))
        next

      coefs <- summary(fit)$Coef

      res <- data.frame(

        sample = s,

        outcome = outcome,

        senescence = sen,

        term = rownames(coefs),

        estimate = coefs[,1],

        std_error = coefs[,2],

        z_value = coefs[,3],

        p_value = coefs[,4],

        rho = fit$rho,

        stringsAsFactors = FALSE

      )

   

      sample_list[[length(sample_list) + 1]] <- res

      spatial_list[[length(spatial_list) + 1]] <- res

    }
  }



  sample_results <- bind_rows(sample_list)

  write.csv(
    sample_results,
    file.path(
      data_dir,
      paste0(
        s,
        "_SPATIAL_LAG_RESULTS.csv"
      )
    ),
    row.names = FALSE
  )

  cat(
    "Saved:",
    paste0(
      s,
      "_SPATIAL_LAG_RESULTS.csv"
    ),
    "\n"
  )
}



spatial_results <- bind_rows(spatial_list)

write.csv(
  spatial_results,
  file.path(
    data_dir,
    "SPATIAL_LAG_REGRESSION_ALL_SAMPLES.csv"
  ),
  row.names = FALSE
)



sig_interactions <- spatial_results %>%
  filter(
    grepl(":", term),
    p_value < 0.05
  ) %>%
  arrange(
    outcome,
    senescence,
    p_value
  )

write.csv(
  sig_interactions,
  file.path(
    data_dir,
    "SPATIAL_LAG_SIGNIFICANT_INTERACTIONS.csv"
  ),
  row.names = FALSE
)

#############ploting##########################################

library(dplyr)
library(ggplot2)
library(patchwork)
data_dir <- "C:/Users/mlnpa/Downloads/GSE194329_RAW"

plot_df <- read.csv(
  file.path(
    data_dir,
    "SPATIAL_LAG_SIGNIFICANT_INTERACTIONS.csv"
  )
)



plot_df$senescence <- recode(
  plot_df$senescence,
  "Sen_senmayo1"      = "SenMayo",
  "Sen_p16_top1501"   = "p16",
  "Sen_IR_top1501"    = "IR-induced",
  "Sen_CELLAGE_SEN1"  = "CELLAGE",
  "Sen_repli_top1501" = "Replicative",
  "Sen_ras_top1501"   = "RAS-induced"
)

plot_df$senescence <- factor(
  plot_df$senescence,
  levels = rev(c(
    "SenMayo",
    "p16",
    "IR-induced",
    "CELLAGE",
    "Replicative",
    "RAS-induced"
  ))
)



sen_colors <- c(
  "SenMayo"      = "#2e6f8e",
  "p16"          = "#29af7f",
  "IR-induced"   = "#bddf26",
  "CELLAGE"      = "#6ec5b8",
  "Replicative"  = "#DCCB8F",
  "RAS-induced"  = "#C97B3C"
)

###################################################
# Cytotoxic panel
###################################################

p_cyt <- ggplot(
  subset(plot_df, outcome=="Cytotox1"),
  aes(
    x = estimate,
    y = senescence,
    fill = senescence
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey40",
    size = 0.8
  ) +

  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    colour = "black",
    size = 0.7
  ) +

  geom_jitter(
  shape = 21,
  aes(fill = senescence),
  colour = "black",      # black outline
  stroke = 0.6,          # border thickness
  size = 3.8,            # larger points
  alpha = 0.6,
  width = 0,
  height = 0.12
) +

  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 4
  ) +

  scale_fill_manual(values = sen_colors) +

  theme_classic(base_size = 14) +

  theme(
    axis.text = element_text(face="bold"),
    axis.title = element_text(face="bold"),
    axis.line = element_line(size=1.1),
    axis.ticks = element_line(size=1.1),
    axis.ticks.length = unit(0.35,"cm"),
    legend.position = "none",
    plot.title = element_text(
      face="bold",
      hjust=0.5
    )
  ) +

  labs(
    title = "Cytotoxic niche",
    x = expression(beta~"(IFN"*gamma*" × Senescence)"),
    y = NULL
  )

###################################################
# Immune suppression panel
###################################################

p_sup <- ggplot(
  subset(plot_df, outcome=="ImmuneSupp1"),
  aes(
    x = estimate,
    y = senescence,
    fill = senescence
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey40",
    size = 0.8
  ) +

  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    colour = "black",
    size = 0.7
  ) +

  geom_jitter(
  shape = 21,
  aes(fill = senescence),
  colour = "black",      # black outline
  stroke = 0.6,          # border thickness
  size = 3.8,            # larger points
  alpha = 0.6,
  width = 0,
  height = 0.12
) +

  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 4
  ) +

  scale_fill_manual(values = sen_colors) +

  theme_classic(base_size = 14) +

  theme(
    axis.text = element_text(face="bold"),
    axis.title = element_text(face="bold"),
    axis.line = element_line(size=1.1),
    axis.ticks = element_line(size=1.1),
    axis.ticks.length = unit(0.35,"cm"),
    legend.position = "none",
    plot.title = element_text(
      face="bold",
      hjust=0.5
    )
  ) +

  labs(
    title = "Immunosuppressive niche",
    x = expression(beta~"(IFN"*gamma*" × Senescence)"),
    y = NULL
  )

###################################################
# Combine
###################################################

ggsave(
  file.path(data_dir, "spatial_lag_cyto.png"),
  p_cyt,
  width = 5,
  height = 5
)

ggsave(
  file.path(data_dir, "spatial_lag_suppress.png"),
  p_sup,
  width = 5,
  height = 5
)



############################################################################################

##Cellstate specific spatial transcriptomics analysis

rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmnet)
library(FNN)
library(spdep)

data_dir <- "C:/Users/mlnpa/Downloads/GSE194329_RAW"

heatmap_dir <- file.path(data_dir,"VISIUM_HEATMAPS")
if(!dir.exists(heatmap_dir)) dir.create(heatmap_dir)

samples <- c(
  "GBM1_spaceranger_out",
  "GBM2_spaceranger_out",
  "GBM3_spaceranger_out",
  "GBM4_spaceranger_out",
  "GBM5_1_spaceranger_out",
  "GBM5_2_spaceranger_out"
)


results <- data.frame()
elastic_results <- data.frame()
all_moran_results <- data.frame()


for(s in samples){

  cat("\n===========================\n")
  cat("Processing sample:", s, "\n")
  cat("===========================\n")


  visium <- Load10X_Spatial(
    data.dir = file.path(data_dir, s)
  )

  visium <- SCTransform(visium, assay="Spatial", verbose=FALSE)


  beta <- read.table(
    "genesets_milan/state_up_markers.txt",
    sep="\t",
    header=TRUE,
    stringsAsFactors=FALSE
  )

  beta_gene <- list()
  for(m in colnames(beta)){
    genes <- unique(beta[[m]])
    genes <- genes[genes != ""]
    beta_gene[[m]] <- genes
  }


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

  # -------------------------
  # Module scores
  # -------------------------
  for(gs in names(beta_gene)){
    genes <- intersect(beta_gene[[gs]], rownames(visium))

    if(length(genes) >= 3){
      visium <- AddModuleScore(
        visium,
        features = list(genes),
        name = gs
      )
    }
  }

  visium <- AddModuleScore(
    visium,
    features = list(intersect(cytotoxic_genes, rownames(visium))),
    name = "Cytotox"
  )

  visium <- AddModuleScore(
    visium,
    features = list(intersect(immunosuppressive_genes, rownames(visium))),
    name = "ImmuneSupp"
  )

  # -------------------------
  # Build dataframe
  # -------------------------
  coords <- GetTissueCoordinates(visium)
  df <- cbind(coords, visium@meta.data)

  # -------------------------
  # STANDARDIZE SCORES
  # -------------------------
  score_cols <- grep("1$", colnames(df), value=TRUE)
  df[score_cols] <- scale(df[score_cols])

  # -------------------------
  # DEFINE KEY FEATURES
  # -------------------------
  ST2 <- "stup21"
  ST5 <- "stup51"
  IFN_IFNG <- "IFN_IFNG_gbm1"

  stopifnot(all(c(ST2, ST5, IFN_IFNG) %in% colnames(df)))

  CYT <- "Cytotox1"
  SUP <- "ImmuneSupp1"

  # -------------------------
  # INTERACTIONS (ONLY STATE × IFN)
  # -------------------------
  df$ST2_IFN <- df[[ST2]] * df[[IFN_IFNG]]
  df$ST5_IFN <- df[[ST5]] * df[[IFN_IFNG]]

  # -------------------------
  # MODEL MATRIX
  # -------------------------
  predictors <- c(
    ST2,
    ST5,
    IFN_IFNG,
    "ST2_IFN",
    "ST5_IFN"
  )

  X <- scale(as.matrix(df[, predictors]))

  y_cyt <- df[[CYT]]
  y_sup <- df[[SUP]]

  # -------------------------
  # ELASTIC NET
  # -------------------------
  cv_cyt <- cv.glmnet(X, y_cyt, alpha=0.5)
  cv_sup <- cv.glmnet(X, y_sup, alpha=0.5)

  fit_cyt <- glmnet(X, y_cyt, alpha=0.5, lambda=cv_cyt$lambda.min)
  fit_sup <- glmnet(X, y_sup, alpha=0.5, lambda=cv_sup$lambda.min)

  pred_cyt <- as.numeric(predict(fit_cyt, X))
  pred_sup <- as.numeric(predict(fit_sup, X))

  df$Cytotox_pred <- pred_cyt
  df$Supp_pred <- pred_sup

  visium$Cytotox_pred <- df$Cytotox_pred
  visium$Supp_pred <- df$Supp_pred

  # -------------------------
  # R2
  # -------------------------
  R2_cyt <- cor(pred_cyt, y_cyt)^2
  R2_sup <- cor(pred_sup, y_sup)^2

  # -------------------------
  # SPATIAL MORAN
  # -------------------------
  coords_mat <- as.matrix(df[,c("imagecol","imagerow")])
  neighbors <- knn2nb(knearneigh(coords_mat, k=6))
  weights <- listw <- nb2listw(neighbors)

  moran_features <- c(ST2, ST5, IFN_IFNG, CYT, SUP)

  moran_results <- do.call(rbind, lapply(moran_features, function(f){

    res <- moran.test(df[[f]], weights)

    data.frame(
      sample = s,
      feature = f,
      Moran_I = unname(res$estimate[1]),
      Moran_p = res$p.value
    )
  }))

  all_moran_results <- rbind(all_moran_results, moran_results)

  # -------------------------
  # SUMMARY
  # -------------------------
  results <- rbind(results, data.frame(
    sample = s,
    Moran_mean = mean(moran_results$Moran_I),
    Moran_IFN = mean(moran_results$Moran_I[grepl(IFN_IFNG, moran_results$feature)]),
    R2_cyt = R2_cyt,
    R2_sup = R2_sup
  ))

  elastic_results <- rbind(elastic_results, data.frame(
    sample = s,
    R2_cytotoxic = R2_cyt,
    R2_suppressive = R2_sup
  ))


  write.csv(df,
            file.path(data_dir, paste0(s,"_spot_cellstate_results.csv")),
            row.names=FALSE)

  write.csv(moran_results,
            file.path(data_dir, paste0(s,"_moran_cellstate_results.csv")),
            row.names=FALSE)

}


write.csv(results,
          file.path(data_dir,"SPATIAL_STATISTICS_cellstate_RESULTS.csv"),
          row.names=FALSE)

write.csv(elastic_results,
          file.path(data_dir,"ELASTIC_NET_cellstate_RESULTS.csv"),
          row.names=FALSE)

write.csv(all_moran_results,
          file.path(data_dir,"ALL_MORAN_cell_state_RESULTS.csv"),
          row.names=FALSE)

cat("\nAnalysis complete\n")


############################
# Spatial Lag Regression

############################

rm(list = ls())

library(spdep)
library(spatialreg)
library(dplyr)


data_dir <- "C:/Users/mlnpa/Downloads/GSE194329_RAW"

samples <- c(
  "GBM1_spaceranger_out",
  "GBM2_spaceranger_out",
  "GBM3_spaceranger_out",
  "GBM4_spaceranger_out",
  "GBM5_1_spaceranger_out",
  "GBM5_2_spaceranger_out"
)

outcomes <- c(
  "Cytotox1",
  "ImmuneSupp1"
)


states <- c("stup21", "stup51")


spatial_list <- list()



for(s in samples){

  cat("\n====================================\n")
  cat("Processing:", s, "\n")
  cat("====================================\n")



  df <- read.csv(
    file.path(data_dir, paste0(s, "_spot_cellstate_results.csv"))
  )



  required_cols <- c(
    states,
    "IFN_IFNG_gbm1",
    "imagecol",
    "imrow"
  )

  
  required_cols <- c(
    states,
    "IFN_IFNG_gbm1",
    "imagecol",
    "imagerow"
  )

  stopifnot(all(required_cols %in% colnames(df)))

  ##########################################################
  # Spatial neighborhood graph
  ##########################################################

  coords_mat <- as.matrix(df[, c("imagecol", "imagerow")])

  neighbors <- knn2nb(
    knearneigh(coords_mat, k = 6)
  )

  lw <- nb2listw(
    neighbors,
    style = "W",
    zero.policy = TRUE
  )



  sample_list <- list()

  ##########################################################
  # SPATIAL LAG MODELS
  ##########################################################

  for(outcome in outcomes){

    cat("\nRunning outcome:", outcome, "\n")

    for(state in states){

      cat("  State:", state, "\n")

      formula_lag <- as.formula(
        paste0(
          outcome,
          " ~ ",
          state,
          " * IFN_IFNG_gbm1"
        )
      )

      fit <- tryCatch(
        lagsarlm(
          formula_lag,
          data = df,
          listw = lw,
          zero.policy = TRUE
        ),
        error = function(e){
          cat("FAILED:", outcome, state, "\n")
          return(NULL)
        }
      )

      if(is.null(fit)) next

      coefs <- summary(fit)$Coef

      res <- data.frame(
        sample = s,
        outcome = outcome,
        state = state,
        term = rownames(coefs),
        estimate = coefs[,1],
        std_error = coefs[,2],
        z_value = coefs[,3],
        p_value = coefs[,4],
        rho = fit$rho,
        stringsAsFactors = FALSE
      )

      sample_list[[length(sample_list) + 1]] <- res
      spatial_list[[length(spatial_list) + 1]] <- res
    }
  }



  sample_results <- bind_rows(sample_list)

  write.csv(
    sample_results,
    file.path(data_dir, paste0(s, "_SPATIAL_LAG_STATE_IFN.csv")),
    row.names = FALSE
  )

  cat("Saved sample:", s, "\n")
}



spatial_results <- bind_rows(spatial_list)

write.csv(
  spatial_results,
  file.path(data_dir, "SPATIAL_LAG_STATE_IFN_ALL_SAMPLES.csv"),
  row.names = FALSE
)



sig_interactions <- spatial_results %>%
  filter(
    grepl("IFN_IFNG_gbm1", term),
    p_value < 0.05
  ) %>%
  arrange(outcome, state, p_value)

write.csv(
  sig_interactions,
  file.path(data_dir, "SPATIAL_LAG_SIGNIFICANT_STATE_IFN.csv"),
  row.names = FALSE
)

cat("\nSpatial lag analysis complete\n")

##ploting

library(dplyr)
library(ggplot2)
library(patchwork)

plot_df <- read.csv(
  file.path(
    data_dir,
    "SPATIAL_LAG_STATE_IFN_ALL_SAMPLES.csv"
  )
)
		  
plot_df <- plot_df %>%
  filter(
    p_value < 0.05) %>%
  filter(
    term %in% c(
      "stup21",
      "stup51",
      "stup21:IFN_IFNG_gbm1",
      "stup51:IFN_IFNG_gbm1"
    )
  )
  
  
  plot_df$effect <- dplyr::recode(
  plot_df$term,
  "stup21" = "State2",
  "stup51" = "State5",
  "stup21:IFN_IFNG_gbm1" = "State2 × IFN",
  "stup51:IFN_IFNG_gbm1" = "State5 × IFN"
)

plot_df$effect <- factor(
  plot_df$effect,
  levels = c(
    "State2",
    "State5",
    "State2 × IFN",
    "State5 × IFN"
  )
)

sen_colors <- c(
  "State5"        = "#2e6f8e",
  "State2"        = "#29af7f",
  "State5 × IFN"  = "#bddf26",
  "State2 × IFN"  = "#6ec5b8"
)

plot_df$effect <- factor(
  plot_df$effect,
  levels = c("State2 × IFN", "State5 × IFN", "State2", "State5"))


p_cyt <- ggplot(
  subset(plot_df, outcome == "Cytotox1"),
  aes(
    x = estimate,
    y = effect,
    fill = effect
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey40",
    size = 0.8
  ) +

  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    colour = "black",
    size = 0.7
  ) +

  geom_jitter(
    shape = 21,
    aes(fill = effect),
    colour = "black",
    stroke = 0.6,
    size = 3.8,
    alpha = 0.6,
    width = 0,
    height = 0.12
  ) +

  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 4
  ) +

  scale_fill_manual(values = sen_colors) +

  theme_classic(base_size = 14) +

  theme(
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(size = 1.1),
    axis.ticks = element_line(size = 1.1),
    axis.ticks.length = unit(0.35, "cm"),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +

  labs(
    title = "Cytotoxic niche",
    x = expression(beta~"(IFN"*gamma*" × State effect)"),
    y = NULL
  )
  
  
  p_sup <- ggplot(
  subset(plot_df, outcome == "ImmuneSupp1"),
  aes(
    x = estimate,
    y = effect,
    fill = effect
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey40",
    size = 0.8
  ) +

  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    colour = "black",
    size = 0.7
  ) +

  geom_jitter(
    shape = 21,
    aes(fill = effect),
    colour = "black",
    stroke = 0.6,
    size = 3.8,
    alpha = 0.6,
    width = 0,
    height = 0.12
  ) +

  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 4
  ) +

  scale_fill_manual(values = sen_colors) +

  theme_classic(base_size = 14) +

  theme(
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(size = 1.1),
    axis.ticks = element_line(size = 1.1),
    axis.ticks.length = unit(0.35, "cm"),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +

  labs(
    title = "Immunosuppressive niche",
    x = expression(beta~"(IFN"*gamma*" × State effect)"),
    y = NULL
  )
  
  ggsave(
  file.path(data_dir, "spatial_lag_cyto.png"),
  p_cyt,
  width = 5,
  height = 5
)

ggsave(
  file.path(data_dir, "spatial_lag_suppress.png"),
  p_sup,
  width = 5,
  height = 5
)



data_dir <- "C:/Users/mlnpa/Downloads/GSE194329_RAW"

library(dplyr)
library(tidyr)
library(ggplot2)


state_stats <- read.csv(
  file.path(data_dir, "ALL_MORAN_cell_state_RESULTS.csv")
)



df_state_moran <- state_stats %>%
  filter(feature %in% c("stup21", "stup51")) %>%
  select(sample, feature, Moran_I)



df_state_moran$state <- recode(
  df_state_moran$feature,
  "stup21" = "Senescence State 2",
  "stup51" = "Inflammatory State 5"
)

df_state_moran$state <- factor(
  df_state_moran$state,
  levels = c("Senescence State 2", "Inflammatory State 5")
)



fig4 <- ggplot(df_state_moran, aes(x = state, y = Moran_I, fill = state)) +

  geom_boxplot(
    color = "black",
    outlier.shape = NA,
    width = 0.65,
    size = 0.7
  ) +

  geom_jitter(
    shape = 21,
    colour = "black",
    stroke = 0.6,
    size = 2.7,
    alpha = 0.6,
    width = 0,
    height = 0.12
  ) +

  scale_fill_manual(values = c(
    "Senescence State 2" = "#bddf26",
    "Inflammatory State 5" = "#29af7f"
  )) +

  scale_y_continuous(
    limits = c(0, 0.85),
    breaks = seq(0, 0.8, 0.2)
  ) +

  theme_classic(base_size = 14) +

  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(size = 1.05),
    axis.ticks = element_line(size = 1.05),
    axis.ticks.length = unit(0.5, "cm"),
    legend.position = "none"
  ) +

  labs(
    title = "Spatial autocorrelation of senescence states",
    x = NULL,
    y = "Moran's I"
  )



ggsave(
  file.path(data_dir, "Fig4_State_Moran_UPDATED.png"),
  fig4,
  width = 3.2,
  height = 5
)

cat("\nFig4 Moran plot updated and saved\n")


