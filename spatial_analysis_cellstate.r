#Cellstate spatial genomics analysis

############################################################
#Train model using ONLY shared features
############################################################

rm(list = ls())

library(Seurat)
library(dplyr)
library(glmnet)
library(FNN)
library(spdep)
library(spatialreg)
library(ggplot2)
library(patchwork)

load("C:/Users/mlnpa/Documents/glioblastoma_scrnaseq/gbm_GSE182109_filtered_scores_heterogeneity.RData")

shared_features <- c(
  "IFN_I_response.score1",
  "CELLAGE_SEN.score1",
  "repli_top150.score1",
  "ras_top150.score1",
  "IR_top150.score1",
  "SASP.score1",
  "ISG_nature.score1",
  "Hallmark_IFNg_resp.score1",
  "Hallmark_IFNa_resp.score1",
  "p16_top150.score1",
  "senmayo.score1"
)

meta_sc <- se@meta.data
shared_features <- intersect(shared_features, colnames(meta_sc))

X_sc <- as.matrix(meta_sc[, shared_features])
y_sc <- as.factor(meta_sc$CellState)

############################################################
# IMPORTANT: store scaling parameters
############################################################

sc_center <- colMeans(X_sc, na.rm = TRUE)
sc_scale  <- apply(X_sc, 2, sd, na.rm = TRUE)

X_sc_scaled <- scale(X_sc, center = sc_center, scale = sc_scale)

fit_state <- cv.glmnet(
  X_sc_scaled,
  y_sc,
  family = "multinomial",
  alpha = 0.5
)

saveRDS(fit_state, "fit_state_model_shared.rds")

###############################################################################

############################################################
# CELL STATE SPATIAL GENOMICS ANALYSIS (FIXED PIPELINE)
############################################################

rm(list = ls())

library(Seurat)
library(dplyr)
library(glmnet)
library(FNN)
library(spdep)
library(spatialreg)
library(ggplot2)
library(patchwork)

############################################################
# 1. TRAIN MODEL (scRNA-seq)
############################################################

load("C:/Users/mlnpa/Documents/glioblastoma_scrnaseq/gbm_GSE182109_filtered_scores_heterogeneity.RData")

shared_features <- c(
  "IFN_I_response.score1",
  "CELLAGE_SEN.score1",
  "repli_top150.score1",
  "ras_top150.score1",
  "IR_top150.score1",
  "SASP.score1",
  "ISG_nature.score1",
  "Hallmark_IFNg_resp.score1",
  "Hallmark_IFNa_resp.score1",
  "p16_top150.score1",
  "senmayo.score1"
)

meta_sc <- se@meta.data
shared_features <- intersect(shared_features, colnames(meta_sc))

X_sc <- as.matrix(meta_sc[, shared_features])
y_sc <- as.factor(meta_sc$CellState)

############################################################
# CRITICAL FIX: store TRAINING scaling parameters
############################################################

sc_center <- colMeans(X_sc, na.rm = TRUE)
sc_scale  <- apply(X_sc, 2, sd, na.rm = TRUE)

X_sc_scaled <- scale(X_sc, center = sc_center, scale = sc_scale)

fit_state <- cv.glmnet(
  X_sc_scaled,
  y_sc,
  family = "multinomial",
  alpha = 0.5
)

saveRDS(
  list(model = fit_state,
       center = sc_center,
       scale = sc_scale),
  "fit_state_model_shared.rds"
)

############################################################
# 2. SPATIAL PROJECTION
############################################################

rm(list = ls())

library(Seurat)
library(dplyr)
library(glmnet)
library(FNN)
library(spdep)
library(spatialreg)
library(ggplot2)
library(patchwork)

data_dir <- "C:/Users/mlnpa/Downloads/GSE194329_RAW"

samples <- c(
  "GBM1_spaceranger_out",
  "GBM2_spaceranger_out",
  "GBM3_spaceranger_out",
  "GBM4_spacer_out",
  "GBM5_1_spaceranger_out",
  "GBM5_2_spaceranger_out"
)

spatial_features <- c(
  "IFN_I_response1",
  "Sen_CELLAGE_SEN1",
  "Sen_repli_top1501",
  "Sen_ras_top1501",
  "Sen_IR_top1501",
  "SASP1",
  "ISG_nature1",
  "IFN_Hallmark_IFNg_resp1",
  "IFN_Hallmark_IFNa_resp1",
  "Sen_p16_top1501",
  "Sen_senmayo1"
)

############################################################
# LOAD MODEL + SCALING
############################################################

obj <- readRDS("fit_state_model_shared.rds")
fit_state <- obj$model
sc_center <- obj$center
sc_scale  <- obj$scale

############################################################
# HELPERS
############################################################

get_p <- function(model, var){
  coefs <- summary(model)$coefficients
  if(var %in% rownames(coefs)) coefs[var,4] else NA
}

get_spatial_coef <- function(model, var){
  coefs <- summary(model)$Coef
  if(var %in% rownames(coefs)){
    c(beta = coefs[var,1], p = coefs[var,4])
  } else {
    c(beta = NA, p = NA)
  }
}

state_results <- data.frame()

############################################################
# 3. MAIN LOOP
############################################################

for(s in samples){

  cat("\nProcessing:", s, "\n")

  visium <- Load10X_Spatial(data.dir = file.path(data_dir, s))
  visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE)

  df <- read.csv(file.path(data_dir, paste0(s, "_spot_results.csv")))

  ##########################################################
  # FEATURE MATRIX (strict alignment)
  ##########################################################

  X_spatial <- matrix(0, nrow = nrow(df), ncol = length(spatial_features))
  colnames(X_spatial) <- spatial_features

  for(f in spatial_features){
    if(f %in% colnames(df)){
      X_spatial[, f] <- df[[f]]
    }
  }

  ##########################################################
  # CRITICAL FIX: APPLY SAME SCALING AS TRAINING
  ##########################################################

  X_spatial_scaled <- scale(
    X_spatial,
    center = sc_center,
    scale  = sc_scale
  )

  ##########################################################
  # PREDICTION (NOW CORRECTLY SCALED)
  ##########################################################

  state_prob <- predict(
    fit_state,
    newx = X_spatial_scaled,
    s = "lambda.min",
    type = "response"
  )

  state_prob_mat <- state_prob[, , 1]

  df$State2_prob <- state_prob_mat[, "State2"]
  df$State5_prob <- state_prob_mat[, "State5"]

  visium$State2_prob <- df$State2_prob
  visium$State5_prob <- df$State5_prob

  ##########################################################
  # CLASS PREDICTION
  ##########################################################

  state_class <- predict(
    fit_state,
    newx = X_spatial_scaled,
    s = "lambda.min",
    type = "class"
  )

  df$PredictedState <- as.factor(state_class[,1])
  visium$PredictedState <- df$PredictedState

  ##########################################################
  # SPATIAL MAPS (SAVE)
  ##########################################################

  p1 <- SpatialFeaturePlot(visium, "State2_prob")
  p2 <- SpatialFeaturePlot(visium, "State5_prob")
  p3 <- SpatialDimPlot(visium, group.by = "PredictedState")

  ggsave(
    file.path(data_dir, paste0(s, "_STATE_SPATIAL.png")),
    p1 + p2 + p3,
    width = 12,
    height = 4
  )

  ##########################################################
  # SPATIAL GRAPH
  ##########################################################

  coords <- as.matrix(df[, c("imagecol", "imagerow")])
  nb <- knn2nb(knearneigh(coords, k = 6))
  lw <- nb2listw(nb)

  ##########################################################
  # MORAN
  ##########################################################

  moran_s2 <- moran.test(df$State2_prob, lw)
  moran_s5 <- moran.test(df$State5_prob, lw)

  ##########################################################
  # REGRESSION
  ##########################################################

  fit_cyt <- lm(Cytotox1 ~ State2_prob + State5_prob, data = df)
  fit_sup <- lm(ImmuneSupp1 ~ State2_prob + State5_prob, data = df)

  ##########################################################
  # SPATIAL LAG MODEL
  ##########################################################

  fit_spatial <- lagsarlm(
    ImmuneSupp1 ~ State5_prob + SASP1,
    data = df,
    listw = lw
  )

  sp_res <- get_spatial_coef(fit_spatial, "State5_prob")

  ##########################################################
  # STORE RESULTS
  ##########################################################

  state_results <- rbind(state_results, data.frame(
    sample = s,

    Moran_State2 = moran_s2$estimate[1],
    Moran_State5 = moran_s5$estimate[1],

    Moran_State2_p = moran_s2$p.value,
    Moran_State5_p = moran_s5$p.value,

    Cyt_State2_p = get_p(fit_cyt, "State2_prob"),
    Cyt_State5_p = get_p(fit_cyt, "State5_prob"),

    Supp_State2_p = get_p(fit_sup, "State2_prob"),
    Supp_State5_p = get_p(fit_sup, "State5_prob"),

    Spatial_beta = sp_res["beta"],
    Spatial_p = sp_res["p"]
  ))

  write.csv(
    df,
    file.path(data_dir, paste0(s, "_STATE_ENRICHED.csv")),
    row.names = FALSE
  )
}

############################################################
# FINAL OUTPUT
############################################################

write.csv(
  state_results,
  file.path(data_dir, "STATE_SPATIAL_RESULTS_UPDATED.csv"),
  row.names = FALSE
)

cat("\nState integration analysis complete.\n")

