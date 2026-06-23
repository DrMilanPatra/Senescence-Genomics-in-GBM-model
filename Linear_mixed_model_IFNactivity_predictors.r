## Mixed model for IFN response and senescence programs

rm(list = ls())


library(Seurat)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(gplots)
library(ggsci)
library(ggplot2)
library(gridExtra)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(performance)
library(grid)


load(file = "C:/Users/mlnpa/Documents/glioblastoma_scrnaseq/gbm_GSE182109_filtered_scores.RData")


df <- se@meta.data[, c(
  "Patient",
  "Hallmark_IFNg_resp.score1",
  "Hallmark_IFNa_resp.score1",
  "ISG_nature.score1",
  "p16_top150.score1",
  "IR_top150.score1",
  "repli_top150.score1",
  "ras_top150.score1",
  "senmayo.score1"
)]

df <- na.omit(df)

## Scale all numeric variables
df_scaled <- df
df_scaled[, -1] <- scale(df_scaled[, -1])

## Patient as factor
df_scaled$Patient <- as.factor(df_scaled$Patient)



ifn_programs <- c(
  "Hallmark_IFNg_resp.score1",
  "Hallmark_IFNa_resp.score1",
  "ISG_nature.score1"
)

##############################################################################
## Fit mixed models
##############################################################################

models <- lapply(ifn_programs, function(x) {

  formula <- as.formula(
    paste(
      x,
      "~ p16_top150.score1 +",
      "IR_top150.score1 +",
      "repli_top150.score1 +",
      "ras_top150.score1 +",
      "senmayo.score1 +",
      "(1|Patient)"
    )
  )

  lmer(formula, data = df_scaled)

})

names(models) <- ifn_programs



summary(models[[1]])
summary(models[[2]])
summary(models[[3]])

##############################################################################
## Extract coefficients
##############################################################################

results <- lapply(models, broom.mixed::tidy)

results_df <- bind_rows(results, .id = "IFN_program")

results_fixed <- results_df %>%
  filter(effect == "fixed", term != "(Intercept)")

print(results_fixed)

##############################################################################
## Extract p-values
##############################################################################

pval_table <- results_fixed %>%
  select(
    IFN_program,
    term,
    estimate,
    std.error,
    statistic,
    p.value
  )

print(pval_table)

write.csv(
  pval_table,
  "mixed_model_coefficients_pvalues.csv",
  row.names = FALSE
)

##############################################################################
## Calculate marginal and conditional R²
## Using Nakagawa method
##############################################################################

r2_results <- lapply(models, performance::r2_nakagawa)

r2_df <- lapply(r2_results, function(x) {

  data.frame(
    R2m = x$R2_marginal,
    R2c = x$R2_conditional
  )

}) %>%
  bind_rows(.id = "IFN_program")

print(r2_df)

write.csv(
  r2_df,
  "mixed_model_R2_values.csv",
  row.names = FALSE
)

##############################################################################
## Collinearity check
##############################################################################

collinearity_results <- lapply(models, check_collinearity)

print(collinearity_results)

##############################################################################
## Plotting data
##############################################################################

## Reorder senescence programs by mean absolute effect size
term_order <- results_fixed %>%
  group_by(term) %>%
  summarize(mean_abs_est = mean(abs(estimate))) %>%
  arrange(desc(mean_abs_est)) %>%
  pull(term)

plot_data <- results_fixed %>%
  mutate(
    IFN_program = factor(
      IFN_program,
      levels = c(
        "Hallmark_IFNg_resp.score1",
        "Hallmark_IFNa_resp.score1",
        "ISG_nature.score1"
      )
    ),

    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error,

    term = factor(
      term,
      levels = rev(term_order)
    )
  )



r2_labels <- r2_df %>%
  mutate(
    label = paste0(
      "R²m = ",
      round(R2m, 2),
      "\nR²c = ",
      round(R2c, 2)
    )
  )



fill_colors <- c(
  "#fdf2b8",
  "#e88200",
  "#cb2800"
)



tiff(
  "forestplot.tiff",
  units = "in",
  width = 12,
  height = 5,
  res = 600,
  compression = "lzw"
)

ggplot(
  plot_data,
  aes(
    x = estimate,
    y = term,
    fill = IFN_program
  )
) +

  geom_point(
    shape = 21,
    position = position_dodge(width = 0.6),
    size = 4,
    stroke = 0.8,
    color = "black"
  ) +

  geom_errorbarh(
    aes(
      xmin = lower,
      xmax = upper
    ),
    position = position_dodge(width = 0.6),
    height = 0.3,
    color = "black",
    size = 0.8
  ) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +

  labs(
    x = "Effect size (β)",
    y = "Senescence program",
    fill = "IFN program"
  ) +

  scale_fill_manual(values = fill_colors) +

  theme_classic(base_size = 12) +

  theme(
    axis.text.y = element_text(
      face = "bold",
      color = "black",
      size = 14
    ),

    axis.text.x = element_text(
      color = "black",
      size = 14
    ),

    axis.line = element_line(
      color = "black",
      size = 1
    ),

    axis.ticks = element_line(
      color = "black",
      size = 1
    ),

    axis.ticks.length = unit(0.35, "cm"),

    legend.position = "right"
  ) +

  ## Add R² labels
  geom_text(
    data = r2_labels,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.5,
    size = 4
  )

dev.off()

##############################################################################
## Save plotting table
##############################################################################

write.csv(
  plot_data,
  "mixed_model_forestplot_data.csv",
  row.names = FALSE
)


