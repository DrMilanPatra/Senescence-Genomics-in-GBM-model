

library(glmnet)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(Hmisc)
library(tidyr)
library(reshape2)

set.seed(123)

data <- read.csv("C:/Users/mlnpa/Downloads/ssgsea_gbm.csv")
data1 <- read.csv("C:/Users/mlnpa/Downloads/idh_quanti_palbo_classif.csv")

data1_clean <- data1 %>%
  select(Sample, palbo_sig) %>%
  distinct()

data <- data %>%
  left_join(data1_clean, by = c("Sample_ID" = "Sample"))


data$OS_time  <- as.numeric(data$Overall.Survival..Months.)
data$OS_event <- ifelse(data$Overall.Survival.Status == "1:DECEASED", 1, 0)

data <- data %>%
  filter(!is.na(OS_time), !is.na(OS_event))


sen <- c(
  "Sen_senmayo","Sen_p16_top150","Sen_ras_top150",
  "Sen_IR_top150","Sen_repli_top150","Sen_CELLAGE_SEN",
  "palbo_sig"
)

ifn <- c(
  "SASP","ISG_up_AVIV","IFN_I_response",
  "ISG_nature","IFN_IFNG_gbm",
  "IFN_Hallmark_IFNg_resp","IFN_Hallmark_IFNa_resp"
)

immune <- c("Cytotoxic","Immunosuppresive")

all_programs <- c(sen, ifn, immune)


data[, all_programs] <- lapply(data[, all_programs], as.numeric)
data[, all_programs] <- scale(data[, all_programs])

data2 <- na.omit(data[, c(all_programs, "OS_time", "OS_event")])

x <- as.matrix(data2[, all_programs])
y <- Surv(data2$OS_time, data2$OS_event)


cor_matrix <- cor(x, method = "spearman")


write.csv(cor_matrix, "correlation_matrix.csv")




# =========================================================
#PAIRWISE FUNCTIONAL COUPLING
# =========================================================

pairwise_results <- data.frame()

for(s in sen){
  for(i in c(ifn, immune)){

    model <- summary(lm(data2[[i]] ~ data2[[s]]))

    pairwise_results <- rbind(
      pairwise_results,
      data.frame(
        Senescence = s,
        Axis = i,
        Beta = model$coefficients[2,1],
        P = model$coefficients[2,4],
        R2 = model$r.squared
      )
    )
  }
}

write.csv(pairwise_results, "pairwise_coupling.csv", row.names = FALSE)

# =========================================================
#DOMINANCE STRUCTURE (IMMUNE VS CYTOTOXIC)
# =========================================================

dominance <- pairwise_results %>%
  filter(Axis %in% c("Cytotoxic","Immunosuppresive")) %>%
  group_by(Senescence, Axis) %>%
  summarise(
    mean_beta = mean(Beta),
    mean_R2 = mean(R2),
    .groups = "drop"
  )

write.csv(dominance, "immune_dominance.csv", row.names = FALSE)

# =========================================================
#RIDGE COX (GLOBAL PROGNOSTIC MODEL)
# =========================================================
ridge_fit <- cv.glmnet(
  x, y,
  family = "cox",
  alpha = 0,
  nfolds = 10
)

risk_score <- predict(ridge_fit, x, s = "lambda.min", type = "link")

cindex <- concordance(Surv(OS_time, OS_event) ~ as.vector(risk_score),
                      data = data2)

print(cindex)

# =========================================================
#STABILITY SELECTION (ROBUSTNESS MAP)
# =========================================================
B <- 200
stab <- matrix(0, nrow = ncol(x), ncol = B)
rownames(stab) <- colnames(x)

for(b in 1:B){

  idx <- sample(1:nrow(x), replace = TRUE)

  fit <- cv.glmnet(x[idx,], y[idx], family = "cox", alpha = 0)

  coef_b <- coef(fit, s = "lambda.min")

  stab[,b] <- as.numeric(abs(coef_b) > 1e-5)
}

stability <- rowMeans(stab)
write.csv(stability, "stability.csv")

# =========================================================
#UNIVARIATE COX
# =========================================================
uni <- data.frame()

for(v in all_programs){

  fit <- coxph(as.formula(paste("Surv(OS_time, OS_event) ~", v)),
               data = data2)

  s <- summary(fit)

  uni <- rbind(uni, data.frame(
    Variable = v,
    HR = s$coefficients[1,"exp(coef)"],
    P = s$coefficients[1,"Pr(>|z|)"]
  ))
}

write.csv(uni, "univariate_cox.csv", row.names = FALSE)

# =========================================================
#FULL MULTIVARIATE COX
# =========================================================

multi_results <- data.frame()

model <- coxph(Surv(OS_time, OS_event) ~ ., data = data2)
s <- summary(model)

coef_table <- s$coefficients

multi_results <- data.frame(
  Variable = rownames(coef_table),
  Beta = coef_table[, "coef"],
  HR = coef_table[, "exp(coef)"],
  P = coef_table[, "Pr(>|z|)"]
)

write.csv(multi_results, "multivariable_cox_full.csv", row.names = FALSE)

# =========================================================
#R² STRUCTURE HEATMAP
# =========================================================

r2_matrix <- pairwise_results %>%
  select(Senescence, Axis, R2) %>%
  pivot_wider(names_from = Axis, values_from = R2)

write.csv(r2_matrix, "R2_matrix.csv", row.names = FALSE)


############################################
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(grid)
library(gridExtra)




# =========================
#Pairwise coupling (R2 map)
# =========================

tiff("pairwise_coupling.tiff",
     units = "in", width = 10, height = 7, res = 600, compression = "lzw")

ggplot(pairwise_results,
       aes(x = Axis, y = Senescence)) +
  geom_point(aes(size = abs(Beta), fill = R2),
             shape = 21, color = "black") +
  scale_fill_gradient(low = "white", high = "darkred") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks = element_line(color="black", size=1),
    axis.ticks.length = unit(0.3, "cm"),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid = element_blank()) +
  labs(title = "Functional coupling between senescence and immune programs",
       size = "|β|", fill = "R²")

dev.off()



# =========================
#Immunosuppressive vs Cytotoxic dominance
# =========================

dom_plot <- dominance %>%
  mutate(Type = Axis)

tiff("immune_dominance.tiff",
     units = "in", width = 8, height = 6, res = 600, compression = "lzw")

ggplot(dom_plot,
       aes(x = Senescence, y = mean_R2, fill = Axis)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Immunosuppressive bias of senescence programs",
       y = "Mean R²")

dev.off()

library(ggplot2)

tiff("univariate_cox.tiff",
     units = "in", width = 9, height = 6, res = 600, compression = "lzw")

ggplot(uni, aes(x = reorder(Variable, HR), y = HR)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = HR * 0.85, ymax = HR * 1.15), width = 0.2) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(title = "Univariate Cox hazard ratios",
       x = "",
       y = "Hazard Ratio")

dev.off()


tiff("multivariable_cox.tiff",
     units = "in", width = 8, height = 6, res = 600, compression = "lzw")

ggplot(multi_results,
       aes(x = reorder(Variable, HR), y = HR)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(title = "Multivariable Cox model (collinearity-adjusted)",
       x = "",
       y = "Hazard Ratio")

dev.off()


tiff("stability.tiff",
     units = "in", width = 8, height = 5, res = 600, compression = "lzw")

barplot(stability,
        las = 2,
        col = "steelblue",
        main = "Stability of senescence programs",
        ylab = "Selection frequency")

dev.off()

