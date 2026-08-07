## =============================================================================
## DSD SIMULATION — GRAPHICAL ANALYSIS: CORRELATION, IMPORTANCE & PREDICTIVE
## VALUE OF DESIGN FACTORS ON EVOLUTIONARY-RESCUE OUTCOMES
## =============================================================================
## Input   : dsd_simulation_design.csv  (Definitive Screening Design, 9 regimes
##           [tipo_plantilla x modo_temporal] x 27 runs, 13 numeric factors at
##           3 levels each, 2 categorical factors)
## Output  : PNG figures in ./output  +  CSV summary tables in ./output
##
## Sections
##   0. Setup & data preparation
##   1. Correlation analysis        (numeric predictors  -> responses)
##   2. Categorical association     (tipo_plantilla, modo_temporal -> responses)
##   3. Variable importance         (Random Forest %IncMSE / IncNodePurity)
##   4. Predictive value            (OOB R^2, k-fold CV R^2, predicted vs actual)
##   5. Effect / partial-dependence plots for the top predictors per response
## =============================================================================

## ---------------------------------------------------------------------------
## 0. SETUP
## ---------------------------------------------------------------------------

required_pkgs <- c("ggplot2", "reshape2", "corrplot", "randomForest",
                    "RColorBrewer", "gridExtra", "GGally")
invisible(lapply(required_pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed. Run: install.packages('%s')", p, p))
  }
}))

library(ggplot2)
library(reshape2)
library(corrplot)
library(randomForest)
library(RColorBrewer)
library(gridExtra)

set.seed(42)     # reproducibility for RF and CV folds
options(scipen = 999)  # avoid scientific notation on plot axes

## --- paths -------------------------------------------------------------
input_csv  <- "dsd_simulation_design.csv"
output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## --- read data -----------------------------------------------------------
d <- read.csv(input_csv, stringsAsFactors = FALSE)

## --- variable definitions --------------------------------------------------
responses <- c(
  "final_g_bar",
  "final_z_bar",
  "final_var_g",
  "final_N",
  "persistence_time"
)

predictors_num <- c(
  "A_max", "omega2", "theta1", "K_theta",
  "g_inicial", "N0", "p_max",
  "paso_introduccion", "mutacion_sd", "sigma_e", "D", "k", "periodo_temporal"
)

predictors_cat <- c("tipo_plantilla", "modo_temporal")

predictors <- c(predictors_num, predictors_cat)

## --- basic type handling ---------------------------------------------------
d[predictors_cat] <- lapply(d[predictors_cat], as.factor)
d$extinct <- as.logical(d$extinct)

## final_N is strongly right-skewed (0 to >20,000) -> add a log1p version used
## only for visualization; models below use the raw scale (tree methods are
## scale/monotonic-transform invariant, so this does not affect RF results).
d$log_final_N <- log1p(d$final_N)

## sanity check: keep only rows with complete predictor + response data
complete_rows <- complete.cases(d[, c(predictors, responses)])
if (any(!complete_rows)) {
  message(sprintf("Dropping %d rows with missing predictor/response values.",
                   sum(!complete_rows)))
}
d <- d[complete_rows, ]

cat(sprintf("Loaded %d simulation runs | %d numeric factors | %d categorical factors | %d responses\n",
            nrow(d), length(predictors_num), length(predictors_cat), length(responses)))

## nice display labels (used throughout for axis/legend text)
resp_labels <- c(
  final_g_bar      = "Final mean genetic value (g_bar)",
  final_z_bar      = "Final mean phenotype (z_bar)",
  final_var_g      = "Final genetic variance (Var[g])",
  final_N          = "Final population size (N)",
  persistence_time = "Persistence time"
)

## ---------------------------------------------------------------------------
## 1. CORRELATION ANALYSIS  (numeric predictors vs. responses)
## ---------------------------------------------------------------------------
## Spearman is used as the primary metric because the design uses only three
## levels per factor and several responses are skewed; Spearman captures
## monotonic (not necessarily linear) relationships robustly. Pearson is
## shown alongside for reference/comparison.

cor_pearson  <- cor(d[predictors_num], d[responses], method = "pearson",  use = "pairwise.complete.obs")
cor_spearman <- cor(d[predictors_num], d[responses], method = "spearman", use = "pairwise.complete.obs")

melt_cor <- function(mat, method_name) {
  m <- melt(mat, varnames = c("predictor", "response"), value.name = "correlation")
  m$method <- method_name
  m
}
cor_long <- rbind(melt_cor(cor_pearson, "Pearson"), melt_cor(cor_spearman, "Spearman"))
cor_long$response <- factor(cor_long$response, levels = responses,
                             labels = resp_labels[responses])
# order predictors by max |Spearman correlation| across responses (readability)
pred_order <- names(sort(apply(abs(cor_spearman), 1, max)))
cor_long$predictor <- factor(cor_long$predictor, levels = pred_order)

p_corr_heatmap <- ggplot(cor_long, aes(x = response, y = predictor, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 0, limits = c(-1, 1), name = "r") +
  facet_wrap(~method) +
  labs(title = "Correlation of numeric design factors with simulation outcomes",
       subtitle = "Pearson (linear) vs. Spearman (monotonic rank) correlation",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "01_correlation_heatmap.png"), p_corr_heatmap,
       width = 11, height = 7, dpi = 150)

write.csv(dcast(cor_long[cor_long$method == "Spearman", ], predictor ~ response, value.var = "correlation"),
          file.path(output_dir, "01_correlation_spearman_table.csv"), row.names = FALSE)

## Full predictor-predictor correlation matrix (checks collinearity in the design)
png(file.path(output_dir, "01b_predictor_collinearity.png"), width = 1000, height = 1000, res = 150)
corrplot(cor(d[predictors_num], method = "spearman"), method = "color", type = "upper",
         addCoef.col = "black", number.cex = 0.55, tl.col = "black", tl.srt = 45,
         title = "Predictor-predictor collinearity (Spearman)", mar = c(0, 0, 2, 0))
dev.off()

## ---------------------------------------------------------------------------
## 2. CATEGORICAL ASSOCIATION  (tipo_plantilla, modo_temporal vs. responses)
## ---------------------------------------------------------------------------
## Association strength quantified with eta-squared (proportion of variance
## in the response explained by group membership, one-way ANOVA), which is
## the categorical analogue of the correlation coefficients above.

eta_squared <- function(y, g) {
  fit <- aov(y ~ g)
  ss <- summary(fit)[[1]][["Sum Sq"]]
  ss[1] / sum(ss)
}

eta_tab <- expand.grid(predictor = predictors_cat, response = responses,
                        stringsAsFactors = FALSE)
eta_tab$eta_squared <- mapply(function(p, r) eta_squared(d[[r]], d[[p]]),
                               eta_tab$predictor, eta_tab$response)
eta_tab$response <- factor(eta_tab$response, levels = responses, labels = resp_labels[responses])

p_eta <- ggplot(eta_tab, aes(x = response, y = eta_squared, fill = predictor)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", eta_squared)),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c(tipo_plantilla = "#4477AA", modo_temporal = "#EE6677"),
                     name = "Categorical factor") +
  labs(title = expression(paste("Association of categorical factors with outcomes (", eta^2, ")")),
       subtitle = "Share of response variance explained by group membership (one-way ANOVA)",
       x = NULL, y = expression(eta^2)) +
  ylim(0, max(eta_tab$eta_squared) * 1.25) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "02_categorical_eta_squared.png"), p_eta, width = 9, height = 6, dpi = 150)

## Boxplots: response distribution by category, for every response x categorical pair
d_box <- melt(d[, c(predictors_cat, responses)], id.vars = predictors_cat,
              variable.name = "response", value.name = "value")
d_box$response <- factor(d_box$response, levels = responses, labels = resp_labels[responses])

p_box_tipo <- ggplot(d_box, aes(x = tipo_plantilla, y = value, fill = tipo_plantilla)) +
  geom_boxplot(alpha = 0.85, outlier.size = 0.8) +
  facet_wrap(~response, scales = "free_y") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Outcome distributions by spatial template (tipo_plantilla)", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "02b_boxplots_tipo_plantilla.png"), p_box_tipo, width = 10, height = 7, dpi = 150)

p_box_modo <- ggplot(d_box, aes(x = modo_temporal, y = value, fill = modo_temporal)) +
  geom_boxplot(alpha = 0.85, outlier.size = 0.8) +
  facet_wrap(~response, scales = "free_y") +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  labs(title = "Outcome distributions by temporal regime (modo_temporal)", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "02c_boxplots_modo_temporal.png"), p_box_modo, width = 10, height = 7, dpi = 150)

## ---------------------------------------------------------------------------
## 3. VARIABLE IMPORTANCE  (Random Forest, one model per response)
## ---------------------------------------------------------------------------
## Random Forest importance captures nonlinear effects and interactions that
## the correlation/eta-squared summaries above cannot, and works natively
## with the mix of numeric and categorical predictors.

rf_models <- list()
importance_list <- list()

for (r in responses) {
  form <- as.formula(paste(r, "~", paste(predictors, collapse = " + ")))
  rf <- randomForest(form, data = d[, c(predictors, r)],
                      ntree = 1000, importance = TRUE, mtry = max(1, floor(length(predictors) / 3)))
  rf_models[[r]] <- rf

  imp <- as.data.frame(importance(rf))
  imp$predictor <- rownames(imp)
  imp$response  <- r
  importance_list[[r]] <- imp
}

imp_all <- do.call(rbind, importance_list)
rownames(imp_all) <- NULL
imp_all$response <- factor(imp_all$response, levels = responses, labels = resp_labels[responses])

# order predictors *within each response facet* by %IncMSE (manual "reorder
# within group", since ggplot's reorder() only supports one global order)
imp_all <- imp_all[order(imp_all$response, imp_all$`%IncMSE`), ]
imp_all$predictor_ord <- factor(paste(imp_all$response, imp_all$predictor, sep = "___"),
                                 levels = paste(imp_all$response, imp_all$predictor, sep = "___"))

p_importance <- ggplot(imp_all, aes(x = predictor_ord, y = `%IncMSE`)) +
  geom_col(aes(fill = `%IncMSE` > 0), width = 0.7) +
  coord_flip() +
  scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
  facet_wrap(~response, scales = "free", ncol = 2) +
  scale_fill_manual(values = c(`TRUE` = "#2E86AB", `FALSE` = "grey70"), guide = "none") +
  labs(title = "Random Forest variable importance by response",
       subtitle = "%IncMSE: increase in mean squared error when a predictor is permuted (higher = more important)",
       x = NULL, y = "% Increase in MSE") +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "03_variable_importance.png"), p_importance, width = 11, height = 10, dpi = 150)

imp_out <- imp_all[order(imp_all$response, -imp_all$`%IncMSE`), setdiff(names(imp_all), "predictor_ord")]
write.csv(imp_out, file.path(output_dir, "03_variable_importance_table.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 4. PREDICTIVE VALUE  (OOB R^2, 5-fold CV R^2, predicted vs. actual)
## ---------------------------------------------------------------------------

r_squared <- function(actual, predicted) {
  1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
}

k_fold_cv_r2 <- function(response, k = 5) {
  n <- nrow(d)
  folds <- sample(rep(1:k, length.out = n))
  preds <- numeric(n)
  for (i in 1:k) {
    train <- d[folds != i, c(predictors, response)]
    test  <- d[folds == i, c(predictors, response)]
    form <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
    rf_i <- randomForest(form, data = train, ntree = 500)
    preds[folds == i] <- predict(rf_i, newdata = test)
  }
  r_squared(d[[response]], preds)
}

perf_list <- list()
pred_vs_actual_list <- list()

for (r in responses) {
  rf <- rf_models[[r]]
  oob_pred <- predict(rf)                      # out-of-bag predictions
  oob_r2   <- r_squared(d[[r]], oob_pred)
  cv_r2    <- k_fold_cv_r2(r, k = 5)

  perf_list[[r]] <- data.frame(response = r, metric = c("OOB R2", "5-fold CV R2"),
                                value = c(oob_r2, cv_r2))
  pred_vs_actual_list[[r]] <- data.frame(response = r, actual = d[[r]], predicted = oob_pred)
}

perf_all <- do.call(rbind, perf_list)
perf_all$response <- factor(perf_all$response, levels = responses, labels = resp_labels[responses])

p_perf <- ggplot(perf_all, aes(x = response, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", value)),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("OOB R2" = "#2E86AB", "5-fold CV R2" = "#F18F01"),
                     name = NULL) +
  labs(title = "Predictive value of the design factors (Random Forest)",
       subtitle = "Out-of-bag R\u00b2 vs. independent 5-fold cross-validated R\u00b2",
       x = NULL, y = expression(R^2)) +
  ylim(min(0, min(perf_all$value) * 1.2), max(perf_all$value) * 1.2) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "04_predictive_r2.png"), p_perf, width = 9, height = 6, dpi = 150)
write.csv(perf_all, file.path(output_dir, "04_predictive_r2_table.csv"), row.names = FALSE)

pva_all <- do.call(rbind, pred_vs_actual_list)
pva_all$response <- factor(pva_all$response, levels = responses, labels = resp_labels[responses])

p_pva <- ggplot(pva_all, aes(x = actual, y = predicted)) +
  geom_point(alpha = 0.55, color = "#2E86AB", size = 1.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_smooth(method = "lm", se = FALSE, color = "#B2182B", linewidth = 0.6) +
  facet_wrap(~response, scales = "free", ncol = 2) +
  labs(title = "Out-of-bag predicted vs. actual values",
       subtitle = "Dashed line = perfect prediction; red line = fitted trend",
       x = "Actual", y = "OOB predicted") +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "04b_predicted_vs_actual.png"), p_pva, width = 11, height = 9, dpi = 150)

## ---------------------------------------------------------------------------
## 5. EFFECT / PARTIAL-DEPENDENCE PLOTS for the top predictors per response
## ---------------------------------------------------------------------------
## Shows the marginal (average) effect of the single most important numeric
## predictor for each response, holding all other predictors at their
## observed distribution (Breiman's partial dependence, via randomForest::partialPlot).

top_n_effects <- 3  # top predictors per response to visualize

pdp_list <- list()
for (r in responses) {
  imp_r <- imp_all[imp_all$response == resp_labels[[r]], ]
  imp_r <- imp_r[imp_r$predictor %in% predictors_num, ]  # partialPlot below handles numeric cleanly
  top_preds <- imp_r$predictor[order(-imp_r$`%IncMSE`)][1:min(top_n_effects, nrow(imp_r))]

  for (p in top_preds) {
    pp <- do.call(partialPlot, list(x = rf_models[[r]], pred.data = d[, c(predictors, r)],
                                     x.var = p, plot = FALSE))
    pdp_list[[length(pdp_list) + 1]] <- data.frame(
      response = r, predictor = p, x = pp$x, y = pp$y
    )
  }
}
pdp_all <- do.call(rbind, pdp_list)
pdp_all$response <- factor(pdp_all$response, levels = responses, labels = resp_labels[responses])
pdp_all$panel <- paste0(pdp_all$response, "\n(", pdp_all$predictor, ")")

p_pdp <- ggplot(pdp_all, aes(x = x, y = y)) +
  geom_line(color = "#2E86AB", linewidth = 0.9) +
  geom_point(color = "#2E86AB", size = 1.3) +
  facet_wrap(~panel, scales = "free", ncol = 3) +
  labs(title = "Partial dependence: top predictors per response",
       subtitle = paste0("Marginal effect of each factor, averaging over all other predictors (top ",
                          top_n_effects, " per response)"),
       x = "Predictor value", y = "Predicted response (partial effect)") +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "05_partial_dependence.png"), p_pdp, width = 13, height = 12, dpi = 150)

## ---------------------------------------------------------------------------
## DONE
## ---------------------------------------------------------------------------
cat("\nAnalysis complete. Figures and tables written to:", normalizePath(output_dir), "\n")
cat(list.files(output_dir), sep = "\n")

