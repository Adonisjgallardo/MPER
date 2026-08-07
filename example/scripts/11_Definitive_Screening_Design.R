## ============================================================
## 11_Definitive_Screening_Design.R
## Definitive Screening Design (DSD) with log-scaled parameters
## Generates a full experimental matrix across 9 spatio-temporal regimes
## ============================================================

# ------------------------------------------------------------------
# 0. Setup environment and locate the project root
# ------------------------------------------------------------------
find_project_root <- function(start_dir = getwd()) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  while (!is.na(current_dir) && nzchar(current_dir)) {
    if (file.exists(file.path(current_dir, "scripts", "05_motor.R"))) {
      return(current_dir)
    }
    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      break
    }
    current_dir <- parent_dir
  }
  return(getwd())
}

run_definitive_screening_design <- function() {
  project_root <- find_project_root()
  source(file.path(project_root, "scripts", "05_motor.R"))
  simulation_hca <- if (exists("simular_hca", mode = "function")) simular_hca else NULL

  if (requireNamespace("DoE.wrapper", quietly = TRUE)) {
    suppressPackageStartupMessages(library(DoE.wrapper))
  } else {
    message("DoE.wrapper is not available; using a built-in 3-level screening fallback design.")
  }

  if (requireNamespace("tidyverse", quietly = TRUE)) {
    suppressPackageStartupMessages(library(tidyverse))
  } else {
    message("tidyverse is not available; using base R operations for the design workflow.")
  }

    if (requireNamespace("openxlsx", quietly = TRUE)) {
    suppressPackageStartupMessages(library(openxlsx))
  } else {
    message("openxlsx is not available; writing CSV/RDS output only.")
  }

  set.seed(42)

  # Create output directories for this experiment (project-root relative)
  exp_data_dir <- file.path(project_root, "example", "data", "experiments", "definitive_screening")
  exp_plots_dir <- file.path(project_root, "example", "plots", "experiments", "definitive_screening")
  dir.create(exp_data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(exp_plots_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# 1. Define factor levels and scale types
# ------------------------------------------------------------------
# Linear factors use the three provided levels directly.
# Log-scaled factors use the low/high range and the center is computed
# as the geometric mean in physical units.
factors_config <- list(
  A_max = list(type = "linear", levels = c(0.15, 0.55, 0.99)),
  omega2 = list(type = "linear", levels = c(0.01, 0.05, 0.09)),
  theta1 = list(type = "linear", levels = c(0.40, 0.80, 1.20)),
  K_theta = list(type = "linear", levels = c(0.10, 0.35, 0.60)),
  g_inicial = list(type = "linear", levels = c(0.02, 0.16, 0.30)),
  N0 = list(type = "linear", levels = c(0.30, 1.15, 2.00)),
  p_max = list(type = "linear", levels = c(0.10, 0.50, 0.90)),
  paso_introduccion = list(type = "linear", levels = c(100, 250, 400)),

  mutacion_sd = list(type = "log10", levels = c(0.001, 0.10)),
  sigma_e = list(type = "log10", levels = c(0.001, 0.10)),
  D = list(type = "log10", levels = c(0.01, 0.49)),
  k = list(type = "log10", levels = c(0.01, 0.30)),
  periodo_temporal = list(type = "log10", levels = c(5, 20))
)

spatio_temporal_grid <- expand.grid(
  tipo_plantilla = c("uniforme", "lineal", "radial"),
  modo_temporal = c("constante", "pulsos", "sinusoidal"),
  stringsAsFactors = FALSE
)

continuous_factors <- names(factors_config)

# ------------------------------------------------------------------
# 2. Resolve low/center/high levels in physical space
# ------------------------------------------------------------------
factor_levels_table <- do.call(rbind, lapply(continuous_factors, function(f_name) {
  cfg <- factors_config[[f_name]]

  if (cfg$type == "linear") {
    low <- cfg$levels[1]
    center <- cfg$levels[2]
    high <- cfg$levels[3]
  } else if (cfg$type == "log10") {
    low <- cfg$levels[1]
    high <- cfg$levels[2]
    center <- 10^((log10(low) + log10(high)) / 2)
  }

  data.frame(
    Factor = f_name,
    Scale = cfg$type,
    Low = low,
    Center = round(center, 4),
    High = high,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}))

print("Resolved factor levels (including log-geometric centers):")
print(factor_levels_table)

# ------------------------------------------------------------------
# 3. Generate the base DSD matrix in coded and decoded form
# ------------------------------------------------------------------
# Generate normalized 3-level screening matrix (-1, 0, +1).
# If DoE.wrapper::dsd is available we use it; otherwise we fall back to a
# compact 3-level matrix that still covers the center and both extremes.
generate_dsd_matrix <- function(nfactors) {
  if (!is.numeric(nfactors) || length(nfactors) != 1 || nfactors < 1) {
    stop("nfactors must be a single positive integer")
  }

  n_runs <- 2 * nfactors + 1L
  design <- matrix(0, nrow = n_runs, ncol = nfactors)
  for (j in seq_len(nfactors)) {
    design[2 * j - 1, j] <- 1
    design[2 * j, j] <- -1
  }
  design
}

if (exists("dsd", mode = "function") && requireNamespace("DoE.wrapper", quietly = TRUE)) {
  dsd_base_coded <- as.data.frame(dsd(nfactors = length(continuous_factors)))
} else {
  dsd_base_coded <- as.data.frame(generate_dsd_matrix(length(continuous_factors)))
}
colnames(dsd_base_coded) <- continuous_factors

# Decode from [-1, 0, +1] to actual simulation values.
# The mapping is linear for regular factors and linear in log10-space
# for log-scaled variables.
dsd_base_natural <- as.data.frame(dsd_base_coded)

for (f_name in continuous_factors) {
  cfg <- factors_config[[f_name]]
  coded_vals <- dsd_base_coded[[f_name]]

  if (cfg$type == "linear") {
    low <- cfg$levels[1]
    center <- cfg$levels[2]
    high <- cfg$levels[3]

    dsd_base_natural[[f_name]] <- ifelse(
      coded_vals == -1, low,
      ifelse(coded_vals == 0, center, high)
    )
  } else if (cfg$type == "log10") {
    log_low <- log10(cfg$levels[1])
    log_high <- log10(cfg$levels[2])
    log_center <- (log_low + log_high) / 2

    dsd_base_natural[[f_name]] <- ifelse(
      coded_vals == -1, 10^log_low,
      ifelse(coded_vals == 0, 10^log_center, 10^log_high)
    )
  }
}

dsd_base_natural$dsd_run_id <- seq_len(nrow(dsd_base_natural))

# ------------------------------------------------------------------
# 4. Expand the design across spatio-temporal regimes
# ------------------------------------------------------------------
full_simulation_design_list <- lapply(seq_len(nrow(spatio_temporal_grid)), function(i) {
  regime_row <- spatio_temporal_grid[i, , drop = FALSE]
  regime_row_df <- as.data.frame(regime_row, stringsAsFactors = FALSE, row.names = NULL)
  out <- cbind(
    regime_id = i,
    regime_row_df,
    dsd_base_natural,
    stringsAsFactors = FALSE
  )
  out$sim_id <- seq_len(nrow(out))
  out
})

full_simulation_design <- do.call(rbind, full_simulation_design_list)
rownames(full_simulation_design) <- NULL
full_simulation_design <- full_simulation_design[, c(
  "sim_id", "regime_id", "dsd_run_id", "tipo_plantilla", "modo_temporal", continuous_factors
)]

# ------------------------------------------------------------------
# 5. Simulate each design row and attach response metrics
# ------------------------------------------------------------------
if (is.null(simulation_hca)) {
  stop("simulation_hca() is not available. Make sure scripts/05_motor.R defines simular_hca() and is sourced.")
}

extract_simulation_responses <- function(sim) {
  last_row <- tail(sim$resumen, 1)
  persistence_time <- if (is.na(sim$paso_extincion)) {
    max(sim$resumen$paso, na.rm = TRUE) - sim$paso_introduccion
  } else {
    sim$paso_extincion - sim$paso_introduccion
  }

  data.frame(
    final_g_bar = last_row$g_bar,
    final_z_bar = last_row$z_bar,
    final_var_g = last_row$var_g,
    final_N = last_row$n_ocupados,
    paso_extincion = sim$paso_extincion,
    persistence_time = persistence_time,
    extinct = !is.na(sim$paso_extincion),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

run_simulation_row <- function(row) {
  params <- list(
    nx = 150,
    ny = 150,
    n_steps = 500,
    dt = 0.25,
    N_umbral = 0.15,
    N_half = 0.3,
    D_A = 0.12,
    delta_A = 0.02,
    tasa_dosificacion = 0.15,
    mort_base = 0.01,
    mort_estres = 0.9,
    K_mort = 0.3,
    Nc = 200,
    guardar_cada = 25,
    seed = 1000 + as.integer(row$sim_id)
  )

  params <- c(params,
              as.list(row[c("tipo_plantilla", "modo_temporal", continuous_factors)]),
              list(theta0 = 0.1))
  sim <- do.call(simulation_hca, params)
  extract_simulation_responses(sim)
}

response_rows <- lapply(seq_len(nrow(full_simulation_design)), function(i) {
  run_simulation_row(full_simulation_design[i, , drop = FALSE])
})

response_df <- do.call(rbind, response_rows)
full_simulation_design <- cbind(full_simulation_design, response_df)

# ------------------------------------------------------------------
# 6. Export the final design matrix
# ------------------------------------------------------------------
output_csv <- file.path(exp_data_dir, "dsd_simulation_design_log_scaled.csv")
output_xlsx <- file.path(exp_data_dir, "dsd_simulation_design_log_scaled.xlsx")
output_rds <- file.path(exp_data_dir, "dsd_simulation_design_log_scaled.rds")

write.csv(full_simulation_design, output_csv, row.names = FALSE)
saveRDS(full_simulation_design, output_rds)

if (file.exists(output_xlsx)) {
  file.remove(output_xlsx)
}
if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(full_simulation_design, output_xlsx, rowNames = FALSE)
  cat(" -", normalizePath(output_xlsx), "\n")
} else {
  message("openxlsx is not installed; skipped Excel export.")
}

cat("Saved design with", nrow(full_simulation_design), "runs to:\n")
cat(" -", normalizePath(output_csv), "\n")
cat(" -", normalizePath(output_rds), "\n")

# ------------------------------------------------------------------
# 6. Regression-analysis template for later modelling
# ------------------------------------------------------------------
analyze_dsd_results <- function(
    results_df,
    responses = c("final_g_bar", "final_z_bar", "final_var_g", "final_N", "persistence_time"),
    categorical = c("tipo_plantilla", "modo_temporal"),
    output_data_dir = NULL,
    output_plots_dir = NULL,
    run_sobol = FALSE,
    sobol_n = 500,
    shap_nsim = 50,
    install_missing = FALSE
) {
  # ------------------------------------------------------------------
  # 0. Setup and package checks
  # ------------------------------------------------------------------
  required_pkgs <- c(
    "tidyverse", "car", "mgcv", "randomForest", "pdp",
    "fastshap", "shapviz", "DiceKriging", "sensitivity", "lhs",
    "caret", "ggplot2", "gridExtra"
  )

  if (install_missing) {
    for (pkg in required_pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
      }
    }
  }

  # Load packages quietly
  suppressPackageStartupMessages({
    library(tidyverse)
    library(car)
    library(mgcv)
    library(randomForest)
    library(pdp)
    library(fastshap)
    library(shapviz)
    library(DiceKriging)
    library(sensitivity)
    library(lhs)
    library(caret)
    library(ggplot2)
    library(gridExtra)
  })

  # Detect log‑scaled factors from global factors_config if available
  if (exists("factors_config", envir = .GlobalEnv)) {
    log_factors <- names(factors_config)[sapply(factors_config, function(x) x$type == "log10")]
  } else {
    log_factors <- NULL
    warning("'factors_config' not found; no log‑transformations will be applied.")
  }

  # If output directories not provided, use global variables
  if (is.null(output_data_dir)) {
    if (exists("exp_data_dir", envir = .GlobalEnv)) {
      output_data_dir <- exp_data_dir
    } else {
      output_data_dir <- file.path(getwd(), "example", "data", "experiments", "definitive_screening")
    }
  }
  if (is.null(output_plots_dir)) {
    if (exists("exp_plots_dir", envir = .GlobalEnv)) {
      output_plots_dir <- exp_plots_dir
    } else {
      output_plots_dir <- file.path(getwd(), "example", "plots", "experiments", "definitive_screening")
    }
  }

  # Create response-specific subdirectories
  for (resp in responses) {
    dir.create(file.path(output_data_dir, resp), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_plots_dir, resp), recursive = TRUE, showWarnings = FALSE)
  }

  # ------------------------------------------------------------------
  # 1. Prepare data
  # ------------------------------------------------------------------
  model_df <- results_df %>%
    mutate(across(all_of(categorical), as.factor))

  # Add log‑transformed columns for log‑scaled factors
  for (f in log_factors) {
    if (f %in% names(model_df)) {
      model_df[[paste0("log_", f)]] <- log10(model_df[[f]])
    }
  }

  # Define numeric predictors (original + log‑transformed)
  numeric_predictors <- setdiff(names(model_df), c(responses, categorical, "sim_id", "regime_id", "dsd_run_id"))
  # Also include log‑transformed ones if they exist
  log_predictors <- paste0("log_", log_factors)
  numeric_predictors <- unique(c(numeric_predictors, log_predictors))

  # Remove any non‑numeric columns from numeric_predictors
  numeric_predictors <- numeric_predictors[sapply(model_df[, numeric_predictors, drop = FALSE], is.numeric)]

  # Final predictor set: numeric + categorical
  all_predictors <- c(numeric_predictors, categorical)

  # ------------------------------------------------------------------
  # 2. Helper functions for model fitting and evaluation
  # ------------------------------------------------------------------
  rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))

  # ------------------------------------------------------------------
  # 3. Loop over responses
  # ------------------------------------------------------------------
  results_list <- list()

  for (response in responses) {

    cat("\n====================\n")
    cat("Analyzing response:", response, "\n")
    cat("====================\n")

    # Remove rows with missing response
    data_clean <- model_df[!is.na(model_df[[response]]), ]

    # Formula for linear model (full quadratic in numeric, plus factors)
    formula_lm <- as.formula(
      paste(
        response, "~",
        paste(numeric_predictors, collapse = " + "), " + ",
        paste(paste0("I(", numeric_predictors, "^2)"), collapse = " + "), " + ",
        paste(categorical, collapse = " * ")
      )
    )

    # Formula for GAM (smooth for numeric, factor terms)
    smooth_terms <- paste0("s(", numeric_predictors, ")")
    formula_gam <- as.formula(
      paste(
        response, "~",
        paste(smooth_terms, collapse = " + "), " + ",
        paste(categorical, collapse = " + ")
      )
    )

    # ---------- Linear Model ----------
    # ---------- Linear Model ----------
    lm_fit <- try(lm(formula_lm, data = data_clean), silent = TRUE)
    if (!inherits(lm_fit, "try-error")) {
      sink(file.path(output_data_dir, response, "lm_summary.txt"))
      cat("Linear model summary\n")
      print(summary(lm_fit))
      cat("\nType II ANOVA (car::Anova)\n")
      
      # Safe evaluation of Anova
      anova_res <- try(car::Anova(lm_fit, type = 2), silent = TRUE)
      if (!inherits(anova_res, "try-error")) {
        print(anova_res)
      } else {
        cat("\nCould not calculate ANOVA (e.g., model saturated or RSS = 0):\n")
        cat(anova_res, "\n")
      }
      sink()
      pred_lm <- predict(lm_fit, newdata = data_clean)
    } else {
      pred_lm <- NA
      warning("Linear model failed for ", response)
    }

    # ---------- GAM ----------
    gam_fit <- try(gam(formula_gam, data = data_clean, method = "REML"), silent = TRUE)
    if (!inherits(gam_fit, "try-error")) {
      sink(file.path(output_data_dir, response, "gam_summary.txt"))
      cat("GAM summary\n")
      print(summary(gam_fit))
      sink()
      # Save GAM plots (smooth terms)
      png(file.path(output_plots_dir, response, "gam_smooth_terms.png"), width = 1000, height = 800)
      plot(gam_fit, pages = 1, se = TRUE, shade = TRUE)
      dev.off()
      pred_gam <- predict(gam_fit, newdata = data_clean)
    } else {
      pred_gam <- NA
      warning("GAM failed for ", response)
    }

    # ---------- Random Forest ----------
    rf_fit <- try(randomForest(formula_lm, data = data_clean, importance = TRUE, ntree = 500), silent = TRUE)
    if (!inherits(rf_fit, "try-error")) {
      sink(file.path(output_data_dir, response, "rf_importance.txt"))
      cat("Random Forest variable importance\n")
      print(importance(rf_fit))
      sink()
      varImpPlot(rf_fit, main = paste("Variable Importance -", response))
      ggsave(file.path(output_plots_dir, response, "rf_var_importance.png"), width = 8, height = 6)
      pred_rf <- predict(rf_fit, newdata = data_clean)

      # SHAP values for RF (using fastshap)
      if (requireNamespace("fastshap", quietly = TRUE)) {
        X <- data_clean[, all_predictors, drop = FALSE]
        shap <- try(fastshap::explain(
          rf_fit,
          X = X,
          pred_wrapper = function(object, newdata) predict(object, newdata),
          nsim = min(shap_nsim, nrow(X))
        ), silent = TRUE)
        if (!inherits(shap, "try-error")) {
          saveRDS(shap, file.path(output_data_dir, response, "shap_values.rds"))
          # Mean absolute SHAP
          mean_shap <- colMeans(abs(shap), na.rm = TRUE)
          write.csv(data.frame(variable = names(mean_shap), mean_abs_shap = mean_shap),
                    file.path(output_data_dir, response, "mean_shap.csv"), row.names = FALSE)
          # SHAP summary plot using shapviz
          if (requireNamespace("shapviz", quietly = TRUE)) {
            sv <- shapviz(shap, X = X)
            png(file.path(output_plots_dir, response, "shap_summary.png"), width = 800, height = 600)
            print(plot(sv))
            dev.off()
          }
        }
      }

      # Partial dependence plots for top 3 predictors (by importance)
      imp <- importance(rf_fit)
      top_vars <- names(sort(imp[, ncol(imp)], decreasing = TRUE))[1:min(3, ncol(imp)-1)]
      for (var in top_vars) {
        p <- try(partial(rf_fit, pred.var = var, plot = TRUE, plot.engine = "ggplot2"), silent = TRUE)
        if (!inherits(p, "try-error") && !is.null(p)) {
          ggsave(file.path(output_plots_dir, response, paste0("pdp_", var, ".png")), p$plot, width = 6, height = 4)
        }
        # Also try two-way partial for top two variables if they exist
        if (length(top_vars) >= 2) {
          p2 <- try(partial(rf_fit, pred.var = top_vars[1:2], plot = TRUE, plot.engine = "ggplot2"), silent = TRUE)
          if (!inherits(p2, "try-error") && !is.null(p2)) {
            ggsave(file.path(output_plots_dir, response, paste0("pdp_", top_vars[1], "_", top_vars[2], ".png")),
                   p2$plot, width = 6, height = 5)
          }
        }
      }
    } else {
      pred_rf <- NA
      warning("Random Forest failed for ", response)
    }

    # ---------- Gaussian Process (DiceKriging) ----------
    gp_fit <- NULL
    pred_gp <- NA
    if (requireNamespace("DiceKriging", quietly = TRUE)) {
      # Use numeric predictors only (km can't handle factors directly; we convert factors to dummies or use them as is if we want)
      # For simplicity, we'll use only numeric predictors for the GP, ignoring categorical ones.
      # Alternatively, we could create dummy variables, but that complicates interpretation.
      # We'll limit to numeric predictors and warn.
      numeric_preds_only <- setdiff(numeric_predictors, log_predictors) # use original numeric, not log? But we can include log.
      # We'll keep all numeric (including log).
      num_data <- data_clean[, numeric_predictors, drop = FALSE]
      if (ncol(num_data) > 0 && nrow(num_data) > ncol(num_data)) {
        gp_fit <- try(km(design = num_data, response = data_clean[[response]], covtype = "gauss"), silent = TRUE)
        if (!inherits(gp_fit, "try-error")) {
          sink(file.path(output_data_dir, response, "gp_summary.txt"))
          cat("Gaussian Process (DiceKriging) summary\n")
          print(summary(gp_fit))
          sink()
          pred_gp <- try(predict(gp_fit, newdata = num_data, type = "SK")$mean, silent = TRUE)
          if (inherits(pred_gp, "try-error")) pred_gp <- NA
        } else {
          warning("GP failed for ", response)
        }
      } else {
        warning("Not enough numeric data or too many predictors for GP.")
      }
    }

    # ---------- Model comparison ----------
    preds <- list(
      LM = pred_lm,
      GAM = pred_gam,
      RF = pred_rf,
      GP = pred_gp
    )
    rmse_vals <- sapply(preds, function(p) {
      if (length(p) == 1 && is.na(p)) return(NA)
      else rmse(data_clean[[response]], p)
    })
    rmse_df <- data.frame(Model = names(preds), RMSE = rmse_vals, stringsAsFactors = FALSE)
    write.csv(rmse_df, file.path(output_data_dir, response, "model_comparison_rmse.csv"), row.names = FALSE)

    # Create a simple plot of predicted vs observed for each model
    plot_list <- list()
    for (model_name in names(preds)) {
      if (!is.na(preds[[model_name]][1])) {
        p <- ggplot(data_clean, aes(x = .data[[response]], y = preds[[model_name]])) +
          geom_point(alpha = 0.6) +
          geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
          labs(x = "Observed", y = "Predicted", title = paste(model_name, "-", response)) +
          theme_minimal()
        plot_list[[model_name]] <- p
      }
    }
    if (length(plot_list) > 0) {
      grid_plot <- gridExtra::grid.arrange(grobs = plot_list, ncol = 2)
      ggsave(file.path(output_plots_dir, response, "pred_vs_obs.png"), grid_plot, width = 10, height = 8)
    }

    # ---------- Sobol sensitivity indices (optional, using GP or RF surrogate) ----------
    if (run_sobol && (exists("gp_fit") && !inherits(gp_fit, "try-error"))) {
      # Use GP as surrogate for Sobol
      sobol_model <- function(X) {
        predict(gp_fit, newdata = as.data.frame(X), type = "SK")$mean
      }
      # Generate Latin Hypercube samples for numeric predictors only
      n_p <- length(numeric_predictors)
      if (n_p >= 2) {
        X1 <- randomLHS(sobol_n, n_p)
        X2 <- randomLHS(sobol_n, n_p)
        colnames(X1) <- colnames(X2) <- numeric_predictors
        # Scale to the ranges of the data (important for GP)
        ranges <- apply(data_clean[, numeric_predictors, drop = FALSE], 2, range, na.rm = TRUE)
        # Map [0,1] to actual ranges
        X1_scaled <- sweep(X1, 2, ranges[1,], "+") * sweep(1, 2, diff(ranges), "*") # careful: scaling
        X2_scaled <- sweep(X2, 2, ranges[1,], "+") * sweep(1, 2, diff(ranges), "*")
        sobol_est <- try(sobolSalt(model = NULL, X1 = X1_scaled, X2 = X2_scaled, nboot = 100), silent = TRUE)
        if (!inherits(sobol_est, "try-error")) {
          # Evaluate model
          Y <- rbind(X1_scaled, X2_scaled) %>% apply(1, sobol_model)
          tell(sobol_est, Y)
          # Save results
          sobol_df <- data.frame(
            parameter = numeric_predictors,
            first_order = sobol_est$S,
            total_order = sobol_est$T
          )
          write.csv(sobol_df, file.path(output_data_dir, response, "sobol_indices.csv"), row.names = FALSE)
          # Plot
          sobol_plot <- ggplot(sobol_df, aes(x = parameter)) +
            geom_point(aes(y = first_order, color = "First order")) +
            geom_point(aes(y = total_order, color = "Total order")) +
            labs(y = "Sensitivity index", title = paste("Sobol indices -", response)) +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
          ggsave(file.path(output_plots_dir, response, "sobol_indices.png"), sobol_plot, width = 8, height = 6)
        }
      } else {
        warning("Not enough numeric predictors for Sobol analysis (need at least 2).")
      }
    } else if (run_sobol && (!exists("gp_fit") || inherits(gp_fit, "try-error"))) {
      message("Skipping Sobol because GP failed.")
    }

    # Collect results for this response
    results_list[[response]] <- list(
      lm = if (exists("lm_fit") && !inherits(lm_fit, "try-error")) lm_fit else NULL,
      gam = if (exists("gam_fit") && !inherits(gam_fit, "try-error")) gam_fit else NULL,
      rf = if (exists("rf_fit") && !inherits(rf_fit, "try-error")) rf_fit else NULL,
      gp = if (exists("gp_fit") && !inherits(gp_fit, "try-error")) gp_fit else NULL,
      rmse = rmse_df
    )

    cat("Analysis for", response, "complete.\n")
  }

  # ------------------------------------------------------------------
  # 4. Return invisible list
  # ------------------------------------------------------------------
  invisible(results_list)
}

  cat("\nScript finished successfully.\n")
}

run_definitive_screening_design()

