## ============================================================
## 10_implement.R
## Full‑factorial core experimental design
## Factors: A_max (3) × omega2 (3) × mutacion_sd (3)
## Replicates per combination: 3
## Response: extinction time, final N, final g_bar, final var_g
## ============================================================

# ------------------------------------------------------------------
# 0. Setup environment and source core functions
# ------------------------------------------------------------------
# Assumes this script can be run either from the project root or from inside `scripts/`.
# Detect a sensible project root and construct paths relative to it so outputs are
# written to the repository's top-level `data/` and `plots/` folders regardless of
# the current working directory used to run the script.
if (file.exists("scripts/05_motor.R")) {
  project_root <- "."
} else if (file.exists("../scripts/05_motor.R")) {
  project_root <- ".."
} else {
  # Fallback: assume current working directory is the project root
  project_root <- "."
}

source(file.path(project_root, "scripts", "05_motor.R"))
source(file.path(project_root, "scripts", "07_analisis.R"))   # only needed for extra post‑proc, but harmless

# Create output directories for this experiment (project-root relative)
exp_data_dir <- file.path(project_root, "data", "experiments", "factorial_core")
exp_plots_dir <- file.path(project_root, "plots", "experiments", "factorial_core")
dir.create(exp_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(exp_plots_dir, recursive = TRUE, showWarnings = FALSE)

# If results already exist, load them and skip the heavy simulation loop
results_rds_path <- file.path(exp_data_dir, "factorial_results_full.rds")
if (file.exists(results_rds_path)) {
  message("Detected existing results RDS; loading and skipping simulations.")
  results_df <- readRDS(results_rds_path)
  if (!"rescued" %in% names(results_df)) results_df$rescued <- !results_df$went_extinct
  skip_simulations <- TRUE
} else {
  skip_simulations <- FALSE
}
# ------------------------------------------------------------------
# 1. Define factor levels (3 × 3 × 3 = 27 combinations)
# ------------------------------------------------------------------
# A_max       : stress severity  (Low, Medium, High)
# omega2      : selection width  (Sharp, Medium, Wide)  -> narrower = stronger selection
# mutacion_sd : evolvability     (Low, Medium, High)    -> mutation step size
#
# These ranges are chosen to cover the biologically relevant space
# while ensuring numerical stability.
design <- expand.grid(
  A_max       = c(0.5, 1.0, 1.5),
  omega2      = c(0.02, 0.05, 0.10),
  mutacion_sd = c(0.005, 0.02, 0.08),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------
# 2. Fixed simulation parameters (common to all runs)
# ------------------------------------------------------------------
# Use a moderate grid (150×150) and 500 steps to keep runtime
# manageable for 81 simulations (27 combos × 3 reps).
# Increase these for production runs if you have more computing power.
base_params <- list(
  nx = 150,
  ny = 150,
  n_steps = 500,
  D = 0.18,
  k = 0.13,
  dt = 1,
  N0 = 1.0,
  N_umbral = 0.15,
  p_max = 0.6,
  N_half = 0.3,
  g_inicial = 0.1,
  sigma_e = 0.05,
  paso_introduccion = 250,        # antibiotic introduced halfway
  tipo_plantilla = "lineal",      # spatial gradient
  modo_temporal = "constante",    # constant stress post‑shock
  D_A = 0.12,
  delta_A = 0.02,
  tasa_dosificacion = 0.15,
  theta0 = 0.1,
  theta1 = 0.8,
  K_theta = 0.3,
  mort_base = 0.01,
  mort_estres = 0.9,
  K_mort = 0.3,
  Nc = 200,                       # critical threshold for "quasi‑extinction"
  guardar_cada = 25               # snapshot frequency (you can increase to save memory)
)

# ------------------------------------------------------------------
# 3. Helper: extract summary statistics from a replicate run
# ------------------------------------------------------------------
extract_replicate_summary <- function(replica_output, combo_id, factor_values) {
  # replica_output is one element from the list returned by simular_hca_replicas()
  res <- replica_output$resumen
  
  # If the simulation went extinct, the last row is the extinction step.
  # Otherwise, it's the final step of the simulation.
  last_row <- tail(res, 1)
  
  data.frame(
    combo_id          = combo_id,
    A_max             = factor_values$A_max,
    omega2            = factor_values$omega2,
    mutacion_sd       = factor_values$mutacion_sd,
    replica           = NA,  # will be filled in the loop
    paso_extincion    = replica_output$paso_extincion,   # NA if never crossed Nc
    final_N           = last_row$n_ocupados,
    final_g_bar       = last_row$g_bar,
    final_var_g       = last_row$var_g,
    final_z_bar       = last_row$z_bar,
    went_extinct      = !is.na(replica_output$paso_extincion),
    N_at_intro        = res$n_ocupados[res$paso == base_params$paso_introduccion][1],
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------
# 4. Launch the factorial experiment
# ------------------------------------------------------------------
if (!skip_simulations) {
  set.seed(12345)   # overall seed for reproducibility (seeds inside replicates are offset)

  all_results <- list()
  total_combos <- nrow(design)
  total_runs   <- total_combos * 3   # 3 replicates per combo

  cat("Starting full‑factorial experiment:\n")
  cat(sprintf("  %d combinations × %d replicates = %d simulations\n",
              total_combos, 3, total_runs))
  cat("Estimated time: ~", round(total_runs * 0.5, 1), "minutes (depending on your machine)\n\n")

  pb <- txtProgressBar(min = 0, max = total_runs, style = 3)
  run_counter <- 0

  for (i in 1:total_combos) {

    # Build parameter list for this combination
    params_i <- base_params
    params_i$A_max       <- design$A_max[i]
    params_i$omega2      <- design$omega2[i]
    params_i$mutacion_sd <- design$mutacion_sd[i]

    # Run 3 independent replicates with a unique base seed per combo
    semilla_base <- 10000 + i * 10
    replicas_out <- simular_hca_replicas(
      n_replicas = 3,
      semilla_base = semilla_base,
      nx = params_i$nx,
      ny = params_i$ny,
      n_steps = params_i$n_steps,
      D = params_i$D,
      k = params_i$k,
      dt = params_i$dt,
      N0 = params_i$N0,
      N_umbral = params_i$N_umbral,
      p_max = params_i$p_max,
      N_half = params_i$N_half,
      g_inicial = params_i$g_inicial,
      sigma_e = params_i$sigma_e,
      mutacion_sd = params_i$mutacion_sd,
      paso_introduccion = params_i$paso_introduccion,
      tipo_plantilla = params_i$tipo_plantilla,
      modo_temporal = params_i$modo_temporal,
      periodo_temporal = params_i$periodo_temporal,  # not used for 'constante'
      A_max = params_i$A_max,
      D_A = params_i$D_A,
      delta_A = params_i$delta_A,
      tasa_dosificacion = params_i$tasa_dosificacion,
      theta0 = params_i$theta0,
      theta1 = params_i$theta1,
      omega2 = params_i$omega2,
      K_theta = params_i$K_theta,
      mort_base = params_i$mort_base,
      mort_estres = params_i$mort_estres,
      K_mort = params_i$K_mort,
      Nc = params_i$Nc,
      guardar_cada = params_i$guardar_cada
    )

    # Extract summary for each replicate
    for (j in 1:length(replicas_out)) {
      summary_row <- extract_replicate_summary(
        replica_output = replicas_out[[j]],
        combo_id = i,
        factor_values = list(
          A_max = design$A_max[i],
          omega2 = design$omega2[i],
          mutacion_sd = design$mutacion_sd[i]
        )
      )
      summary_row$replica <- j
      all_results[[length(all_results) + 1]] <- summary_row

      run_counter <- run_counter + 1
      setTxtProgressBar(pb, run_counter)
    }
  }

  close(pb)
  cat("\nAll simulations completed!\n")
}

# ------------------------------------------------------------------
# 5. Combine results and save to disk (or use existing RDS)
# ------------------------------------------------------------------
if (!skip_simulations) {
  results_df <- do.call(rbind, all_results)
  rownames(results_df) <- NULL

  # Add derived columns for easier plotting
  results_df$rescued <- !results_df$went_extinct
  # Compute whether it stayed above Nc at the very end
  results_df$final_above_Nc <- results_df$final_N >= base_params$Nc

  # Save the full dataset
  write.csv(results_df,
            file.path(exp_data_dir, "factorial_results_full.csv"),
            row.names = FALSE)

  # Also save as RDS (preserves data types and is faster to load)
  saveRDS(results_df,
          file.path(exp_data_dir, "factorial_results_full.rds"))

} else {
  # results_df was loaded from existing RDS earlier. Ensure derived columns exist
  if (!"rescued" %in% names(results_df)) results_df$rescued <- !results_df$went_extinct
  if (!"final_above_Nc" %in% names(results_df)) {
    if ("final_N" %in% names(results_df)) {
      results_df$final_above_Nc <- results_df$final_N >= base_params$Nc
    } else {
      results_df$final_above_Nc <- NA
    }
  }
}

# ------------------------------------------------------------------
# 6b. Extended plotting & modelling (safe: will load RDS if results_df missing)
# ------------------------------------------------------------------
{
  # ensure results_df exists (in case we loaded earlier)
  if (!exists("results_df")) {
    if (file.exists(results_rds_path)) {
      results_df <- readRDS(results_rds_path)
    } else stop("results_df not found and RDS missing")
  }

  # Load useful packages (fail gracefully)
  pkgs <- c("ggplot2", "dplyr", "emmeans", "ggeffects", "cowplot", "scales", "car")
  have <- sapply(pkgs, require, quietly = TRUE, character.only = TRUE)
  if (!all(have[c("ggplot2","dplyr")] )) {
    warning("Please install 'ggplot2' and 'dplyr' to create plots")
  }

  # Prepare factors
  results_df <- results_df %>%
    mutate(
      A_max_f = factor(A_max),
      omega2_f = factor(omega2),
      mutacion_sd_f = factor(mutacion_sd)
    )

  # 1) Extinction probability (grouped summary)
  ext_prob <- results_df %>%
    group_by(A_max_f, omega2_f, mutacion_sd_f) %>%
    summarise(prop_extinct = mean(!rescued), n = n(), .groups = "drop")

  p1 <- ggplot(ext_prob, aes(x = A_max_f, y = prop_extinct, color = mutacion_sd_f, group = mutacion_sd_f)) +
    geom_point(position = position_dodge(width = 0.4), size = 3) +
    geom_line(position = position_dodge(width = 0.4)) +
    facet_grid(. ~ omega2_f, labeller = label_both) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
    labs(x = "A_max", y = "Extinction probability", color = "mutacion_sd") +
    theme_minimal()

  ggsave(file.path(exp_plots_dir, "ext_prob_by_Amax_extended.png"), p1, width = 10, height = 5, dpi = 150)

  # 2) Fit models (if packages available)
  car_avail <- isTRUE(have["car"])

  if (requireNamespace("stats", quietly = TRUE)) {
    glm_bin <- try(stats::glm(rescued ~ A_max_f * omega2_f * mutacion_sd_f, data = results_df, family = binomial), silent = TRUE)
    if (!inherits(glm_bin, "try-error")) {
      # Save summary to text
      capture.output(summary(glm_bin), file = file.path(exp_data_dir, "glm_bin_summary.txt"))
      if (car_avail) capture.output(car::Anova(glm_bin, type = 2), file = file.path(exp_data_dir, "glm_bin_anova.txt"))

      if (have["emmeans"]) {
        emm_bin <- try(emmeans::emmeans(glm_bin, ~ A_max_f * mutacion_sd_f | omega2_f, type = "response"), silent = TRUE)
        if (!inherits(emm_bin, "try-error")) {
          emm_df <- try(as.data.frame(summary(emm_bin, infer = TRUE)), silent = TRUE)
          if (!inherits(emm_df, "try-error")) {
            if (!all(c("lower.CL", "upper.CL") %in% names(emm_df)) &&
                all(c("asymp.LCL", "asymp.UCL") %in% names(emm_df))) {
              emm_df$lower.CL <- emm_df$asymp.LCL
              emm_df$upper.CL <- emm_df$asymp.UCL
            }
            if (all(c("prob", "lower.CL", "upper.CL") %in% names(emm_df))) {
              p_emm <- ggplot(emm_df, aes(x = A_max_f, y = prob, color = mutacion_sd_f, group = mutacion_sd_f)) +
                geom_point() + geom_line() +
                geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
                facet_wrap(~ omega2_f) + labs(y = "Predicted survival probability") + theme_minimal()
              try(ggsave(file.path(exp_plots_dir, "emm_glm.png"), p_emm, width = 10, height = 6, dpi = 150), silent = TRUE)
            } else {
              warning("Skipping emmeans plot because the computed emmeans result does not include lower.CL/upper.CL columns.")
            }
          }
        }
      }
    }
  }

  # 3) Continuous response example: final_N
  p2 <- ggplot(results_df, aes(x = A_max_f, y = final_N, fill = mutacion_sd_f)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    geom_jitter(aes(color = mutacion_sd_f), width = 0.15, alpha = 0.6, size = 1) +
    facet_wrap(~ omega2_f) + theme_minimal()
  ggsave(file.path(exp_plots_dir, "finalN_by_factors.png"), p2, width = 10, height = 6, dpi = 150)

  # 4) Save a small model for final_N
  lm_N <- try(lm(final_N ~ A_max_f * omega2_f * mutacion_sd_f, data = results_df), silent = TRUE)
  if (!inherits(lm_N, "try-error")) capture.output(summary(lm_N), file = file.path(exp_data_dir, "lm_finalN_summary.txt"))

  message("Plots and model summaries saved to:")
  message(" - ", normalizePath(exp_plots_dir))
  message(" - ", normalizePath(exp_data_dir))
}

# ------------------------------------------------------------------
# 6. Quick preliminary visualisation / ANOVA (optional)
# ------------------------------------------------------------------
cat("\nQuick summary of extinction rates per factor:\n")
print(aggregate(rescued ~ A_max + omega2 + mutacion_sd,
                data = results_df, FUN = mean))

# Save a simple plot of extinction probability vs A_max, faceted by omega2
if (require(ggplot2, quietly = TRUE)) {
  p <- ggplot(results_df, aes(x = factor(A_max), fill = factor(mutacion_sd))) +
    geom_bar(aes(y = after_stat(prop), group = factor(mutacion_sd)),
             position = position_dodge(), stat = "count") +
    facet_grid(~ omega2, labeller = label_both) +
    labs(title = "Extinction probability by factor combination",
         x = "Antibiotic max (A_max)", y = "Proportion extinct",
         fill = "Mutation sd") +
    theme_minimal()
  
  ggsave(file.path(exp_plots_dir, "extinction_probability.png"),
         plot = p, width = 10, height = 6, dpi = 150)
}

cat("\nAll results saved to:\n")
cat(" -", normalizePath(exp_data_dir), "\n")
cat("Plots saved to:\n")
cat(" -", normalizePath(exp_plots_dir), "\n")

# ------------------------------------------------------------------
# 7. Example: load and run a linear model on the results
# ------------------------------------------------------------------
# Uncomment this block if you want to test the main effects immediately
#
# model <- glm(rescued ~ A_max * omega2 * mutacion_sd,
#              data = results_df, family = binomial)
# summary(model)
#
# cat("\nANOVA‑type summary (Type II Wald tests):\n")
# print(car::Anova(model))

cat("\nScript finished successfully.\n")


