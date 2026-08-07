## ============================================================
## 12_Replicated_Definitive_Screening_Design.R
## Definitive Screening Design (DSD) with replicated simulations
## Generates a full experimental matrix across 9 spatio-temporal regimes
## and runs N_REPLICAS independent replicates for each design point.
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

run_replicated_definitive_screening_design <- function() {
  project_root <- find_project_root()
  source(file.path(project_root, "scripts", "05_motor.R"))

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

  exp_data_dir <- file.path(project_root, "example", "data", "experiments", "definitive_screening")
  exp_plots_dir <- file.path(project_root, "example", "plots", "experiments", "definitive_screening")
  dir.create(exp_data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(exp_plots_dir, recursive = TRUE, showWarnings = FALSE)

  # ------------------------------------------------------------------
  # 1. Define factor levels and scale types
  # ------------------------------------------------------------------
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
  # 5. Create replicate-aware simulation helpers
  # ------------------------------------------------------------------
  if (!exists("simular_hca_replicas", mode = "function")) {
    stop("simular_hca_replicas() is not available. Make sure scripts/05_motor.R defines it.")
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

  run_replicated_design_row <- function(row, n_replicas = 3) {
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

    params <- c(
      params,
      as.list(row[c("tipo_plantilla", "modo_temporal", continuous_factors)]),
      list(theta0 = 0.1)
    )

    replicas_out <- simular_hca_replicas(
      n_replicas = n_replicas,
      semilla_base = params$seed,
      nx = params$nx,
      ny = params$ny,
      n_steps = params$n_steps,
      D = params$D,
      k = params$k,
      dt = params$dt,
      N0 = params$N0,
      N_umbral = params$N_umbral,
      p_max = params$p_max,
      N_half = params$N_half,
      g_inicial = params$g_inicial,
      sigma_e = params$sigma_e,
      mutacion_sd = params$mutacion_sd,
      paso_introduccion = params$paso_introduccion,
      tipo_plantilla = params$tipo_plantilla,
      modo_temporal = params$modo_temporal,
      periodo_temporal = params$periodo_temporal,
      A_max = params$A_max,
      D_A = params$D_A,
      delta_A = params$delta_A,
      tasa_dosificacion = params$tasa_dosificacion,
      theta0 = params$theta0,
      theta1 = params$theta1,
      omega2 = params$omega2,
      K_theta = params$K_theta,
      mort_base = params$mort_base,
      mort_estres = params$mort_estres,
      K_mort = params$K_mort,
      Nc = params$Nc,
      guardar_cada = params$guardar_cada
    )

    do.call(rbind, lapply(seq_along(replicas_out), function(j) {
      resp <- extract_simulation_responses(replicas_out[[j]])
      cbind(
        row[c("sim_id", "regime_id", "dsd_run_id", "tipo_plantilla", "modo_temporal", continuous_factors)],
        replica_id = j,
        resp
      )
    }))
  }

  # ------------------------------------------------------------------
  # 6. Run the replicated design and export the results
  # ------------------------------------------------------------------
  N_REPLICAS <- 3

  response_list <- vector("list", nrow(full_simulation_design))
  total_runs <- nrow(full_simulation_design) * N_REPLICAS

  cat("Starting replicated definitive-screening experiment:\n")
  cat(sprintf("  %d design points × %d replicates = %d simulations\n",
              nrow(full_simulation_design), N_REPLICAS, total_runs))
  cat("Estimated time: ~", round(total_runs * 0.5, 1), "minutes (depending on your machine)\n\n")

  pb <- txtProgressBar(min = 0, max = total_runs, style = 3)
  run_counter <- 0

  for (i in seq_len(nrow(full_simulation_design))) {
    response_list[[i]] <- run_replicated_design_row(full_simulation_design[i, , drop = FALSE], n_replicas = N_REPLICAS)
    run_counter <- run_counter + N_REPLICAS
    setTxtProgressBar(pb, run_counter)
  }

  close(pb)

  response_df <- do.call(rbind, response_list)
  rownames(response_df) <- NULL

  output_csv <- file.path(exp_data_dir, "dsd_simulation_design_replicated.csv")
  output_xlsx <- file.path(exp_data_dir, "dsd_simulation_design_replicated.xlsx")
  output_rds <- file.path(exp_data_dir, "dsd_simulation_design_replicated.rds")

  write.csv(response_df, output_csv, row.names = FALSE)
  saveRDS(response_df, output_rds)

  if (file.exists(output_xlsx)) {
    file.remove(output_xlsx)
  }
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    openxlsx::write.xlsx(response_df, output_xlsx, rowNames = FALSE)
    cat(" -", normalizePath(output_xlsx), "\n")
  } else {
    message("openxlsx is not installed; skipped Excel export.")
  }

  cat("Saved replicated design results with", nrow(response_df), "rows to:\n")
  cat(" -", normalizePath(output_csv), "\n")
  cat(" -", normalizePath(output_rds), "\n")

  invisible(response_df)
}

# Run the workflow when the file is sourced directly.
if (sys.nframe() == 0) {
  run_replicated_definitive_screening_design()
}
