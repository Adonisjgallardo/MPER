## ============================================================
## 09_profile.R
## Performance profiling for the HCA evolutionary rescue simulation.
## Uses base R profiling with Rprof() and optionally profvis if installed.
## ============================================================

# Load simulation engine
source("scripts/05_motor.R")

# Example configuration for a moderate test run.
params <- list(
  nx = 200,
  ny = 200,
  n_steps = 100,
  D = 0.18,
  k = 0.13,
  dt = 1,
  N0 = 1.0,
  N_umbral = 0.15,
  p_max = 0.6,
  N_half = 0.3,
  g_inicial = 0.1,
  sigma_e = 0.05,
  sigma_g = 0.01,
  mutacion_sd = 0.02,
  paso_introduccion = 40,
  tipo_plantilla = "lineal",
  modo_temporal = "constante",
  periodo_temporal = 10,
  A_max = 1.0,
  D_A = 0.12,
  delta_A = 0.02,
  tasa_dosificacion = 0.15,
  theta0 = 0.1,
  theta1 = 0.8,
  omega2 = 0.05,
  K_theta = 0.3,
  mort_base = 0.01,
  mort_estres = 0.9,
  K_mort = 0.3,
  Nc = 200,
  guardar_cada = 10,
  seed = 123,
  guardar_historial_en_disco = TRUE,
  historial_dir = "profile_historial"
)

profile_file <- "profile_simular_hca.out"

cat("Profiling simulation with Rprof() to ", profile_file, "\n")
Rprof(profile_file, interval = 0.01)
res <- do.call(simular_hca, params)
Rprof(NULL)
cat("Profiling complete. Use summaryRprof(\"", profile_file, "\") to inspect.\n", sep = "")

if (requireNamespace("profvis", quietly = TRUE)) {
  cat("profvis is installed. Generating profvis object...\n")
  library(profvis)
  pv <- profvis({
    res <- do.call(simular_hca, params)
  })
  print(pv)
} else {
  cat("profvis is not installed. Install it to generate an interactive profile.\n")
}
