## ============================================================
## 08_simulate.R nb
## Script de ejecución de simulación, exportación de datos y
## generación de gráficas con métricas.
## The goal is to simulate evolutionary rescue under an abrupt antibiotic shock, 
## with quantitative genetics (heritable genetic value g, non‑heritable environmental deviation e, phenotype z = g + e) 
##and stabilising selection around a local optimum θ(x,y,t) that depends on local antibiotic concentration
## ============================================================

output_data_dir <- file.path("example", "data", "simulation")
output_plots_dir <- file.path("example", "plots", "simulation")

dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_plots_dir, recursive = TRUE, showWarnings = FALSE)

source("scripts/05_motor.R")
source("scripts/07_analisis.R")

# progress bar
steps <- c(
  "Ejecutar simulación",
  "Guardar resúmenes y parámetros",
  "Calcular dimensión fractal",
  "Comparar contra analítica",
  "Generar gráficos"
)
pb <- txtProgressBar(min = 0, max = length(steps), style = 3)
setTxtProgressBar(pb, 0)
flush.console()

# Parámetros de ejemplo. Ajusta según necesites.
#
# Experimental factors are grouped by mechanistic role.
# Non-experimental parameters are treated as cofactors or simulation controls.
params <- list(
  # --- Spatio-temporal heterogeneity (refuges and fluctuation) ---
  tipo_plantilla = "lineal",     # uniform / lineal / radial
  modo_temporal = "constante",   # constante / pulsos / sinusoidal
  periodo_temporal = 10,
  paso_introduccion = 250,

  # --- Selection landscape (treatment factors) ---
  A_max = 1.0,
  theta0 = 0.1,
  theta1 = 0.8,
  omega2 = 0.05,
  K_theta = 0.3,

  # --- Evolvability & genetic architecture (raw material) ---
  g_inicial = 0.1,
  sigma_e = 0.05,
  mutacion_sd = 0.02,

  # --- Demographic context (carrying capacity & drift) ---
  N0 = 1.0,
  D = 0.18,
  k = 0.13,
  p_max = 0.6,

  # --- Cofactors and simulation controls ---
  nx = 600,
  ny = 600,
  n_steps = 500,
  dt = 1,
  N_umbral = 0.15,
  N_half = 0.3,
  D_A = 0.12,
  delta_A = 0.02,
  tasa_dosificacion = 0.15,
  mort_base = 0.01,
  mort_estres = 0.9,
  K_mort = 0.3,
  Nc = 200,
  guardar_cada = 10,
  seed = 123
)

# Experimental factors to vary in a design of experiments.
factors <- list(
  selection_landscape = list(A_max = params$A_max,
                             omega2 = params$omega2,
                             theta1 = params$theta1,
                             K_theta = params$K_theta),
  evolvability = list(mutacion_sd = params$mutacion_sd,
                      sigma_e = params$sigma_e,
                      g_inicial = params$g_inicial),
  demographic = list(N0 = params$N0,
                     D = params$D,
                     k = params$k,
                     p_max = params$p_max),
  spatio_temporal = list(tipo_plantilla = params$tipo_plantilla,
                         modo_temporal = params$modo_temporal,
                         periodo_temporal = params$periodo_temporal,
                         paso_introduccion = params$paso_introduccion)
)

cofactors <- list(nx = params$nx,
                  ny = params$ny,
                  n_steps = params$n_steps,
                  dt = params$dt,
                  N_umbral = params$N_umbral,
                  N_half = params$N_half,
                  D_A = params$D_A,
                  delta_A = params$delta_A,
                  tasa_dosificacion = params$tasa_dosificacion,
                  mort_base = params$mort_base,
                  mort_estres = params$mort_estres,
                  K_mort = params$K_mort,
                  Nc = params$Nc,
                  guardar_cada = params$guardar_cada,
                  seed = params$seed)

# Resistance responses of interest:
# - g_bar: mean genetic value (adaptive response)
# - z_bar: mean realized phenotype (actual resistance phenotype)
# - var_g: genetic variance (evolvability and standing variation)
# - persistence: population persistence / extinction timing
#
# These are the primary response variables for rescue analysis.

# Paso 1: ejecución de la simulación
execution_time <- system.time({
  sim <- do.call(simular_hca, params)
})
print(execution_time)
setTxtProgressBar(pb, 1)
flush.console()

# Paso 2: guardado de resúmenes, parámetros y respuestas
write.csv(sim$resumen,
          file.path(output_data_dir, "resumen.csv"),
          row.names = FALSE)
write.csv(data.frame(name = names(params), value = unlist(params)),
          file.path(output_data_dir, "parametros.csv"),
          row.names = FALSE)
write.csv(data.frame(
  factor = names(factors),
  values = vapply(factors, function(x) paste0(names(x), "=", unlist(x), collapse = "; "), "")),
  file.path(output_data_dir, "factores_experimentales.csv"),
  row.names = FALSE)
write.csv(data.frame(
  cofactor = names(cofactors),
  value = unlist(cofactors)
),
file.path(output_data_dir, "cofactores_y_controles.csv"),
row.names = FALSE)

persistence_time <- if (is.na(sim$paso_extincion)) {
  max(sim$resumen$paso) - sim$paso_introduccion
} else {
  sim$paso_extincion - sim$paso_introduccion
}

response_metrics <- data.frame(
  response = c("final_g_bar", "final_z_bar", "final_var_g", "persistence_time", "extinct"),
  value = c(
    tail(sim$resumen$g_bar, 1),
    tail(sim$resumen$z_bar, 1),
    tail(sim$resumen$var_g, 1),
    persistence_time,
    !is.na(sim$paso_extincion)
  )
)
write.csv(response_metrics,
          file.path(output_data_dir, "respuestas.csv"),
          row.names = FALSE)
setTxtProgressBar(pb, 2)
flush.console()

# Paso 3: cálculo y guardado de dimensión fractal
fractal_por_paso <- do.call(rbind, lapply(sim$historial, function(s) {
  df_fractal <- calcular_dimension_fractal(s$estado)
  data.frame(paso = s$paso, dimension_fractal = df_fractal$dimension)
}))
write.csv(fractal_por_paso,
          file.path(output_data_dir, "dimension_fractal_por_paso.csv"),
          row.names = FALSE)
fractal_final <- calcular_dimension_fractal(sim$estado_final)
write.csv(data.frame(dimension_fractal_final = fractal_final$dimension),
          file.path(output_data_dir, "dimension_fractal_final.csv"),
          row.names = FALSE)
setTxtProgressBar(pb, 3)
flush.console()

# Paso 4: comparaciones analíticas
idx_intro <- which(sim$resumen$paso == sim$paso_introduccion)
sigma_g2_empirico <- if (length(idx_intro) == 1) sim$resumen$var_g[idx_intro] else sim$resumen$var_g[1]
sigma_e2 <- params$sigma_e^2
comparacion_Nt <- comparar_Nt_caja3(sim,
                                     Wmax = 1,
                                     sigma_g2 = sigma_g2_empirico,
                                     sigma_e2 = sigma_e2)
write.csv(comparacion_Nt,
          file.path(output_data_dir, "comparacion_Nt_caja3.csv"),
          row.names = FALSE)
comparacion_gbar <- comparar_gbar_caja3(sim,
                                         k = 0,
                                         sigma_g2 = sigma_g2_empirico,
                                         sigma_e2 = sigma_e2)
write.csv(comparacion_gbar,
          file.path(output_data_dir, "comparacion_gbar_caja3.csv"),
          row.names = FALSE)
setTxtProgressBar(pb, 4)
flush.console()

# Paso 5: crear gráficos
graficar_estado_hca(sim$estado_final,
                     sim$N_final,
                     sim$resist_final,
                     archivo = file.path(output_plots_dir, "estado_final.png"),
                     titulo = "Estado final de la simulación",
                     variable = "resist",
                     mostrar_fractal = TRUE)
animar_historial_hca(sim$historial,
                     archivo_gif = file.path(output_plots_dir, "evolucion.gif"),
                     archivo_mp4 = file.path(output_plots_dir, "evolucion.mp4"),
                     variable = "resist",
                     fps = 10,
                     width = 600,
                     height = 600,
                     mostrar_fractal = TRUE)
setTxtProgressBar(pb, 5)
flush.console()

close(pb)

cat("Simulación completada:\n")
cat(" - Datos guardados en:", normalizePath(output_data_dir), "\n")
cat(" - Gráficos guardados en:", normalizePath(output_plots_dir), "\n")
