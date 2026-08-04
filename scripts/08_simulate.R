## ============================================================
## 08_simulate.R
## Script de ejecución de simulación, exportación de datos y
## generación de gráficas con métricas.
## ============================================================

output_data_dir <- file.path("data", "simulation")
output_plots_dir <- file.path("plots", "simulation")

dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_plots_dir, recursive = TRUE, showWarnings = FALSE)

source("../scripts/05_motor.R")
source("../scripts/07_analisis.R")

# Parámetros de ejemplo. Ajusta según necesites.
params <- list(
  nx = 1200,
  ny = 1200,
  n_steps = 1000,
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
  paso_introduccion = 400,
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
  seed = 123
)

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

# Paso 1: ejecución de la simulación
sim <- do.call(simular_hca, params)
setTxtProgressBar(pb, 1)
flush.console()

# Paso 2: guardado de resúmenes y parámetros
write.csv(sim$resumen,
          file.path(output_data_dir, "resumen.csv"),
          row.names = FALSE)
write.csv(data.frame(name = names(params), value = unlist(params)),
          file.path(output_data_dir, "parametros.csv"),
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
sigma_g2 <- params$sigma_g^2
sigma_e2 <- params$sigma_e^2
comparacion_Nt <- comparar_Nt_caja3(sim,
                                     Wmax = 1,
                                     sigma_g2 = sigma_g2,
                                     sigma_e2 = sigma_e2)
write.csv(comparacion_Nt,
          file.path(output_data_dir, "comparacion_Nt_caja3.csv"),
          row.names = FALSE)
comparacion_gbar <- comparar_gbar_caja3(sim,
                                         k = 0,
                                         sigma_g2 = sigma_g2,
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
