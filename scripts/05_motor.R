## ============================================================
## 05_motor.R
## Motor principal: acopla PDE (nutriente) <-> CA (crecimiento)
## Modelo HCA nutriente-limitado PURO (sin antibiótico todavía)
## ============================================================

source_scripts <- c("01_campos.R", "02_colonia.R",
                     "03_pde_difusion.R", "04_ca_crecimiento.R")
for (s in source_scripts) source(file.path("scripts", s))

#' Ejecuta la simulación HCA completa
#'
#' @param nx,ny dimensiones de la retícula
#' @param n_steps número de pasos de tiempo
#' @param D coeficiente de difusión del nutriente
#' @param k tasa máxima de consumo por sitio ocupado
#' @param dt paso de tiempo (verificado contra estabilidad FTCS)
#' @param N0 nutriente inicial homogéneo
#' @param N_umbral umbral de nutriente para fenotipo proliferativo
#' @param p_max prob. máxima de división por paso
#' @param mutacion_sd sd de la mutación del rasgo resistencia por división
#' @param guardar_cada cada cuántos pasos se guarda un snapshot (memoria)
#' @param seed semilla
#'
#' @return lista: `historial` (snapshots), `resumen` (data.frame por paso)
simular_hca <- function(nx = 100, ny = 100,
                         n_steps = 300,
                         D = 0.20, k = 0.15, dt = 1,
                         N0 = 1.0,
                         N_umbral = 0.15, p_max = 0.6, N_half = 0.3,
                         mutacion_sd = 0.02,
                         guardar_cada = 10,
                         seed = 123) {

  verificar_estabilidad_dt(D, dt)

  N <- crear_campo_nutriente(nx, ny, N0 = N0)
  colonia <- crear_colonia_inicial(nx, ny, resist_inicial = 0.1, seed = seed)
  estado <- colonia$estado
  resist <- colonia$resist

  historial <- list()
  resumen <- data.frame(paso = integer(0), n_ocupados = integer(0),
                         n_frente = integer(0), N_medio = numeric(0),
                         resist_media = numeric(0))

  for (paso in seq_len(n_steps)) {

    N <- actualizar_nutriente(N, estado, D = D, k = k, dt = dt, N_half = N_half)

    res_ca <- crecer_colonia(estado, resist, N,
                              N_umbral = N_umbral, p_max = p_max,
                              N_half = N_half, mutacion_sd = mutacion_sd)
    estado <- res_ca$estado
    resist <- res_ca$resist

    n_ocupados <- sum(estado)
    n_frente   <- sum(res_ca$fenotipo == "proliferativo", na.rm = TRUE)

    resumen <- rbind(resumen, data.frame(
      paso = paso,
      n_ocupados = n_ocupados,
      n_frente = n_frente,
      N_medio = mean(N),
      resist_media = mean(resist, na.rm = TRUE)
    ))

    if (paso %% guardar_cada == 0 || paso == n_steps) {
      historial[[length(historial) + 1]] <- list(
        paso = paso, estado = estado, N = N,
        resist = resist, fenotipo = res_ca$fenotipo
      )
    }
  }

  list(historial = historial, resumen = resumen,
       estado_final = estado, N_final = N, resist_final = resist)
}

# ------------------------------------------------------------
# Ejecución de ejemplo (descomentar / adaptar al pipeline)
# ------------------------------------------------------------
sim <- simular_hca(nx = 100, ny = 100, n_steps = 300)
saveRDS(sim, "data/hca_nutriente_puro.rds")
