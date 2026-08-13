## ============================================================
## utils.R
## Funciones auxiliares compartidas por la app Shiny (no forman
## parte del motor de simulacion en si).
## ============================================================

#' Recorta un valor al rango [lo, hi]
clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

#' Tema ggplot compartido por todos los graficos de la app, para
#' mantener consistencia visual entre pestañas.
tema_hca <- function() {
  ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey40"),
      legend.position = "right"
    )
}

#' Construye una curva de supervivencia estilo Kaplan-Meier a partir
#' de tiempos de persistencia empiricos (data.frame con columnas
#' `tiempo_persistencia` y `extinguio`, tal como devuelve
#' `tiempos_persistencia_empiricos()` en 06_genetica_cuantitativa.R).
#'
#' No se usa el paquete `survival` (para no añadir una dependencia
#' mas); se calcula el estimador producto-limite "a mano", que con
#' datos sin empates ni riesgos competitivos se reduce a una simple
#' decreciente escalonada.
#'
#' @param tp data.frame con columnas tiempo_persistencia, extinguio
#' @return data.frame con columnas tiempo, superviviencia (0-1),
#'         listo para geom_step()
curva_supervivencia_km <- function(tp) {
  tp <- tp[order(tp$tiempo_persistencia), ]
  n_total <- nrow(tp)
  if (n_total == 0) {
    return(data.frame(tiempo = numeric(0), supervivencia = numeric(0)))
  }

  tiempos_evento <- sort(unique(tp$tiempo_persistencia[tp$extinguio]))
  s <- 1
  filas <- list(data.frame(tiempo = 0, supervivencia = 1))
  en_riesgo <- n_total

  for (te in tiempos_evento) {
    en_riesgo <- sum(tp$tiempo_persistencia >= te)
    d_i <- sum(tp$tiempo_persistencia == te & tp$extinguio)
    s <- s * (1 - d_i / en_riesgo)
    filas[[length(filas) + 1]] <- data.frame(tiempo = te, supervivencia = s)
  }

  do.call(rbind, filas)
}

#' Formatea un numero con separador de miles, o "-" si es NA
fmt_num <- function(x, digits = 2) {
  if (length(x) == 0 || is.na(x)) return("-")
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

#' Empaqueta todos los parametros de simular_hca() leidos desde los
#' inputs de la UI en una lista lista para usar con do.call().
#' Centralizar esto evita repetir el mismo bloque largo en varios
#' observeEvent().
recolectar_parametros_sim <- function(input) {
  list(
    nx = input$nx, ny = input$nx,
    n_steps = input$n_steps,
    D = input$D, k = input$k, N0 = input$N0,
    N_umbral = input$N_umbral, p_max = input$p_max, N_half = input$N_half,
    g_inicial = input$g_inicial, sigma_e = input$sigma_e,
    mutacion_sd = input$mutacion_sd,
    paso_introduccion = input$paso_introduccion,
    tipo_plantilla = input$tipo_plantilla,
    modo_temporal = input$modo_temporal,
    periodo_temporal = input$periodo_temporal,
    A_max = input$A_max, D_A = input$D_A, delta_A = input$delta_A,
    tasa_dosificacion = input$tasa_dosificacion,
    theta0 = input$theta0, theta1 = input$theta1,
    omega2 = input$omega2, K_theta = input$K_theta,
    mort_base = input$mort_base, mort_estres = input$mort_estres,
    K_mort = input$K_mort,
    Nc = input$Nc,
    guardar_cada = if (isTRUE(input$mantener_historial)) input$guardar_cada else input$n_steps + 1L,
    seed = sample.int(1e6, 1)
  )
}
