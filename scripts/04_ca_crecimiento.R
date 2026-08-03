## ============================================================
## 04_ca_crecimiento.R
## Componente CA: división celular limitada por nutriente local
## + zonación fenotípica DERIVADA del campo (no heredada)
## + herencia/mutación del rasgo de resistencia (por-celda)
## 
## Vecindad: von Neumann (4 vecinos), fronteras PERIÓDICAS.
## ============================================================

#' Índices de los 4 vecinos von Neumann con envoltura periódica (toroide)
#'
#' @return matriz 4x2 de coordenadas (fila, columna)
vecinos_von_neumann <- function(f, c, nx, ny) {
  rbind(
    c(if (f == 1) nx else f - 1, c),   # arriba
    c(if (f == nx) 1 else f + 1, c),   # abajo
    c(f, if (c == 1) ny else c - 1),   # izquierda
    c(f, if (c == ny) 1 else c + 1)    # derecha
  )
}

#' Ejecuta un paso de crecimiento de la colonia (asíncrono, orden aleatorio)
#'
#' @param estado matriz 0/1 de ocupación
#' @param resist matriz de resistencia por celda (NA si vacío)
#' @param N campo de nutriente actual
#' @param N_umbral nutriente mínimo para que una celda sea "proliferativa"
#'        (por debajo => fenotipo "dormante", no divide)
#' @param p_max probabilidad máxima de división por paso (nutriente saturante)
#' @param N_half constante de saturación de la probabilidad de división
#'        (misma forma funcional que el consumo, ligada al flujo local)
#' @param mutacion_sd sd del ruido gaussiano aditivo sobre `resist` heredado
#'
#' @return lista con `estado`, `resist` actualizados y `fenotipo` derivado
crecer_colonia <- function(estado, resist, N,
                            N_umbral   = 0.15,
                            p_max      = 0.6,
                            N_half     = 0.3,
                            mutacion_sd = 0.02) {

  nx <- nrow(estado); ny <- ncol(estado)
  ocupados <- which(estado == 1L, arr.ind = TRUE)

  # Fenotipo derivado del campo local (NO es un rasgo heredado)
  fenotipo <- matrix(NA_character_, nx, ny)
  if (nrow(ocupados) > 0) {
    idx <- ocupados
    n_local <- N[idx]
    fenotipo[idx] <- ifelse(n_local >= N_umbral, "proliferativo", "dormante")
  }

  if (nrow(ocupados) == 0) {
    return(list(estado = estado, resist = resist, fenotipo = fenotipo))
  }

  # Orden aleatorio de actualización (evita sesgo de barrido tipo raster)
  orden <- sample(seq_len(nrow(ocupados)))

  for (i in orden) {
    f <- ocupados[i, 1]; c <- ocupados[i, 2]

    n_local <- N[f, c]
    if (n_local < N_umbral) next   # dormante: no intenta dividir

    vecinos <- vecinos_von_neumann(f, c, nx, ny)
    libres  <- vecinos[estado[vecinos] == 0L, , drop = FALSE]
    if (nrow(libres) == 0) next    # sin espacio local: no divide

    p_div <- p_max * (n_local / (n_local + N_half))
    if (runif(1) > p_div) next

    destino <- libres[sample.int(nrow(libres), 1), ]

    # Herencia con mutación del rasgo de resistencia
    hijo_resist <- resist[f, c] + rnorm(1, 0, mutacion_sd)
    hijo_resist <- min(max(hijo_resist, 0), 1)   # truncado a [0,1]

    estado[destino[1], destino[2]] <- 1L
    resist[destino[1], destino[2]] <- hijo_resist
  }

  list(estado = estado, resist = resist, fenotipo = fenotipo)
}
