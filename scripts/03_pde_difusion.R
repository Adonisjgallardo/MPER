## ============================================================
## 03_pde_difusion.R
## Componente PDE / reacción-difusión: dinámica del nutriente
##
##   dN/dt = D * Laplaciano(N)  -  k * estado * N/(N + N_half)
##
## Difusión: 2do orden centrado (FTCS), fronteras periódicas.
## Consumo: cinética tipo Michaelis-Menten (evita N < 0 y
##          reproduce saturación del consumo a alta concentración).
## ============================================================

#' Un paso de integración explícita del campo de nutriente
#'
#' @param N matriz de concentración de nutriente
#' @param estado matriz 0/1 de ocupación (consume nutriente donde ==1)
#' @param D coeficiente de difusión
#' @param k tasa máxima de consumo por sitio ocupado
#' @param dt paso de tiempo
#' @param dx paso espacial de la retícula
#' @param N_half constante de semisaturación (Michaelis-Menten)
actualizar_nutriente <- function(N, estado, D, k, dt,
                                  dx = 1, N_half = 0.3) {
  lap     <- laplaciano_periodico(N, dx)
  consumo <- k * estado * (N / (N + N_half))
  N_new   <- N + dt * (D * lap - consumo)
  pmax(N_new, 0)
}

#' Verifica la condición de estabilidad de von Neumann para FTCS 2D
#' (condición necesaria, no analiza el término de reacción)
#'
#' @return TRUE/FALSE + mensaje informativo
verificar_estabilidad_dt <- function(D, dt, dx = 1) {
  criterio <- D * dt / dx^2
  ok <- criterio <= 0.25
  if (!ok) {
    warning(sprintf(
      "Posible inestabilidad numérica: D*dt/dx^2 = %.4f > 0.25. Reduce dt o aumenta dx.",
      criterio))
  }
  invisible(ok)
}
