## ============================================================
## 01_campos.R nb
## Inicialización de campos continuos sobre la retícula (HCA)
## Modelo nutriente-limitado puro (antibiótico se añade después)
## Condiciones de frontera: PERIÓDICAS (toroide)
## ============================================================

#' Crea el campo inicial de nutriente
#'
#' @param nx,ny dimensiones de la retícula
#' @param N0 concentración inicial homogénea de nutriente
#' @param ruido_sd ruido espacial opcional (0 = campo perfectamente homogéneo)
crear_campo_nutriente <- function(nx, ny, N0 = 1.0, ruido_sd = 0) {
  N <- matrix(N0, nrow = nx, ncol = ny)
  if (ruido_sd > 0) {
    N <- N + matrix(rnorm(nx * ny, 0, ruido_sd), nx, ny)
    N <- pmax(N, 0)
  }
  N
}

#' Laplaciano discreto 2D con condiciones de frontera PERIÓDICAS
#' (vectorizado, sin loops, usando desplazamiento circular de índices)
laplaciano_periodico <- function(M, dx = 1) {
  nx <- nrow(M); ny <- ncol(M)
  arriba  <- M[c(2:nx, 1), , drop = FALSE]
  abajo   <- M[c(nx, 1:(nx - 1)), , drop = FALSE]
  derecha <- M[, c(2:ny, 1), drop = FALSE]
  izq     <- M[, c(ny, 1:(ny - 1)), drop = FALSE]
  (arriba + abajo + derecha + izq - 4 * M) / dx^2
}

# ------------------------------------------------------------
# Parámetros globales de la retícula (ajustar según necesidad)
# ------------------------------------------------------------
PARAMS_RETICULA <- list(
  nx = 100,
  ny = 100,
  dx = 1
)

