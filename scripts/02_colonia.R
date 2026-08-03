## ============================================================
## 02_colonia.R
## Siembra de la colonia inicial sobre la retícula (HCA)
##
## A diferencia del modelo agente-continuo anterior, aquí NO se
## generan miles de individuos: la colonia nace de un puñado de
## sitios ocupados (típicamente 1, el "inóculo" central) y crece
## por división hacia sitios vecinos vacíos.
##
## Cada sitio ocupado lleva un único rasgo por-celda:
##   resist: resistencia al antibiótico, en [0,1]
## que se HEREDA con mutación en cada evento de división
## (preparado para la capa de antibiótico que se añadirá después;
## en esta fase nutriente-limitada pura, `resist` es inerte:
## se transmite y muta, pero no afecta todavía al crecimiento).
## ============================================================

#' Crea el estado inicial de la colonia (retícula vacía + inóculo)
#'
#' @param nx,ny dimensiones de la retícula (deben coincidir con el campo)
#' @param inoculo lista de coordenadas (fila,col) iniciales ocupadas;
#'        por defecto, un único sitio en el centro
#' @param resist_inicial valor (o función) de resistencia de fundadores
crear_colonia_inicial <- function(nx, ny,
                                   inoculo = NULL,
                                   resist_inicial = 0.1,
                                   seed = 123) {
  set.seed(seed)

  estado <- matrix(0L, nx, ny)                 # 0 = vacío, 1 = ocupado
  resist <- matrix(NA_real_, nx, ny)            # rasgo por celda (solo si ocupado)

  if (is.null(inoculo)) {
    inoculo <- matrix(c(ceiling(nx / 2), ceiling(ny / 2)), ncol = 2)
  }

  for (i in seq_len(nrow(inoculo))) {
    f <- inoculo[i, 1]; c <- inoculo[i, 2]
    estado[f, c] <- 1L
    resist[f, c] <- resist_inicial
  }

  list(estado = estado, resist = resist)
}
