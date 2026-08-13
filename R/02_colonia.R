## ============================================================
## 02_colonia.R
## Siembra de la colonia inicial (HCA) — versión genética
## cuantitativa (Caja 3 del documento de rescate evolutivo)
##
## El fenotipo de resistencia ya NO es un único número heredado
## directamente. Se descompone como:
##
##     z = g + e
##
##   g: valor genético aditivo (heredable, muta en cada división)
##   e: desviación ambiental (NO heredable, se resortea cada paso
##      de tiempo para TODOS los individuos vivos — ruido de
##      desarrollo/microambiente, independiente entre generaciones)
##   z: fenotipo realizado (= g+e), sobre el que actúa la SELECCIÓN
##
## Esta separación es la que permite luego calcular
## sigma_g^2 (varianza genética aditiva) y sigma_e^2 (varianza
## ambiental) por separado, que es justo lo que exige la Caja 3
## (phi = sigma_g^2 / V, con V = sigma_g^2 + sigma_e^2 + omega^2)
## para comparar la trayectoria simulada contra la predicción
## analítica de genética cuantitativa.
## ============================================================

#' Crea el estado inicial de la colonia (retícula vacía + inóculo)
#'
#' @param nx,ny dimensiones de la retícula
#' @param inoculo coordenadas (fila,col) iniciales ocupadas;
#'        por defecto, un único sitio en el centro
#' @param g_inicial valor genético aditivo de los fundadores
#'        (ANTES del cambio ambiental: representa el fenotipo
#'        "silvestre" bien adaptado al entorno sin antibiótico)
#' @param sigma_e sd de la componente ambiental del fenotipo
crear_colonia_inicial <- function(nx, ny,
                                   inoculo = NULL,
                                   g_inicial = 0.1,
                                   sigma_e = 0.05,
                                   seed = 123) {
  set.seed(seed)

  estado <- matrix(FALSE, nx, ny)      # FALSE = vacío, TRUE = ocupado
  g <- matrix(NA_real_, nx, ny)       # componente genética aditiva (heredable)
  e <- matrix(NA_real_, nx, ny)       # componente ambiental (no heredable)

  if (is.null(inoculo)) {
    inoculo <- matrix(c(ceiling(nx / 2), ceiling(ny / 2)), ncol = 2)
  }

  for (i in seq_len(nrow(inoculo))) {
    f <- inoculo[i, 1]; c <- inoculo[i, 2]
    estado[f, c] <- TRUE
    g[f, c] <- g_inicial
    e[f, c] <- rnorm(1, 0, sigma_e)
  }

  list(estado = estado, g = g, e = e)
}

#' Resortea la componente ambiental `e` de TODOS los sitios ocupados
#' en el paso actual. Modela el supuesto de genética cuantitativa
#' de que e ~ N(0, sigma_e^2) i.i.d. en cada generación, independiente
#' del valor genético y del pasado del individuo.
#'
#' @param estado matriz lógica de ocupación
#' @param e matriz de componente ambiental (se sobreescribe donde
#'        estado == TRUE; se deja NA donde está vacío)
#' @param sigma_e sd de la componente ambiental
resortear_ambiente <- function(estado, e, sigma_e) {
  idx <- which(estado == TRUE)
  if (length(idx) > 0) e[idx] <- rnorm(length(idx), 0, sigma_e)
  e
}

#' Fenotipo realizado z = g + e (solo definido donde estado==1)
fenotipo_z <- function(g, e) g + e
