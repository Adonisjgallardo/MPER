## ============================================================
## 03_pde_difusion.R
## Componente PDE / reacción-difusión: nutriente Y antibiótico
##
## NUTRIENTE:
##   dN/dt = D_N * Laplaciano(N)  -  k * estado * N/(N + N_half)
##
## ANTIBIÓTICO (nuevo):
##   dA/dt = D_A * Laplaciano(A)  -  delta_A * A  + relajación hacia A_objetivo
##   con un CAMBIO AMBIENTAL ABRUPTO: A es idénticamente 0 hasta
##   el paso `paso_introduccion` (Box 3 del documento de rescate
##   evolutivo: cambio ambiental abrupto como la introducción de
##   un contaminante o antibiótico), y heterogéneo tanto en
##   espacio (gradiente inicial) como en tiempo (modulación
##   post-introducción: pulsos / fluctuación sinusoidal).
##
## Difusión: 2do orden centrado (FTCS), fronteras periódicas.
## ============================================================

#' Un paso de integración explícita del campo de nutriente
## (sin cambios respecto a la versión anterior)
actualizar_nutriente <- function(N, estado, D, k, dt,
                                  dx = 1, N_half = 0.3) {
  lap     <- laplaciano_periodico(N, dx)
  consumo <- k * estado * (N / (N + N_half))
  N_new   <- N + dt * (D * lap - consumo)
  pmax(N_new, 0)
}

#' Verifica la condición de estabilidad de von Neumann para FTCS 2D
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

## ------------------------------------------------------------
## Campo de antibiótico: heterogeneidad ESPACIAL inicial
## ------------------------------------------------------------

#' Crea la "plantilla" espacial del antibiótico (antes de escalar
#' por tiempo/magnitud). Inspirada en los diseños de gradiente en
#' placa de la literatura ABM bacteria-antibiótico (ej. gradiente
#' creciente desde 1xMIC hasta 1000xMIC en un extremo de la placa).
#'
#' @param nx,ny dimensiones de la retícula
#' @param tipo "uniforme", "lineal" (gradiente en x) o "radial"
#'        (mínimo en el centro, útil si el inóculo también está
#'        ahí y se quiere estrés creciente hacia afuera)
#' @param foco coordenadas (fila,col) del foco para tipo="radial";
#'        por defecto, el centro de la retícula
crear_plantilla_antibiotico <- function(nx, ny,
                                         tipo = c("uniforme", "lineal", "radial"),
                                         foco = NULL) {
  tipo <- match.arg(tipo)

  if (tipo == "uniforme") {
    plantilla <- matrix(1, nx, ny)
  } else if (tipo == "lineal") {
    grad_x <- seq(0, 1, length.out = nx)
    plantilla <- matrix(rep(grad_x, ny), nrow = nx, ncol = ny)
  } else { # radial
    if (is.null(foco)) foco <- c(ceiling(nx / 2), ceiling(ny / 2))
    fi <- matrix(1:nx, nx, ny)
    fj <- matrix(1:ny, nx, ny, byrow = TRUE)
    dist_foco <- sqrt((fi - foco[1])^2 + (fj - foco[2])^2)
    plantilla <- dist_foco / max(dist_foco)   # 0 en el foco, 1 en el borde
  }
  plantilla
}

#' Concentración de antibiótico "objetivo" en el paso actual: 0
#' antes de `paso_introduccion` (CAMBIO ABRUPTO), y la plantilla
#' espacial escalada por A_max y por una modulación temporal
#' después (heterogeneidad temporal: fluctuación).
#'
#' @param plantilla matriz espacial en `[0,1]` (de `crear_plantilla_antibiotico`)
#' @param paso paso de tiempo actual
#' @param paso_introduccion paso en el que ocurre el CAMBIO ABRUPTO
#' @param A_max concentración máxima post-introducción
#' @param modo_temporal "constante", "pulsos" o "sinusoidal"
#' @param periodo periodo (en pasos) de la fluctuación (pulsos/sinusoidal)
antibiotico_objetivo <- function(plantilla, paso, paso_introduccion,
                                  A_max = 1.0,
                                  modo_temporal = c("constante", "pulsos", "sinusoidal"),
                                  periodo = 40) {
  modo_temporal <- match.arg(modo_temporal)

  if (paso < paso_introduccion) {
    return(matrix(0, nrow(plantilla), ncol(plantilla)))  # SIN antibiótico: cambio abrupto
  }

  t_desde_intro <- paso - paso_introduccion
  factor_t <- switch(modo_temporal,
    constante  = 1,
    pulsos     = as.numeric((t_desde_intro %% periodo) < (periodo / 2)),
    sinusoidal = 0.5 * (1 + sin(2 * pi * t_desde_intro / periodo))
  )

  A_max * factor_t * plantilla
}

#' Un paso de integración del campo de antibiótico: difusión +
#' degradación (natural/enzimática de fondo) + relajación hacia
#' el objetivo espacio-temporal (dosificación externa).
#'
#' Se modela como relajación hacia `A_objetivo` en vez de una
#' fuente/sumidero rígida, para que el campo siga siendo suave
#' (difusivo) aun cuando el objetivo cambie abruptamente en t.
#'
#' @param A campo actual de antibiótico
#' @param A_objetivo campo "deseado" en este paso (de `antibiotico_objetivo`)
#' @param D_A coeficiente de difusión del antibiótico
#' @param delta_A tasa de degradación natural
#' @param tasa_dosificacion qué tan rápido A persigue a A_objetivo
#'        (1 = instantáneo, valores menores = dosificación gradual)
actualizar_antibiotico <- function(A, A_objetivo, D_A, delta_A, dt,
                                    dx = 1, tasa_dosificacion = 0.3) {
  lap <- laplaciano_periodico(A, dx)
  fuente <- tasa_dosificacion * (A_objetivo - A)
  A_new <- A + dt * (D_A * lap - delta_A * A + fuente)
  pmax(A_new, 0)
}

## ------------------------------------------------------------
## Puente genética-cuantitativa <-> campo de antibiótico
## ------------------------------------------------------------

#' Mapea la concentración LOCAL de antibiótico A(x,y,t) al óptimo
#' fenotípico LOCAL theta(x,y,t) de la Caja 3 (genética cuantitativa).
#'
#' Sin antibiótico (A=0) el óptimo es theta0 (el fenotipo silvestre,
#' de bajo costo, es el mejor adaptado). A medida que A crece, el
#' óptimo se desplaza hacia theta1 (alta resistencia), con forma
#' Michaelis-Menten (saturante): el desplazamiento del óptimo no
#' sigue creciendo indefinidamente con dosis cada vez mayores.
#'
#' Esta es la pieza que conecta el "cambio ambiental abrupto" de
#' la Caja 3 (theta salta de theta0 a theta1) con la heterogeneidad
#' ESPACIAL: como A varía en el espacio (gradiente/plantilla),
#' theta también varía en el espacio, dando un óptimo fenotípico
#' que depende de la posición, no solo del tiempo.
#'
#' @param A matriz de concentración local de antibiótico
#' @param theta0 óptimo fenotípico SIN antibiótico (antes del shock,
#'        o donde A~0 localmente)
#' @param theta1 óptimo fenotípico bajo antibiótico saturante
#' @param K_theta constante de semisaturación (misma forma que
#'        Michaelis-Menten del consumo de nutriente)
#' @return matriz theta(x,y) del mismo tamaño que A
optimo_local <- function(A, theta0, theta1, K_theta = 0.3) {
  theta0 + (theta1 - theta0) * (A / (A + K_theta))
}
