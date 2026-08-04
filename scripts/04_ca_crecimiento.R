## ============================================================
## 04_ca_crecimiento.R
## Componente CA: división Y MUERTE celular, ahora gobernadas por
## selección ESTABILIZADORA sobre un fenotipo cuantitativo z=g+e,
## con óptimo local theta(x,y,t) que depende del antibiótico.
##
## Cambios respecto a la versión nutriente-limitada pura:
##
## 1) SELECCIÓN ESTABILIZADORA (Caja 3): la aptitud de un individuo
##    ya no depende solo del nutriente local, sino de qué tan cerca
##    está su fenotipo z del óptimo local theta(x,y,t):
##
##       w(z) = exp( -(theta - z)^2 / (2*omega^2) )
##
##    omega^2 controla el ANCHO de la selección (1/omega^2 es la
##    "fuerza" de la selección estabilizadora, tal como en la Caja 3).
##
## 2) MUERTE (nueva): antes, una bacteria sin espacio o sin nutriente
##    simplemente no se dividía, pero JAMÁS moría → la población
##    nunca podía declinar. El "rescate evolutivo" solo tiene sentido
##    si existe la posibilidad real de declive/extinción. Se añade
##    una probabilidad de muerte que aumenta con la concentración
##    LOCAL de antibiótico y con qué tan lejos está z del óptimo:
##
##       p_mort = mort_base + mort_estres * A/(A+K_mort) * (1 - w(z))
##
##    Sin antibiótico (A=0) la mortalidad es solo la basal (turnover
##    normal); con antibiótico, los fenotipos mal adaptados mueren
##    mucho más rápido que los bien adaptados.
##
## 3) HERENCIA cuantitativa: al dividir, el hijo hereda `g` (valor
##    genético aditivo) del padre + mutación gaussiana (difusión en
##    espacio de rasgos, como se derivó antes). La componente `e`
##    del hijo se sortea de nuevo, independiente — no se hereda.
##
## Vecindad: von Neumann (4 vecinos), fronteras PERIÓDICAS.
## ============================================================

#' Índices de los 4 vecinos von Neumann con envoltura periódica (toroide)
vecinos_von_neumann <- function(f, c, nx, ny) {
  rbind(
    c(if (f == 1) nx else f - 1, c),
    c(if (f == nx) 1 else f + 1, c),
    c(f, if (c == 1) ny else c - 1),
    c(f, if (c == ny) 1 else c + 1)
  )
}

#' Fitness de selección estabilizadora individual (Caja 3, sin
#' integrar sobre la distribución poblacional: aquí se evalúa
#' directamente en el fenotipo REALIZADO z de cada agente)
#'
#' @param theta_local óptimo fenotípico local (matriz, misma forma que z)
#' @param z fenotipo realizado local (matriz)
#' @param omega2 ancho^2 de la selección estabilizadora
fitness_estabilizador <- function(theta_local, z, omega2) {
  exp(-(theta_local - z)^2 / (2 * omega2))
}

#' Ejecuta un paso de crecimiento + muerte de la colonia (asíncrono,
#' orden aleatorio), con selección estabilizadora sobre z=g+e.
#'
#' @param estado matriz 0/1 de ocupación
#' @param g matriz de valor genético aditivo (heredable)
#' @param e matriz de componente ambiental (ya resorteada este paso)
#' @param N campo de nutriente
#' @param A campo de antibiótico (0 = sin estrés)
#' @param theta0,theta1 óptimo fenotípico sin/con antibiótico saturante
#' @param omega2 ancho^2 de la selección estabilizadora
#' @param K_theta semisaturación del mapeo A -> theta_local
#' @param N_umbral,p_max,N_half igual que en la versión nutriente-pura
#' @param mort_base mortalidad basal (turnover, sin relación al antibiótico)
#' @param mort_estres mortalidad adicional máxima por desajuste fenotípico
#'        bajo antibiótico saturante
#' @param K_mort semisaturación de A para la mortalidad por estrés
#' @param mutacion_sd sd de la mutación aditiva sobre `g` por división
#' @param sigma_e sd de la componente ambiental de los hijos
#'
#' @return lista con `estado`,`g`,`e`,`fenotipo` (proliferativo/dormante,
#'         derivado del nutriente) y contadores de nacimientos/muertes
crecer_colonia <- function(estado, g, e, N, A,
                            theta0 = 0.1, theta1 = 0.8, omega2 = 0.05,
                            K_theta = 0.3,
                            N_umbral = 0.15, p_max = 0.6, N_half = 0.3,
                            mort_base = 0.01, mort_estres = 0.9, K_mort = 0.3,
                            mutacion_sd = 0.02, sigma_e = 0.05) {

  nx <- nrow(estado); ny <- ncol(estado)
  z <- fenotipo_z(g, e)
  theta_local <- optimo_local(A, theta0, theta1, K_theta)
  w <- fitness_estabilizador(theta_local, z, omega2)   # NA fuera de la colonia

  # Fenotipo derivado del campo local de NUTRIENTE (igual que antes)
  fenotipo <- matrix(NA_character_, nx, ny)
  ocupados <- which(estado == 1L, arr.ind = TRUE)
  if (nrow(ocupados) > 0) {
    n_local <- N[ocupados]
    fenotipo[ocupados] <- ifelse(n_local >= N_umbral, "proliferativo", "dormante")
  }

  n_nacimientos <- 0L; n_muertes <- 0L

  if (nrow(ocupados) == 0) {
    return(list(estado = estado, g = g, e = e, fenotipo = fenotipo,
                nacimientos = n_nacimientos, muertes = n_muertes))
  }

  orden <- sample(seq_len(nrow(ocupados)))

  for (i in orden) {
    f <- ocupados[i, 1]; c <- ocupados[i, 2]
    if (estado[f, c] == 0L) next   # pudo haber muerto ya en este mismo paso

    # ---- 1) ¿Muere este paso? ----------------------------------
    p_mort <- mort_base + mort_estres * (A[f, c] / (A[f, c] + K_mort)) * (1 - w[f, c])
    if (runif(1) < p_mort) {
      estado[f, c] <- 0L
      g[f, c] <- NA_real_
      e[f, c] <- NA_real_
      n_muertes <- n_muertes + 1L
      next   # un sitio muerto no divide en el mismo paso
    }

    # ---- 2) ¿Divide este paso? ---------------------------------
    n_local <- N[f, c]
    if (n_local < N_umbral) next   # dormante (limitado por nutriente): no divide

    vecinos <- vecinos_von_neumann(f, c, nx, ny)
    libres  <- vecinos[estado[vecinos] == 0L, , drop = FALSE]
    if (nrow(libres) == 0) next

    p_div <- p_max * (n_local / (n_local + N_half)) * w[f, c]
    if (runif(1) > p_div) next

    destino <- libres[sample.int(nrow(libres), 1), ]

    # Herencia: g muta (difusión en espacio de rasgos), e se resortea
    hijo_g <- g[f, c] + rnorm(1, 0, mutacion_sd)
    hijo_g <- min(max(hijo_g, 0), 1)
    hijo_e <- rnorm(1, 0, sigma_e)

    estado[destino[1], destino[2]] <- 1L
    g[destino[1], destino[2]] <- hijo_g
    e[destino[1], destino[2]] <- hijo_e
    n_nacimientos <- n_nacimientos + 1L
  }

  list(estado = estado, g = g, e = e, fenotipo = fenotipo,
       nacimientos = n_nacimientos, muertes = n_muertes)
}
