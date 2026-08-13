## ============================================================
## 07_genetica_cuantitativa.R
## Fórmulas ANALÍTICAS del documento de rescate evolutivo
## (Caja 3: genética cuantitativa: Caja 4: genética de
## poblaciones/probabilidad de rescate; Suplemento S2: tiempos
## de persistencia vía proceso de nacimiento-muerte), más
## utilidades para comparar contra la simulación HCA.
## ============================================================

## ------------------------------------------------------------
## CAJA 3 — Genética cuantitativa: trayectoria analítica de N_t
## bajo un CAMBIO AMBIENTAL ABRUPTO (theta pasa de theta0 a theta1)
## ------------------------------------------------------------

#' Trayectoria analítica de N_t tras un cambio ambiental abrupto
#' (fórmula de la sección "Abrupt environmental change" de la Caja 3,
#' evolución determinista, demografía densidad-independiente)
#'
#' @param t vector de generaciones (t=0,1,2,...)
#' @param N0 tamaño poblacional inicial (en el momento del shock)
#' @param Wmax aptitud máxima posible
#' @param omega2 ancho^2 de la selección estabilizadora
#' @param sigma_g2 varianza genética aditiva (constante, supuesto de la Caja 3)
#' @param sigma_e2 varianza ambiental
#' @param d0 desviación inicial del óptimo: d0 = theta1 - g_bar_0
#'        (en el momento del shock, ANTES de que la población responda)
#' @return vector N_t (misma longitud que `t`)
Nt_analitico_qg <- function(t, N0, Wmax, omega2, sigma_g2, sigma_e2, d0) {
  V   <- sigma_g2 + sigma_e2 + omega2
  phi <- sigma_g2 / V
  N0 * (Wmax * sqrt(omega2 / V))^t *
    exp(-(d0^2) / (2 * V) * (1 - (1 - phi)^(2 * t)) / (1 - (1 - phi)^2))
}

#' Trayectoria analítica del valor genético medio E[g_bar_t]
#' (ecuación B3.2, recursión de Lande con varianza genética constante).
#' Con k=0 (óptimo estático tras el shock) se reduce al acercamiento
#' geométrico de g_bar al óptimo.
#'
#' @param t vector de generaciones
#' @param k tasa de cambio LINEAL del óptimo (0 para shock abrupto
#'          seguido de óptimo estático; >0 para "cambio gradual")
#' @param sigma_g2,sigma_e2,omega2 igual que arriba
#' @param g_bar0 valor genético medio inicial
gbar_analitico_qg <- function(t, k, sigma_g2, sigma_e2, omega2, g_bar0) {
  V   <- sigma_g2 + sigma_e2 + omega2
  phi <- sigma_g2 / V
  k * t - (k / phi) * (1 - (1 - phi)^t) + (1 - phi)^t * g_bar0
}

#' Tasa CRÍTICA de cambio ambiental (ecuación B3.4): la máxima
#' velocidad k a la que el óptimo puede moverse manteniendo aptitud
#' media >= reemplazo en estado estacionario. Útil como referencia
#' aunque tu escenario sea de shock ABRUPTO, no gradual: sirve para
#' saber qué tan "brusco" es tu shock comparado con el límite adaptativo.
tasa_critica_qg <- function(Wmax, omega2, sigma_g2, sigma_e2) {
  V   <- sigma_g2 + sigma_e2 + omega2
  phi <- sigma_g2 / V
  phi * sqrt(2 * V * log(Wmax * sqrt(omega2 / V)))
}

## ------------------------------------------------------------
## CAJA 4 — Genética de poblaciones: probabilidad de rescate
## ------------------------------------------------------------

#' Probabilidad de rescate vía mutación DE NOVO (ecuación B4.1),
#' asumiendo p_est constante en el tiempo.
#'
#' @param lambda_DN tasa esperada de linajes mutantes exitosos
#'        (= integral de mu(t)*N_silvestre(t) en el tiempo; ver
#'        `lambda_DN_exponencial()` para el caso de declive exponencial)
#' @param p_est probabilidad de establecimiento de un linaje mutante
prob_rescate_de_novo <- function(lambda_DN, p_est) {
  1 - exp(-lambda_DN * p_est)
}

#' Caso particular de lambda_DN cuando la población silvestre declina
#' exponencialmente: N_w(t) = N0*exp(-r*t), con tasa de mutación mu
#' constante por individuo. lambda_DN = mu*N0/r (ver Caja 4).
lambda_DN_exponencial <- function(mu, N0, r) mu * N0 / r

#' Probabilidad TOTAL de rescate (ecuación B4.2): combina mutantes
#' preexistentes (varianza genética "standing", lambda_SGV,
#' aproximado Poisson) + mutantes de novo tras el shock.
prob_rescate_total <- function(lambda_SGV, lambda_DN, p_est) {
  1 - exp(-(lambda_SGV + lambda_DN) * p_est)
}

## ------------------------------------------------------------
## SUPLEMENTO S2 — Tiempos de persistencia (proceso de
## nacimiento-muerte / branching process)
## ------------------------------------------------------------

#' P(T_ext <= t) para un proceso de nacimiento-muerte puro con
#' tasa de nacimiento `lambda`, tasa de muerte `mu`, y n0 individuos
#' iniciales (ecuación S3/S6, forma general con n0 arbitrario).
prob_extincion_hasta_t <- function(t, lambda, mu, n0) {
  factor <- (exp(-(lambda - mu) * t) - 1) / (exp(-(lambda - mu) * t) - lambda / mu)
  factor^n0
}

#' P(T_ext <= t | CON evolución): mezcla de m0 mutantes (lambda_m,mu_m)
#' y (N0-m0) silvestres (lambda_w,mu_w), sin nuevas mutaciones
#' (ecuación S5).
prob_extincion_con_evolucion <- function(t, N0, m0, lambda_w, mu_w, lambda_m, mu_m) {
  p_w <- prob_extincion_hasta_t(t, lambda_w, mu_w, N0 - m0)
  p_m <- prob_extincion_hasta_t(t, lambda_m, mu_m, m0)
  p_w * p_m
}

#' Probabilidad de EXCESO de persistencia por evolucionar vs. no
#' evolucionar (ecuación S9): diferencia entre sobrevivir-con-evolución
#' y sobrevivir-sin-evolución, en función del tiempo.
prob_exceso_persistencia <- function(t, N0, m0, lambda_w, mu_w, lambda_m, mu_m) {
  P_evol <- 1 - prob_extincion_con_evolucion(t, N0, m0, lambda_w, mu_w, lambda_m, mu_m)
  P_no_evol <- 1 - prob_extincion_hasta_t(t, lambda_w, mu_w, N0)
  P_evol - P_no_evol
}

#' Probabilidad ASINTÓTICA de rescate (t->Inf), ecuación S11:
#' solo depende de si el linaje mutante tiene lambda_m > mu_m.
prob_rescate_asintotico <- function(lambda_m, mu_m, m0) {
  if (lambda_m <= mu_m) return(0)
  (mu_m / lambda_m)^m0
}

## ------------------------------------------------------------
## Extracción EMPÍRICA de tiempos de persistencia desde réplicas
## de la simulación HCA (para comparar contra las fórmulas de arriba)
## ------------------------------------------------------------

#' A partir de una lista de resúmenes (de `simular_hca_replicas()`),
#' extrae el tiempo (en pasos DESDE la introducción del antibiótico)
#' en que cada réplica cruza por primera vez el umbral crítico Nc.
#' Réplicas que nunca cruzan el umbral devuelven `Inf` (censura por
#' la derecha: "sobrevivió todo el horizonte simulado").
#'
#' @param replicas lista devuelta por `simular_hca_replicas()`
#' @param Nc umbral crítico (si es NULL, usa el guardado en cada corrida)
#' @return data.frame con columnas replica, tiempo_persistencia, extinguio
tiempos_persistencia_empiricos <- function(replicas, Nc = NULL) {
  do.call(rbind, lapply(seq_along(replicas), function(i) {
    r <- replicas[[i]]
    intro <- r$paso_introduccion
    if (is.na(r$paso_extincion)) {
      t_persist <- max(r$resumen$paso) - intro
      extinguio <- FALSE
    } else {
      t_persist <- r$paso_extincion - intro
      extinguio <- TRUE
    }
    data.frame(replica = i, tiempo_persistencia = t_persist, extinguio = extinguio)
  }))
}

#' Fracción de réplicas RESCATADAS (no cruzaron Nc dentro del horizonte
#' simulado) — el análogo empírico directo de `prob_rescate_asintotico()`
#' y `prob_rescate_total()`.
fraccion_rescatada_empirica <- function(replicas) {
  tp <- tiempos_persistencia_empiricos(replicas)
  mean(!tp$extinguio)
}
