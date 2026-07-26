#Archivo 02_bacterias.R

#Ahora agregamos las bacterias, pero todavía sin antibiótico y sin selección fuerte. 
# La idea es construir la condición inicial: una comunidad diversa que existe antes de cualquier perturbación.

#La pregunta que responde esta etapa es:

#¿Existe una reserva fenotípica distribuida en un paisaje heterogéneo antes del disturbio?

library(data.table)

#Paisaje
# E(x,y) 

paisaje <- readRDS("data/paisaje_inicial.rds") |>
  as.data.table()

# Y Comunidad 
# P=(F1,F2,F3,F4:7)

crear_poblacion_inicial <- function(pop, n = 10000, seed = 123) {
  set.seed(seed)
  pop[, F1 := rbeta(n, 2, 2)]                   # Reproducción
  pop[, F2 := rgamma(n, shape = 2, rate = 2)]      # dispersión
  pop[, F3 := rlnorm(n, meanlog = -2, sdlog = 0.5)]  # mutabilidad
  pop[, F4 := plogis(-1 + 1.5 * E1 + rnorm(n, 0, 0.3))]  # resistencia ~ E_antibiótico (¡confirmar cuál E!)
  pop[, F5 := plogis(-1 + 1.5 * E2 + rnorm(n, 0, 0.3))]  # resistencia ~ E_antibiótico 
  pop[, F6 := plogis(-1 + 1.5 * E3 + rnorm(n, 0, 0.3))]  # resistencia ~ E_antibiótico 
  pop[, F7 := plogis(-1 + 1.5 * E4 + rnorm(n, 0, 0.3))]  # resistencia ~ E_antibiótico 
  pop[, `:=`(generacion = 1L, viva = TRUE)] # datos del agente
  pop
}

pop <- crear_poblacion_inicial(paisaje)

rm(list=ls()[-which(ls() %in% "pop")])

# Verificación rápida (no-clones + distribución de fenotipos)

verificar_poblacion <- function(pop) {
  n_unicos <- nrow(unique(pop[, .(F1, F2, F3, F4)]))
  cat("Individuos totales:      ", nrow(pop), "\n")
  cat("Fenotipos únicos (no clon):", n_unicos, "\n")
  cat("\nResumen de fenotipos:\n")
  print(summary(pop[, .(F1, F2, F3, F4)]))
  invisible(pop)
}

verificar_poblacion(pop)

saveRDS(pop, "data/bacterias.rds")
