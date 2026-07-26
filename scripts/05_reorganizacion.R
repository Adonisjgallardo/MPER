#_reorganizacion.R

#Esta es la prueba conceptual

#Si después del disturbio observamos aumento de ciertos fenotipos, podremos decir:

#El ambiente no generó la diversidad; reorganizó una diversidad existente.

###############################################################
# MPER 1.0
# Modelo de Perturbación Ecológica y Reorganización Adaptativa
#
# Archivo 05:
# Evolución después del disturbio
###############################################################

rm(list=ls())

library(ggplot2)
library(dplyr)


set.seed(123)


###############################################################
# Cargar estado perturbado
###############################################################

paisaje <- readRDS(
  "paisaje_perturbado.rds"
)


bacterias <- readRDS(
  "bacterias_pre_disturbio.rds"
)



###############################################################
# Función fitness
###############################################################

calcular_fitness <- function(bacterias,paisaje){

  
  datos <- bacterias %>%
    left_join(
      paisaje,
      by=c("x","y")
    )


  distancia <- sqrt(
    (datos$F1-datos$E1)^2+
    (datos$F2-datos$E2)^2+
    (datos$F3-datos$E3)^2+
    (datos$F4-datos$E4)^2
  )


  datos$fitness <- exp(-distancia)


  return(datos)

}



###############################################################
# Número de generaciones después del disturbio
###############################################################

generaciones <- 100


historial_post <- list()



###############################################################
# Evolución post-disturbio
###############################################################

for(g in 1:generaciones){

  
  bacterias <-
    calcular_fitness(
      bacterias,
      paisaje
    )


  bacterias$generacion <- g

  
  historial_post[[g]] <- bacterias

  #############################################################
  # Selección
  #############################################################

  sobrevivientes <-
    bacterias %>%
    filter(
      runif(n()) < fitness
    )


  if(nrow(sobrevivientes)<100){

    sobrevivientes <-
      bacterias %>%
      slice_sample(
        n=100
      )

  }



  #############################################################
  # Reproducción
  #############################################################

  bacterias <-
    sobrevivientes %>%
    slice_sample(
      n=10000,
      replace=TRUE
    )


  bacterias$id <-
    1:nrow(bacterias)



  #############################################################
  # Mutación
  #############################################################

  bacterias$F1 <-
    bacterias$F1+
    rnorm(10000,0,0.05)


  bacterias$F2 <-
    bacterias$F2+
    rnorm(10000,0,0.05)


  bacterias$F3 <-
    bacterias$F3+
    rnorm(10000,0,0.05)


  bacterias$F4 <-
    bacterias$F4+
    rnorm(10000,0,0.05)


}



###############################################################
# Guardar evolución
###############################################################

saveRDS(
  historial_post,
  "historial_post_disturbio.rds"
)



###############################################################
# Comparación fenotípica
###############################################################

final <-
  historial_post[[generaciones]]



ggplot(
  final,
  aes(
    F1,
    F2
  )
)+

geom_point(
  alpha=0.15
)+

theme_bw()+

labs(
 title="Comunidad después de la perturbación",
 subtitle="Reorganización adaptativa"
)



###############################################################
# Distribución espacial final
###############################################################

ggplot(
  final,
  aes(
    x,
    y
  )
)+

geom_point(
  alpha=0.15
)+

theme_bw()+

labs(
 title="Distribución espacial posterior",
 subtitle="Expansión hacia nuevos nichos"
)05