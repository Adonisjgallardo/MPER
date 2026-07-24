03_dinamica.R

Ahí aparecerá la selección:

cada bacteria lee su microambiente;
se calcula distancia fenotipo-ambiente;
se calcula fitness;
unas aumentan y otras disminuyen;
aparecen patrones espaciales.


La lógica es:

Fenotipo bacteriano⟶Comparacion con ambiente⟶Fitness⟶Cambio poblacional

Todavía no existe antibiótico.

###############################################################
# MPER 1.0
# Modelo de Perturbación Ecológica y Reorganización Adaptativa
#
# Archivo 03:
# Dinámica evolutiva básica
###############################################################

rm(list=ls())

library(ggplot2)
library(dplyr)


set.seed(123)


###############################################################
# Cargar datos
###############################################################

paisaje <- readRDS(
  "paisaje_inicial.rds"
)

bacterias <- readRDS(
  "bacterias_iniciales.rds"
)



###############################################################
# Número de generaciones
###############################################################

generaciones <- 100



###############################################################
# Función de distancia fenotipo-ambiente
###############################################################

calcular_fitness <- function(bacterias,paisaje){

  
  ambiente <- paisaje %>%
    select(
      x,
      y,
      E1,
      E2,
      E3,
      E4
    )

  
  datos <- bacterias %>%
    left_join(
      ambiente,
      by=c("x","y")
    )


  distancia <- sqrt(
    
    (datos$F1-datos$E1)^2 +
    (datos$F2-datos$E2)^2 +
    (datos$F3-datos$E3)^2 +
    (datos$F4-datos$E4)^2
    
  )


  datos$fitness <-
    exp(-distancia)


  return(datos)

}



###############################################################
# Evolución
###############################################################

historial <- list()


for(g in 1:generaciones){

  
  bacterias <-
    calcular_fitness(
      bacterias,
      paisaje
    )

  
  bacterias$generacion <- g

  
  historial[[g]] <- bacterias



  #############################################################
  # Selección
  #############################################################

  sobrevivientes <-
    bacterias %>%
    filter(
      runif(n()) < fitness
    )


  #############################################################
  # Si quedan pocos individuos
  #############################################################

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

  descendencia <-
    sobrevivientes %>%
    slice_sample(
      n=10000,
      replace=TRUE
    )


  descendencia$id <-
    1:nrow(descendencia)



  #############################################################
  # Mutación
  #############################################################

  descendencia$F1 <-
    descendencia$F1 +
    rnorm(
      nrow(descendencia),
      0,
      0.05
    )


  descendencia$F2 <-
    descendencia$F2 +
    rnorm(
      nrow(descendencia),
      0,
      0.05
    )


  descendencia$F3 <-
    descendencia$F3 +
    rnorm(
      nrow(descendencia),
      0,
      0.05
    )


  descendencia$F4 <-
    descendencia$F4 +
    rnorm(
      nrow(descendencia),
      0,
      0.05
    )


  bacterias <-
    descendencia %>%
    select(
      id,
      F1,
      F2,
      F3,
      F4,
      x,
      y
    )

}



###############################################################
# Guardar evolución
###############################################################

saveRDS(
  historial,
  "historial_evolucion.rds"
)



###############################################################
# Visualización final
###############################################################

final <-
  historial[[generaciones]]



ggplot(
  final,
  aes(
    F1,
    F2
  )
)+

geom_point(
  alpha=0.2
)+

theme_bw()+

labs(
 title="Espacio fenotípico después de la selección",
 subtitle="Evolución sin perturbación"
)



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
 title="Distribución espacial final",
 subtitle="Ocupación del paisaje adaptativo"
)

