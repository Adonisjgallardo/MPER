#Archivo 01_paisaje.R

#crea una matriz espacial;
#asigna un óptimo fenotípico a cada microambiente;
#genera heterogeneidad;
#permite visualizar gradientes;
#prepara el terreno para la evolución.

#La idea es que cuando después agreguemos bacterias, ellas no estarán en un espacio vacío: estarán explorando un paisaje con estructura ecológica.

###############################################################
# MPER 1.0
# Modelo de Perturbación Ecológica y Reorganización Adaptativa
#
# Archivo 01:
# Creación del paisaje adaptativo
###############################################################

rm(list=ls())

library(ggplot2)

set.seed(123)


###############################################################
# Tamaño del paisaje
###############################################################

nx <- 100
ny <- 100


###############################################################
# Crear coordenadas
###############################################################

paisaje <- expand.grid(
  x = 1:nx,
  y = 1:ny
)


###############################################################
# Crear gradientes ambientales
###############################################################

paisaje$E1 <- 
  scale(paisaje$x)[,1] +
  rnorm(nrow(paisaje),0,0.15)


paisaje$E2 <- 
  scale(paisaje$y)[,1] +
  rnorm(nrow(paisaje),0,0.15)


paisaje$E3 <-
  sin(paisaje$x/10) *
  cos(paisaje$y/10) +
  rnorm(nrow(paisaje),0,0.15)


paisaje$E4 <-
  sqrt(
    (paisaje$x-nx/2)^2+
    (paisaje$y-ny/2)^2
  )

paisaje$E4 <-
  scale(paisaje$E4)[,1]



###############################################################
# Visualización del paisaje
###############################################################

ggplot(
  paisaje,
  aes(x,y,fill=E1))+
    geom_tile()+
    scale_fill_viridis_c()+
    theme_bw()+
    labs(
    title="Paisaje adaptativo inicial",
    subtitle="Microambientes con diferentes óptimos ecológicos",
    fill="E1")



###############################################################
# Visualizar otro eje ambiental
###############################################################

ggplot(
  paisaje,
  aes(
    x,
    y,
    fill=E2
  )
)+
geom_tile()+
scale_fill_viridis_c()+
theme_bw()+
labs(
 title="Segundo componente del paisaje",
 fill="E2"
)



###############################################################
# Guardar paisaje
###############################################################

saveRDS(paisaje,"../data/paisaje_inicial.rds")
