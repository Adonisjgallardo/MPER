#Archivo 04_disturbio.R

#El antibiótico no crea resistencia. Modifica el paisaje adaptativo.

#No modificamos ninguna bacteria.

#No agregamos genes.

#No aumentamos F1.

#Solo cambiamos el ambiente.

###############################################################
# MPER 1.0
# Modelo de Perturbación Ecológica y Reorganización Adaptativa
#
# Archivo 04:
# Introducción del disturbio ecológico
###############################################################

rm(list=ls())

library(ggplot2)
library(dplyr)


set.seed(123)


###############################################################
# Cargar estado previo
###############################################################

paisaje <- readRDS(
  "paisaje_inicial.rds"
)

historial <- readRDS(
  "historial_evolucion.rds"
)


bacterias <-
  historial[[100]]



###############################################################
# Visualizar paisaje antes del disturbio
###############################################################

ggplot(
  paisaje,
  aes(
    x,
    y,
    fill=E1
  )
)+
geom_tile()+
scale_fill_viridis_c()+
theme_bw()+
labs(
 title="Paisaje antes de la perturbación",
 fill="E1"
)



###############################################################
# Definir zona de perturbación
###############################################################

zona <- paisaje$x > 65 &
        paisaje$x < 90 &
        paisaje$y > 30 &
        paisaje$y < 70



###############################################################
# Aplicar disturbio
#
# El ambiente cambia
# Las bacterias NO cambian
###############################################################


paisaje_disturbio <-
  paisaje


paisaje_disturbio$E1[zona] <-
  paisaje_disturbio$E1[zona] + 5



###############################################################
# Visualizar nuevo paisaje
###############################################################

ggplot(
  paisaje_disturbio,
  aes(
    x,
    y,
    fill=E1
  )
)+
geom_tile()+
scale_fill_viridis_c()+
theme_bw()+
labs(
 title="Paisaje después del disturbio",
 subtitle="Cambio del ambiente adaptativo",
 fill="E1"
)



###############################################################
# Comparar población antes del disturbio
###############################################################

ggplot(
  bacterias,
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
 title="Comunidad antes del disturbio",
 subtitle="La diversidad ya existe"
)

# Difusion fenotipica
# Proceso de orstein - uhleeck

dE = v(E-F)dt + s^2 dW 

# Difusion espacial
# reaccion-difusion Focker-Planck

# Metricas Diversidad
# Shannon
# Distancia 

# Visualizaciones
# Entorno 100x100 celdas

# Si los camios son rapidos no se diferencian en la medicion
# Transferencia Lateral y presion selectiva

# TESIS
# 

###############################################################
# Guardar nuevo ambiente
###############################################################

saveRDS(
  paisaje_disturbio,
  "paisaje_perturbado.rds"
)


saveRDS(
  bacterias,
  "bacterias_pre_disturbio.rds"
)
