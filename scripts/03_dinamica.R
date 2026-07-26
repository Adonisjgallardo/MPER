#03_dinamica.R

#Ahí aparecerá la selección:
# cada bacteria lee su microambiente;
#se calcula distancia fenotipo-ambiente;
#se calcula fitness;
#unas aumentan y otras disminuyen;
#aparecen patrones espaciales.

#La lógica es:

#Fenotipo bacteriano⟶Comparacion con ambiente⟶Fitness⟶Cambio poblacional

#Todavía no existe antibiótico.

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
library(tidyr)

set.seed(123)

###############################################################
# Cargar datos
###############################################################
bacterias <- readRDS("data/bacterias.rds") |> select(-c(E2, E3, E4))

###############################################################
# Función de distancia fenotipo-ambiente
# Primera generacion de bacterias
###############################################################

bacterias$dist <- bacterias |> 
  apply(1, function(row) {
  dist(rbind(row[c("E1","E1","E1","E1")], row[c("F1","F2","F3","F4")]))
})

bacterias$fitness <- exp(-bacterias$dist)

###############################################################
# Evolución
###############################################################
## Número de generaciones
###############################################################

generaciones <- nrow(bacterias)

bacterias$viva <- (runif(generaciones) < bacterias$fitness)

historial <- list()
historial[[1]] <- bacterias

for (g in 1:generaciones) {
  # Current survivors
  sobrevivientes <- historial[[g]] |> filter(viva)
  
  # If too few, replenish with 50 random individuals from the initial pool
  if (nrow(sobrevivientes) < 50) {
    sobrevivientes <- bacterias |> slice_sample(n = 50)
  }
  
  # ---- Diffusion / Reproduction ----
  # 1. Random threshold for phenotypic selection
  threshold <- runif(1, 0, 0.5)
  # 2. Identify individuals that pass the condition
  padres <- sobrevivientes |> filter(F1*F2*F3 < threshold)
  n_offspring <- nrow(padres)
  
  # If no one passes, we still need to create a next generation.
  # Option: create offspring from all survivors (or keep the same population).
  # Here we'll sample from all survivors if n_offspring == 0, but you could adjust.
  if (n_offspring == 0) {
    # Fallback: produce offspring from all survivors (with replacement)
    # to maintain population size (optional)
    padres <- sobrevivientes
    n_offspring <- nrow(padres)
  }
  
  # 3. Generate offspring by sampling from the selected parents (with replacement)
  descendencia <- padres |>
    slice_sample(n = n_offspring, replace = TRUE) |>
    mutate(
      # Spatial random walk (x, y)
      across(x:y, ~ .x + round(runif(n_offspring, -1, 1))),
      # Keep within 1..100 (circular or boundary reflection)
      x = ifelse(x > 100, x - 100, ifelse(x < 1, x + 100, x)),
      y = ifelse(y > 100, y - 100, ifelse(y < 1, y + 100, y)),
      # Increment generation
      generacion = generacion + 1,
      # Phenotypic diffusion (add independent Gaussian noise to each trait)
      F2 = F2 + rnorm(n_offspring, 0, 0.05),
      F3 = F3 + rnorm(n_offspring, 0, 0.05),
      F4 = F4 + rnorm(n_offspring, 0, 0.05),
      # Set viva to TRUE for offspring (or you can later assign based on fitness)
      viva = TRUE
    )
  
  # Store the new generation
  descendencia$dist <- descendencia |> 
  apply(1, function(row) {
    dist(rbind(row[c("E1","E1","E1","E1")], row[c("F1","F2","F3","F4")]))})

descendencia$fitness <- exp(-descendencia$dist)
descendencia$viva <- (runif(n_offspring) < descendencia$fitness)
historial[[g + 1]] <- descendencia
}
rm(list=ls()[-which(ls() %in% "historial")])
###############################################################
# Guardar evolución
###############################################################

saveRDS(historial,"data/historial_evolucion.rds")

###############################################################
# Visualización final
###############################################################

final <- historial[[10000]]

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

# Animation
library(ggplot2)
library(gganimate)
library(dplyr)
library(tidyr)   # for bind_rows

# Combine all data frames from the list into one, adding a generation index
# (if the 'generacion' column is already correct, you can skip the .id step)
df_all <- bind_rows(historial, .id = "gen") %>%
  transmute(
    x = x, 
    y = y,
    gen = as.integer(gen) - 1,          # convert to 0‑based generations
    sum_F = F1*F2*F3              # compute the colour variable
  )

# Optional: for faster rendering, you can sample a fraction of points per generation
# df_all <- df_all %>% group_by(gen) %>% slice_sample(prop = 0.2) %>% ungroup()

# Create the animated plot
p <- ggplot(df_all, aes(x = x, y = y, colour = sum_F)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_colour_viridis_c(option = "plasma", name = "F1(F2)(F3)") +
  coord_fixed() +                       # keep aspect ratio square
  labs(
    title = "Generation: {frame_time}",
    x = "X position",
    y = "Y position"
  ) +
  theme_minimal() +
  transition_time(gen) +
  ease_aes('linear')

# Render the animation
animate(
  p,
  nframes = length(historial),          # one frame per generation
  fps = 10,
  width = 600,
  height = 600,
  renderer = gifski_renderer("plots/evolution.gif")   # saves as GIF
)

animate(
  p,
  nframes = length(historial),          # one frame per generation
  fps = 10,
  width = 600,
  height = 600,
  renderer = av_renderer("plots/evolution.mp4")   # saves as video
)

# To display in RStudio viewer/plot pane, just run `animate(p)`