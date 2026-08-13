## ============================================================
## app.R
## Punto de entrada de la app Shiny "HCA Evolutionary Rescue
## Simulator". Ejecutar con:
##
##   shiny::runApp("ruta/a/HCARescue")
##
## o abrir este archivo en RStudio y pulsar "Run App".
##
## Estructura (ver README.md para más detalle):
##   R/01_campos.R ... R/07_analisis.R   -> motor de simulacion
##   R/utils.R                           -> helpers compartidos
##   R/mod_plot_*.R                      -> modulos Shiny (graficos)
##   R/app_ui.R, R/app_server.R          -> UI y servidor
##   R/run_app.R                         -> lanzador opcional
##   www/styles.css                      -> estilos
## ============================================================

## ---- Paquetes obligatorios ----------------------------------------
paquetes_requeridos <- c("shiny", "ggplot2", "bslib", "shinyjs")
faltantes <- paquetes_requeridos[!vapply(paquetes_requeridos, requireNamespace,
                                          logical(1), quietly = TRUE)]
if (length(faltantes) > 0) {
  stop(
    "Faltan paquetes obligatorios: ", paste(faltantes, collapse = ", "),
    "\nInstalalos con: install.packages(c(", paste0('"', faltantes, '"', collapse = ", "), "))"
  )
}

library(shiny)
library(ggplot2)

## Paquetes opcionales (mejoran la experiencia pero la app funciona sin ellos):
##  - future / promises : simulaciones no bloqueantes (ver app_server.R)
##  - magick / av        : exportar animaciones GIF/MP4 (ver 07_analisis.R)
##  - parallel / foreach / doParallel : paralelizar simular_hca_replicas()

## ---- Cargar el motor de simulacion y la app ----------------------
## Shiny ejecuta app.R con el directorio de trabajo ya puesto en la
## carpeta de la app (HCARescue/), asi que las rutas relativas "R/..."
## funcionan tanto con shiny::runApp() como con el boton "Run App".
archivos_r <- c(
  "01_campos.R", "02_colonia.R", "03_pde_difusion.R", "04_ca_crecimiento.R",
  "05_motor.R", "06_genetica_cuantitativa.R", "07_analisis.R",
  "utils.R",
  "mod_plot_population.R", "mod_plot_genetics.R",
  "mod_plot_spatial.R", "mod_plot_replicates.R",
  "app_ui.R", "app_server.R"
)

for (f in archivos_r) {
  ruta <- file.path("R", f)
  if (!file.exists(ruta)) stop("No se encontro el archivo requerido: ", ruta,
                                "\n(¿Estas corriendo la app desde la carpeta HCARescue/?)")
  source(ruta, local = FALSE)
}

## ---- Lanzar la app -------------------------------------------------
shinyApp(ui = hca_ui(), server = hca_server)
