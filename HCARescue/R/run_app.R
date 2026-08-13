## ============================================================
## run_app.R
## Punto de entrada opcional para lanzar la app desde una sesion
## de R ya abierta en la carpeta del proyecto (alternativa a
## abrir/ejecutar app.R directamente en RStudio).
## ============================================================

run_hca_app <- function(...) {
  shiny::shinyApp(ui = hca_ui(), server = hca_server, ...)
}
