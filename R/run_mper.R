## ============================================================
## run_app.R
## Punto de entrada opcional para lanzar la app desde una sesion
## de R ya abierta en la carpeta del proyecto (alternativa a
## abrir/ejecutar app.R directamente en RStudio).
## ============================================================

#' Launch the Shiny app
#'
#' @name run_mper
#' @return A Shiny app object.
#' @export
run_mper <- function(...) {
  shiny::shinyApp(ui = mper_ui(), server = mper_server, ...)
}
