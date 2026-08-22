## ============================================================
## app.R
## Punto de entrada de la app Shiny "HCA Evolutionary Rescue
## Simulator".
##
## Se mantiene en la ruta de paquete estándar inst/app/ para que
## el proyecto pueda publicarse como paquete de R sin depender de
## una carpeta HCARescue/ temporal.
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

## ---- Lanzar la app -------------------------------------------------
if (exists("run_mper", mode = "function")) {
  run_mper()
} else {
  shiny::shinyApp(ui = mper_ui(), server = mper_server)
}
