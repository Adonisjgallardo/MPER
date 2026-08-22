## ============================================================
## app.R
## Punto de entrada de la app Shiny "MPER: HCA Evolutionary
## Rescue Simulator".
##
## Se mantiene en la ruta de paquete estandar inst/app/ para que
## el proyecto pueda publicarse como paquete de R.
##
## Estrategia de carga (sin tocar globalenv jamas):
##  1) Si el paquete MPER esta instalado, se usa el paquete.
##     Esta es la fuente unica de verdad para apps desplegadas.
##  2) Si no hay instalacion (checkout de fuente), los archivos
##     de R/ se cargan en un entorno PROPIO de esta app. Cargarlos
##     en globalenv provoca mascaras tipo
##     "`mper_server` masks `MPER::mper_server()`" al correr
##     devtools::load_all() despues en la misma sesion.
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
if (requireNamespace("MPER", quietly = TRUE)) {
  MPER::run_mper()
} else {
  ## Localiza la raiz del checkout (directorio con DESCRIPTION y R/)
  ## subiendo desde el directorio de trabajo actual.
  raiz <- local({
    cand <- normalizePath(getwd(), mustWork = FALSE)
    repeat {
      if (file.exists(file.path(cand, "DESCRIPTION")) &&
          dir.exists(file.path(cand, "R"))) return(cand)
      nxt <- dirname(cand)
      if (identical(nxt, cand)) stop("No se encontro la raiz del proyecto MPER.")
      cand <- nxt
    }
  })

  env_app <- new.env(parent = globalenv())
  for (f in sort(list.files(file.path(raiz, "R"), pattern = "[.]R$", full.names = TRUE))) {
    source(f, local = env_app)
  }

  shiny::shinyApp(ui = env_app$mper_ui(), server = env_app$mper_server)
}
