## Carga el codigo del paquete para los tests.
## SIEMPRE se hace source() de R/ (los tests verifican el codigo fuente
## actual, que puede incluir funciones todavia no presentes en ninguna
## instalacion previa del paquete).

..mper_root <- local({
  candidatos <- character(0)
  args <- commandArgs(FALSE)
  farg <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(farg) && nzchar(farg)) {
    candidatos <- c(candidatos, dirname(normalizePath(farg)))
  }
  candidatos <- c(candidatos, normalizePath(getwd()))

  for (ini in candidatos) {
    cand <- ini
    repeat {
      if (file.exists(file.path(cand, "DESCRIPTION")) &&
          dir.exists(file.path(cand, "R"))) {
        return(cand)
      }
      nxt <- dirname(cand)
      if (identical(nxt, cand)) break
      cand <- nxt
    }
  }
  stop("No se encontro la raiz del paquete MPER (directorio con DESCRIPTION y R/).")
})

for (f in sort(list.files(file.path(..mper_root, "R"), pattern = "[.]R$", full.names = TRUE))) {
  source(f, local = globalenv())
}
