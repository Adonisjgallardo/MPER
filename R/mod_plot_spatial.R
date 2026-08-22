## ============================================================
## mod_plot_spatial.R
## Modulo Shiny: vista espacial de la colonia (estado final o
## snapshot elegido del historial), dimension fractal, y
## animacion completa de la simulacion embebida como video MP4
## (con descarga opcional). El MP4 se genera con av::av_encode_video()
## (FFmpeg estatico incluido en el paquete `av`), que es mas
## liviano y mas confiable de visualizar en navegador que un GIF
## ensamblado con magick.
## ============================================================

#' Grafica un snapshot con ejes FIJOS a la reticula completa
#'
#' A diferencia de dejar que ggplot elija los limites a partir de los
#' sitios ocupados, aqui siempre se usa la reticula entera
#' [0.5, nx+0.5] x [0.5, ny+0.5], de modo que la posicion y el tamano
#' del grafico NO cambian entre snapshots ni durante la reproduccion.
#'
#' @param snap lista con `estado`, `N`, `resist`, `paso` (un snapshot
#'        tal como los guarda `simular_hca()`)
#' @param variable "resist" o "N": que variable colorear
#' @param point_size,alpha estetica de los puntos
#'
#' @return un objeto ggplot
plot_snapshot_estatico <- function(snap, variable = c("resist", "N"),
                                    point_size = 0.7, alpha = 0.75) {
  variable <- match.arg(variable)
  dims <- dim(snap$estado)
  nx <- if (is.null(dims)) 100L else as.integer(dims[1])
  ny <- if (is.null(dims)) 100L else as.integer(dims[2])

  df <- snapshot_a_df(snap$estado, snap$N, snap$resist, paso = snap$paso)
  etiqueta <- if (variable == "resist") "Resistencia (z)" else "Nutriente (N)"
  limites_color <- if (nrow(df) > 0) range(df[[variable]], na.rm = TRUE) else c(NA_real_, NA_real_)

  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = .data[[variable]])) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::scale_colour_viridis_c(option = "plasma", name = etiqueta,
                                     limits = limites_color) +
    ggplot2::coord_fixed(xlim = c(0.5, nx + 0.5), ylim = c(0.5, ny + 0.5),
                          expand = FALSE) +
    ggplot2::labs(title = paste("Paso:", snap$paso),
                  x = "Posicion X", y = "Posicion Y") +
    tema_hca()
}

mod_spatial_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(
        width = 4,
        shiny::selectInput(ns("variable"), "Colorear por:",
                            choices = c("Resistencia (z)" = "resist", "Nutriente (N)" = "N"),
                            selected = "resist"),
        shiny::uiOutput(ns("snapshot_slider_ui")),
        shiny::verbatimTextOutput(ns("fractal_info")),
        shiny::hr(),
        shiny::h5("Animacion (video)"),
        shiny::div(style = "color: grey; font-size: 85%;",
          "Genera un video MP4 de la evolucion completa. Requiere el historial ",
          "espacial completo (checkbox \"Guardar historial espacial completo\" en ",
          "el panel lateral), el paquete av instalado y un navegador con soporte HTML5."),
        shiny::actionButton(ns("generar_anim"), "Generar animacion", class = "btn-sm"),
        shiny::downloadButton(ns("download_mp4"), "Descargar MP4", class = "btn-sm")
      ),
      shiny::column(
        width = 8,
        shiny::plotOutput(ns("spatial_plot"), height = "560px"),
        shiny::uiOutput(ns("animacion_ui"))
      )
    )
  )
}

mod_spatial_server <- function(id, sim_reactive) {
  shiny::moduleServer(id, function(input, output, session) {

    ## Directorio y prefijo de recurso por-sesion para servir el MP4.
    ## addResourcePath publica `anim_dir` bajo /<prefijo>/... en la URL
    ## de la app, de modo que <video src> pueda apuntar al archivo.
    anim_dir    <- file.path(tempdir(), sprintf("mper_anim_%s", session$token))
    anim_prefix <- paste0("mper_anim_", session$token)
    anim_archivo <- "evolucion.mp4"
    dir.create(anim_dir, showWarnings = FALSE, recursive = TRUE)
    shiny::addResourcePath(anim_prefix, anim_dir)

    rv_anim <- shiny::reactiveValues(mp4_generado = FALSE, generando = FALSE)

    ## Al correr una nueva simulacion, invalida el video anterior para
    ## no mostrar una animacion desactualizada. Se usa isolate() para NO
    ## crear dependencias sobre rv_anim: si el observador dependiera de
    ## `generando`, se re-ejecutaria justo tras terminar la generacion y
    ## borraria el mp4_generado=TRUE recien asignado.
    shiny::observe({
      sim_reactive()
      shiny::isolate({
        rv_anim$mp4_generado <- FALSE
        if (!rv_anim$generando) shinyjs::disable("download_mp4")
      })
    })

    ## El boton de descarga arranca deshabilitado: solo tiene sentido
    ## tras generar la animacion (asi la descarga es una copia local,
    ## sin regenerar nada ni validar dentro del content del handler).
    shinyjs::disable("download_mp4")

    output$snapshot_slider_ui <- shiny::renderUI({
      ns <- session$ns
      sim <- sim_reactive()
      shiny::req(sim)
      pasos_disponibles <- vapply(sim$historial, function(s) if (is.null(s)) NA_integer_ else s$paso, integer(1))
      pasos_disponibles <- pasos_disponibles[!is.na(pasos_disponibles)]
      if (length(pasos_disponibles) <= 1) {
        return(shiny::div(style = "color: grey; font-size: 85%;",
          "Solo se guardo el estado final (activa \"Guardar historial completo\" para explorar snapshots intermedios)."))
      }
      ## Rango completo: desde el paso 1 hasta la extincion (primer cruce
      ## bajo Nc post-shock) o, si la poblacion persiste, hasta el ultimo
      ## paso real de la corrida. Los pasos sin snapshot guardado muestran
      ## el snapshot mas cercano.
      ultimo_paso <- if (!is.na(sim$paso_extincion)) sim$paso_extincion else max(sim$resumen$paso)
      shiny::sliderInput(ns("paso_snapshot"), "Snapshot (paso):",
                          min = 1L, max = ultimo_paso,
                          value = ultimo_paso, step = 1L,
                          animate = shiny::animationOptions(interval = 400))
    })

    snapshot_actual <- shiny::reactive({
      sim <- sim_reactive()
      shiny::req(sim)
      historial_valido <- Filter(Negate(is.null), sim$historial)
      if (length(historial_valido) == 0) {
        return(list(paso = utils::tail(sim$resumen$paso, 1),
                     estado = sim$estado_final, N = sim$N_final,
                     resist = sim$resist_final))
      }
      if (is.null(input$paso_snapshot)) {
        return(historial_valido[[length(historial_valido)]])
      }
      pasos <- vapply(historial_valido, function(s) s$paso, integer(1))
      historial_valido[[which.min(abs(pasos - input$paso_snapshot))]]
    })

    output$spatial_plot <- shiny::renderPlot({
      snap <- snapshot_actual()
      shiny::req(snap)
      variable <- if (is.null(input$variable)) "resist" else input$variable

      df <- snapshot_a_df(snap$estado, snap$N, snap$resist, paso = snap$paso)
      shiny::validate(shiny::need(nrow(df) > 0, "La colonia esta extinta en este paso: no hay sitios ocupados."))

      plot_snapshot_estatico(snap, variable)
    })

    output$fractal_info <- shiny::renderText({
      snap <- snapshot_actual()
      shiny::req(snap)
      res <- calcular_dimension_fractal(snap$estado)
      if (is.na(res$dimension)) {
        "Dimension fractal: N/A (colonia extinta o insuficientes sitios ocupados)"
      } else {
        sprintf("Dimension fractal (box-counting): %.3f\nSitios ocupados: %d",
                res$dimension, sum(snap$estado, na.rm = TRUE))
      }
    })

    historial_para_animar <- shiny::eventReactive(input$generar_anim, {
      sim <- sim_reactive()
      shiny::validate(shiny::need(!is.null(sim), "Corre una simulacion primero."))
      historial_valido <- Filter(Negate(is.null), sim$historial)
      shiny::validate(shiny::need(length(historial_valido) > 1,
        "Necesitas guardar el historial espacial completo para animar (ver checkbox en el panel lateral)."))
      historial_valido
    }, ignoreNULL = FALSE)

    shiny::observeEvent(input$generar_anim, {
      if (rv_anim$generando) return(invisible(NULL))

      if (!.av_disponible) {
        shiny::showNotification(
          "El paquete 'av' no esta instalado en este entorno R (install.packages('av')).",
          type = "error", duration = NULL)
        return(invisible(NULL))
      }

      hist <- tryCatch(historial_para_animar(), error = function(e) e)
      if (inherits(hist, "error")) {
        shiny::showNotification(conditionMessage(hist), type = "warning", duration = 8)
        return(invisible(NULL))
      }

      variable <- if (is.null(input$variable)) "resist" else input$variable
      destino  <- file.path(anim_dir, anim_archivo)

      rv_anim$generando <- TRUE
      shinyjs::disable("generar_anim"); shinyjs::disable("download_mp4")

      resultado <- tryCatch({
        shiny::withProgress(message = "Generando animacion (MP4)", value = 0.3, {
          unlink(destino)
          ## Ejes fijos a la reticula completa para que el video no haga
          ## zoom conforme crece/decrece la colonia.
          dims_anim <- dim(hist[[1]]$estado)
          animar_historial_hca(hist, archivo_gif = NULL, archivo_mp4 = destino,
                                variable = variable,
                                nx = if (is.null(dims_anim)) NULL else as.integer(dims_anim[1]),
                                ny = if (is.null(dims_anim)) NULL else as.integer(dims_anim[2]))
          shiny::setProgress(value = 1)
        })
        destino
      }, error = function(e) e)

      rv_anim$generando <- FALSE
      shinyjs::enable("generar_anim")

      if (inherits(resultado, "error")) {
        shiny::showNotification(paste("Error generando la animacion:",
                                       conditionMessage(resultado)),
                                 type = "error", duration = NULL)
        return(invisible(NULL))
      }
      if (!file.exists(resultado)) {
        shiny::showNotification("No se pudo escribir el archivo de video.",
                                 type = "error", duration = NULL)
        return(invisible(NULL))
      }

      rv_anim$mp4_generado <- TRUE
      shinyjs::enable("download_mp4")
      shiny::showNotification("Animacion generada.", type = "message", duration = 5)
      invisible(NULL)
    })

    output$animacion_ui <- shiny::renderUI({
      shiny::req(rv_anim$mp4_generado)
      shiny::tags$div(
        style = "margin-top: 12px;",
        shiny::tags$label("Evolucion completa (puede tardar unos segundos en cargar):"),
        shiny::tags$video(
          src = file.path(anim_prefix, anim_archivo),
          type = "video/mp4",
          controls = NA, loop = NA, muted = NA, autoplay = NA,
          preload = "auto",
          style = "width: 100%; max-width: 640px; max-height: 560px; background: #222;"
        )
      )
    })

    output$download_mp4 <- shiny::downloadHandler(
      filename = function() paste0("hca_evolucion_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".mp4"),
      content = function(file) {
        origen <- file.path(anim_dir, anim_archivo)
        if (rv_anim$mp4_generado && file.exists(origen)) {
          file.copy(origen, file, overwrite = TRUE)
        } else {
          sim <- sim_reactive()
          shiny::req(sim)
          historial_valido <- Filter(Negate(is.null), sim$historial)
          if (length(historial_valido) <= 1 || !.av_disponible) {
            stop("Genera primero la animacion (boton \"Generar animacion\").")
          }
          variable <- if (is.null(input$variable)) "resist" else input$variable
          animar_historial_hca(historial_valido, archivo_gif = NULL,
                                archivo_mp4 = file, variable = variable)
        }
      }
    )
  })
}
