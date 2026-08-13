## ============================================================
## mod_plot_spatial.R
## Modulo Shiny: vista espacial de la colonia (estado final o
## snapshot elegido del historial), dimension fractal, y
## exportacion opcional de la animacion completa (GIF/MP4).
## ============================================================

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
        shiny::h5("Exportar animacion"),
        shiny::div(style = "color: grey; font-size: 85%;",
          "Requiere el historial espacial completo (checkbox ",
          "\"Guardar historial espacial completo\" en el panel lateral) ",
          "y los paquetes magick/av instalados."),
        shiny::downloadButton(ns("download_gif"), "Descargar GIF", class = "btn-sm")
      ),
      shiny::column(
        width = 8,
        shiny::plotOutput(ns("spatial_plot"), height = "560px")
      )
    )
  )
}

mod_spatial_server <- function(id, sim_reactive) {
  shiny::moduleServer(id, function(input, output, session) {

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
      shiny::sliderInput(ns("paso_snapshot"), "Snapshot (paso):",
                          min = min(pasos_disponibles), max = max(pasos_disponibles),
                          value = max(pasos_disponibles), step = 1,
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

      etiqueta <- if (variable == "resist") "Resistencia (z)" else "Nutriente (N)"

      ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = .data[[variable]])) +
        ggplot2::geom_point(size = 0.7, alpha = 0.75) +
        ggplot2::scale_colour_viridis_c(option = "plasma", name = etiqueta) +
        ggplot2::coord_fixed() +
        ggplot2::labs(title = paste("Paso:", snap$paso), x = "Posicion X", y = "Posicion Y") +
        tema_hca()
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

    output$download_gif <- shiny::downloadHandler(
      filename = function() paste0("hca_evolucion_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".gif"),
      content = function(file) {
        sim <- sim_reactive()
        shiny::validate(shiny::need(!is.null(sim), "Corre una simulacion primero."))
        historial_valido <- Filter(Negate(is.null), sim$historial)
        shiny::validate(shiny::need(length(historial_valido) > 1,
          "Necesitas guardar el historial espacial completo para animar (ver checkbox en el panel lateral)."))
        shiny::validate(shiny::need(.magick_disponible,
          "El paquete 'magick' no esta instalado en este entorno R."))

        variable <- if (is.null(input$variable)) "resist" else input$variable
        animar_historial_hca(historial_valido, archivo_gif = file, archivo_mp4 = NULL,
                              variable = variable)
      }
    )
  })
}
