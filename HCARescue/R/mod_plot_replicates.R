## ============================================================
## mod_plot_replicates.R
## Modulo Shiny: curva de supervivencia (persistencia) a partir
## de multiples replicas, y estadisticas de rescate evolutivo.
## ============================================================

mod_replicates_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::plotOutput(ns("replicate_plot"), height = "440px"),
    shiny::verbatimTextOutput(ns("rescue_stats"))
  )
}

mod_replicates_server <- function(id, reps_reactive) {
  shiny::moduleServer(id, function(input, output, session) {

    tiempos_persistencia <- shiny::reactive({
      reps <- reps_reactive()
      shiny::req(reps)
      tiempos_persistencia_empiricos(reps)
    })

    output$replicate_plot <- shiny::renderPlot({
      tp <- tiempos_persistencia()
      shiny::req(nrow(tp) > 0)
      km <- curva_supervivencia_km(tp)

      ggplot2::ggplot(km, ggplot2::aes(x = tiempo, y = supervivencia)) +
        ggplot2::geom_step(colour = "#1b6ca8", linewidth = 1.1) +
        ggplot2::geom_point(data = tp[!tp$extinguio, , drop = FALSE],
                             ggplot2::aes(x = tiempo_persistencia, y = 0), shape = 3,
                             colour = "#2a9d8f", size = 2,
                             inherit.aes = FALSE) +
        ggplot2::scale_y_continuous(limits = c(0, 1),
                                     labels = function(x) paste0(round(100 * x), "%")) +
        ggplot2::labs(x = "Pasos desde la introduccion del antibiotico",
                      y = "Fraccion de replicas que persisten",
                      title = "Curva de supervivencia (rescate evolutivo)",
                      subtitle = "Estimador producto-limite (Kaplan-Meier); + = replicas censuradas (sobrevivieron todo el horizonte)") +
        tema_hca()
    })

    output$rescue_stats <- shiny::renderText({
      tp <- tiempos_persistencia()
      shiny::req(nrow(tp) > 0)
      n <- nrow(tp)
      n_rescatadas <- sum(!tp$extinguio)
      frac <- n_rescatadas / n
      mediana_persist <- stats::median(tp$tiempo_persistencia[tp$extinguio], na.rm = TRUE)

      sprintf(paste(
        "Replicas totales: %d",
        "Replicas rescatadas (sobrevivieron todo el horizonte): %d (%.1f%%)",
        "Replicas extinguidas: %d",
        "Mediana del tiempo de persistencia (extinguidas): %s pasos post-choque",
        sep = "\n"
      ), n, n_rescatadas, 100 * frac, n - n_rescatadas, fmt_num(mediana_persist, 1))
    })
  })
}
