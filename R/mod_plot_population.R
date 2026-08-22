## ============================================================
## mod_plot_population.R
## Modulo Shiny: dinamica poblacional (N_total) + comparacion
## contra la prediccion analitica de la Caja 3.
## ============================================================

mod_pop_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::plotOutput(ns("pop_plot"), height = "380px"),
    shiny::hr(),
    shiny::h5("Comparacion con la prediccion analitica (Caja 3)"),
    shiny::plotOutput(ns("analytical_overlay"), height = "320px"),
    shiny::div(style = "color: grey; font-size: 90%; margin-top: 4px;",
      "N_t bajo un cambio ambiental abrupto, evolucion deterministica ",
      "densidad-independiente (formula B3.1 del documento de rescate evolutivo)."
    )
  )
}

mod_pop_server <- function(id, sim_reactive, sigma_g2, sigma_e2, Wmax) {
  shiny::moduleServer(id, function(input, output, session) {

    output$pop_plot <- shiny::renderPlot({
      shiny::req(sim_reactive())
      sim <- sim_reactive()
      df <- sim$resumen
      intro <- sim$paso_introduccion

      ## Marcadores de eventos antibioticos (introduccion + multirresistencia)
      eventos <- if (!is.null(sim$eventos_shock) && nrow(sim$eventos_shock) > 0) {
        subset(sim$eventos_shock, tipo != "introduccion")
      } else {
        NULL
      }

      p <- ggplot2::ggplot(df, ggplot2::aes(x = paso)) +
        ggplot2::geom_line(ggplot2::aes(y = n_ocupados, colour = "Poblacion (N)"), linewidth = 1.1) +
        ggplot2::geom_line(ggplot2::aes(y = n_frente, colour = "Celulas proliferativas"),
                            linewidth = 0.7, linetype = "dashed") +
        ggplot2::geom_vline(xintercept = intro, linetype = "dashed", colour = "firebrick") +
        ggplot2::annotate("text", x = intro, y = max(df$n_ocupados, na.rm = TRUE),
                           label = "  introduccion antibiotico", hjust = 0, vjust = 1,
                           colour = "firebrick", size = 3.2)

      if (!is.null(eventos)) {
        y_tope <- max(df$n_ocupados, na.rm = TRUE)
        p <- p +
          ggplot2::geom_vline(data = eventos,
                              ggplot2::aes(xintercept = paso, linetype = tipo),
                              colour = "darkorange3", linewidth = 0.8) +
          ggplot2::geom_text(data = eventos,
                             ggplot2::aes(x = paso, y = y_tope,
                                          label = paste0(" nivel ", nivel + 1L, " (", tipo, ")")),
                             hjust = -0.05, vjust = 1.4,
                             colour = "darkorange3", size = 3.1, inherit.aes = FALSE)
      }

      p <- p +
        ggplot2::scale_colour_manual(values = c("Poblacion (N)" = "#1b6ca8",
                                                 "Celulas proliferativas" = "#4caf50")) +
        ggplot2::labs(y = "Numero de sitios ocupados", x = "Paso de tiempo",
                      title = "Dinamica poblacional",
                      subtitle = if (!is.na(sim$paso_extincion))
                        paste0("Cruce del umbral critico Nc en el paso ", sim$paso_extincion)
                      else "Poblacion se mantuvo sobre el umbral critico Nc durante toda la simulacion",
                      colour = NULL) +
        tema_hca()
      p
    })

    output$analytical_overlay <- shiny::renderPlot({
      shiny::req(sim_reactive())
      sim <- sim_reactive()

      cmp <- tryCatch(
        comparar_Nt_caja3(sim, Wmax = Wmax(), sigma_g2 = sigma_g2(), sigma_e2 = sigma_e2()),
        error = function(e) NULL
      )
      shiny::validate(shiny::need(!is.null(cmp), "No se pudo calcular la comparacion analitica todavia."))

      cmp_long <- data.frame(
        t = rep(cmp$t, 2),
        N = c(cmp$N_sim, cmp$N_analitico),
        serie = rep(c("Simulacion (HCA)", "Analitico (Caja 3)"), each = nrow(cmp))
      )

      ggplot2::ggplot(cmp_long, ggplot2::aes(x = t, y = N, colour = serie, linetype = serie)) +
        ggplot2::geom_line(linewidth = 1.1) +
        ggplot2::scale_colour_manual(values = c("Simulacion (HCA)" = "#1b6ca8",
                                                 "Analitico (Caja 3)" = "#e07b39")) +
        ggplot2::labs(x = "Pasos desde la introduccion del antibiotico",
                      y = "Tamaño poblacional",
                      title = "N_t simulado vs. prediccion analitica",
                      colour = NULL, linetype = NULL) +
        tema_hca()
    })
  })
}
