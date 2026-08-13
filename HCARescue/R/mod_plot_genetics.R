## ============================================================
## mod_plot_genetics.R
## Modulo Shiny: valor genetico medio (g_bar), fenotipo medio
## (z_bar) y su comparacion con la recursion analitica de Lande
## (Caja 3, ecuacion B3.2).
## ============================================================

mod_genetics_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::plotOutput(ns("genetics_plot"), height = "420px"),
    shiny::hr(),
    shiny::plotOutput(ns("gbar_overlay"), height = "320px"),
    shiny::div(style = "color: grey; font-size: 90%; margin-top: 4px;",
      "g_bar_analitico usa la recursion de Lande con varianza genetica ",
      "aditiva constante (k=0: optimo estatico tras el choque abrupto)."
    )
  )
}

mod_genetics_server <- function(id, sim_reactive, sigma_g2, sigma_e2) {
  shiny::moduleServer(id, function(input, output, session) {

    output$genetics_plot <- shiny::renderPlot({
      shiny::req(sim_reactive())
      sim <- sim_reactive()
      df <- sim$resumen
      intro <- sim$paso_introduccion
      theta0 <- sim$parametros$theta0
      theta1 <- sim$parametros$theta1

      df_long <- data.frame(
        paso = rep(df$paso, 2),
        valor = c(df$g_bar, df$z_bar),
        serie = rep(c("Valor genetico medio (g_bar)", "Fenotipo medio (z_bar)"), each = nrow(df))
      )

      ggplot2::ggplot(df_long, ggplot2::aes(x = paso, y = valor, colour = serie)) +
        ggplot2::geom_line(linewidth = 1.1, na.rm = TRUE) +
        ggplot2::geom_hline(yintercept = theta0, linetype = "dotted", colour = "grey40") +
        ggplot2::geom_hline(yintercept = theta1, linetype = "dotted", colour = "grey40") +
        ggplot2::annotate("text", x = min(df$paso), y = theta0, label = "theta0", vjust = -0.5,
                           hjust = 0, colour = "grey40", size = 3) +
        ggplot2::annotate("text", x = min(df$paso), y = theta1, label = "theta1", vjust = -0.5,
                           hjust = 0, colour = "grey40", size = 3) +
        ggplot2::geom_vline(xintercept = intro, linetype = "dashed", colour = "firebrick") +
        ggplot2::scale_colour_manual(values = c("Valor genetico medio (g_bar)" = "#7b3f9e",
                                                 "Fenotipo medio (z_bar)" = "#2a9d8f")) +
        ggplot2::labs(x = "Paso de tiempo", y = "Valor de rasgo",
                      title = "Respuesta evolutiva del rasgo de resistencia",
                      colour = NULL) +
        tema_hca()
    })

    output$gbar_overlay <- shiny::renderPlot({
      shiny::req(sim_reactive())
      sim <- sim_reactive()

      cmp <- tryCatch(
        comparar_gbar_caja3(sim, k = 0, sigma_g2 = sigma_g2(), sigma_e2 = sigma_e2()),
        error = function(e) NULL
      )
      shiny::validate(shiny::need(!is.null(cmp), "No se pudo calcular g_bar analitico todavia."))

      cmp_long <- data.frame(
        t = rep(cmp$t, 2),
        g = c(cmp$g_bar_sim, cmp$g_bar_analitico),
        serie = rep(c("Simulacion (HCA)", "Analitico (Lande, Caja 3)"), each = nrow(cmp))
      )

      ggplot2::ggplot(cmp_long, ggplot2::aes(x = t, y = g, colour = serie, linetype = serie)) +
        ggplot2::geom_line(linewidth = 1.1, na.rm = TRUE) +
        ggplot2::scale_colour_manual(values = c("Simulacion (HCA)" = "#7b3f9e",
                                                 "Analitico (Lande, Caja 3)" = "#e07b39")) +
        ggplot2::labs(x = "Pasos desde la introduccion del antibiotico",
                      y = "g_bar", title = "g_bar simulado vs. recursion de Lande",
                      colour = NULL, linetype = NULL) +
        tema_hca()
    })
  })
}
