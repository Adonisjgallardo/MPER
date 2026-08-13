## ============================================================
## app_ui.R
## Definicion de la interfaz de usuario de HCARescue.
## ============================================================

hca_ui <- function() {
  shiny::fluidPage(
    theme = bslib::bs_theme(bootswatch = "flatly", primary = "#1b6ca8"),
    shinyjs::useShinyjs(),
    shiny::tags$head(
      shiny::tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),

    shiny::titlePanel("HCA Evolutionary Rescue Simulator",
                       windowTitle = "HCA Rescate Evolutivo"),
    shiny::div(class = "hca-subtitle",
      "Autómata celular híbrido (PDE nutriente/antibiótico + genética cuantitativa) ",
      "para rescate evolutivo bajo estrés antibiótico espacialmente heterogéneo."),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,

        shiny::h5("Retícula y tiempo"),
        shiny::sliderInput("nx", "Tamaño de la retícula (nx = ny)", 30, 200, 100, step = 10),
        shiny::sliderInput("n_steps", "Pasos de simulación", 100, 1200, 600, step = 50),
        shiny::checkboxInput("mantener_historial", "Guardar historial espacial completo (para animar)", value = FALSE),
        shiny::conditionalPanel(
          condition = "input.mantener_historial == true",
          shiny::sliderInput("guardar_cada", "Guardar snapshot cada N pasos", 5, 100, 20, step = 5)
        ),

        shiny::h5("Nutriente"),
        shiny::sliderInput("D", "D_N (difusión)", 0, 1, 0.20, step = 0.01),
        shiny::sliderInput("k", "Tasa de consumo (k)", 0, 1, 0.15, step = 0.01),
        shiny::sliderInput("N0", "Nutriente inicial (N0)", 0.1, 2, 1.0, step = 0.1),
        shiny::sliderInput("N_umbral", "Umbral proliferativo (N_umbral)", 0, 1, 0.15, step = 0.01),
        shiny::sliderInput("N_half", "Semisaturación (N_half)", 0.05, 1, 0.30, step = 0.05),
        shiny::sliderInput("p_max", "Probabilidad máxima de división (p_max)", 0, 1, 0.6, step = 0.05),

        shiny::h5("Choque de antibiótico"),
        shiny::sliderInput("paso_introduccion", "Paso de introducción", 0, 500, 150, step = 10),
        shiny::sliderInput("A_max", "Concentración máxima (A_max)", 0, 5, 1.0, step = 0.1),
        shiny::selectInput("tipo_plantilla", "Gradiente espacial",
                            choices = c("Lineal" = "lineal", "Uniforme" = "uniforme", "Radial" = "radial"),
                            selected = "lineal"),
        shiny::selectInput("modo_temporal", "Modulación temporal (post-choque)",
                            choices = c("Constante" = "constante", "Sinusoidal" = "sinusoidal", "Pulsos" = "pulsos"),
                            selected = "constante"),
        shiny::conditionalPanel(
          condition = "input.modo_temporal != 'constante'",
          shiny::sliderInput("periodo_temporal", "Periodo de fluctuación (pasos)", 5, 200, 40, step = 5)
        ),
        shiny::sliderInput("D_A", "D_A (difusión antibiótico)", 0, 1, 0.15, step = 0.01),
        shiny::sliderInput("delta_A", "Degradación natural (delta_A)", 0, 0.5, 0.02, step = 0.01),
        shiny::sliderInput("tasa_dosificacion", "Tasa de dosificación", 0, 1, 0.30, step = 0.05),

        shiny::h5("Genética cuantitativa (Caja 3)"),
        shiny::sliderInput("g_inicial", "Valor genético inicial (g_inicial)", 0, 1, 0.10, step = 0.01),
        shiny::sliderInput("sigma_e", "SD ambiental (sigma_e)", 0, 0.3, 0.05, step = 0.01),
        shiny::sliderInput("mutacion_sd", "SD de mutación", 0, 0.2, 0.02, step = 0.005),
        shiny::sliderInput("theta0", "Óptimo sin antibiótico (theta0)", 0, 1, 0.10, step = 0.01),
        shiny::sliderInput("theta1", "Óptimo con antibiótico saturante (theta1)", 0, 1, 0.80, step = 0.01),
        shiny::sliderInput("omega2", "Ancho² de selección (omega2)", 0.005, 0.5, 0.05, step = 0.005),
        shiny::sliderInput("K_theta", "Semisaturación óptimo (K_theta)", 0.05, 1, 0.30, step = 0.05),

        shiny::h5("Mortalidad"),
        shiny::sliderInput("mort_base", "Mortalidad basal", 0, 0.2, 0.01, step = 0.005),
        shiny::sliderInput("mort_estres", "Mortalidad máx. por estrés", 0, 1, 0.90, step = 0.05),
        shiny::sliderInput("K_mort", "Semisaturación mortalidad (K_mort)", 0.05, 1, 0.30, step = 0.05),
        shiny::sliderInput("Nc", "Umbral crítico de rescate (Nc)", 0, 1000, 200, step = 10),

        shiny::h5("Comparación analítica (Caja 3)"),
        shiny::sliderInput("Wmax", "Aptitud máxima (Wmax)", 0.5, 1.5, 1.0, step = 0.05),
        shiny::sliderInput("sigma_g2", "Varianza genética aditiva asumida (sigma_g2)", 0.0001, 0.05, 0.002, step = 0.0001),

        shiny::hr(),
        shiny::actionButton("run", "\u25b6 Ejecutar simulación", icon = shiny::icon("play"),
                             class = "btn-primary btn-lg", width = "100%"),
        shiny::br(), shiny::br(),
        shiny::numericInput("n_reps", "Número de réplicas", value = 20, min = 2, max = 200, step = 1),
        shiny::actionButton("run_reps", "\u23f3 Ejecutar réplicas", icon = shiny::icon("clone"),
                             class = "btn-success", width = "100%"),
        shiny::br(), shiny::br(),
        shiny::downloadButton("download_results", "\U0001F4BE Descargar resultados (.rds)", class = "btn-info"),
        shiny::br(), shiny::br(),
        shiny::div(id = "status_msg", class = "hca-status", "Listo.")
      ),

      shiny::mainPanel(
        width = 9,
        shiny::tabsetPanel(
          id = "main_tabs",
          shiny::tabPanel("\U0001F4CA Dinámica poblacional", mod_pop_ui("pop")),
          shiny::tabPanel("\U0001F9EC Genética y fenotipo", mod_genetics_ui("gen")),
          shiny::tabPanel("\U0001F5FA Vista espacial", mod_spatial_ui("spatial")),
          shiny::tabPanel("\U0001F4C8 Réplicas (rescate)", mod_replicates_ui("reps"))
        )
      )
    )
  )
}
