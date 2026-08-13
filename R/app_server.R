## ============================================================
## app_server.R
## Logica del servidor de HCARescue.
##
## Ejecucion no bloqueante: si los paquetes `future` y `promises`
## estan instalados, la simulacion corre en un proceso en segundo
## plano (future::plan(multisession)) y la UI permanece receptiva
## mientras corre. Si no estan instalados, se hace un fallback
## sincrono con una barra de progreso nativa de Shiny (la app sigue
## siendo utilizable, solo que bloquea la sesion mientras corre).
## ============================================================

.async_disponible <- requireNamespace("future", quietly = TRUE) &&
                      requireNamespace("promises", quietly = TRUE)

if (.async_disponible) {
  future::plan(future::multisession)
}

hca_server <- function(input, output, session) {

  rv <- shiny::reactiveValues(
    sim_result   = NULL,
    reps_result  = NULL,
    is_running   = FALSE,
    reps_running = FALSE
  )

  set_status <- function(msg) {
    shinyjs::html("status_msg", msg)
  }

  ## ---- 1) Simulacion unica ------------------------------------
  shiny::observeEvent(input$run, {
    if (rv$is_running) return(invisible(NULL))
    rv$is_running <- TRUE
    shinyjs::disable("run"); shinyjs::disable("run_reps")
    set_status("Ejecutando simulacion...")

    params <- recolectar_parametros_sim(input)

    if (.async_disponible) {
      prom <- promises::future_promise({
        do.call(simular_hca, params)
      }, seed = TRUE)

      promises::then(prom,
        onFulfilled = function(result) {
          rv$sim_result <- result
          rv$is_running <- FALSE
          shinyjs::enable("run"); shinyjs::enable("run_reps")
          set_status("Simulacion completada.")
          shiny::showNotification("Simulacion completada.", type = "message")
        },
        onRejected = function(error) {
          rv$is_running <- FALSE
          shinyjs::enable("run"); shinyjs::enable("run_reps")
          set_status("Error en la simulacion.")
          shiny::showNotification(paste("Error:", conditionMessage(error)), type = "error", duration = NULL)
        }
      )
      NULL
    } else {
      shiny::withProgress(message = "Ejecutando simulacion", value = 0, {
        n_steps <- params$n_steps
        params$on_step <- function(paso, n_steps_total) {
          if (paso %% max(1L, floor(n_steps_total / 50)) == 0) {
            shiny::setProgress(value = paso / n_steps_total,
                                detail = sprintf("paso %d / %d", paso, n_steps_total))
          }
        }
        resultado <- tryCatch(do.call(simular_hca, params),
                               error = function(e) e)
        if (inherits(resultado, "error")) {
          shiny::showNotification(paste("Error:", conditionMessage(resultado)), type = "error", duration = NULL)
          set_status("Error en la simulacion.")
        } else {
          rv$sim_result <- resultado
          set_status("Simulacion completada.")
          shiny::showNotification("Simulacion completada.", type = "message")
        }
      })
      rv$is_running <- FALSE
      shinyjs::enable("run"); shinyjs::enable("run_reps")
    }
  })

  ## ---- 2) Replicas ---------------------------------------------
  shiny::observeEvent(input$run_reps, {
    if (rv$reps_running) return(invisible(NULL))
    rv$reps_running <- TRUE
    shinyjs::disable("run"); shinyjs::disable("run_reps")
    set_status(sprintf("Ejecutando %d replicas...", input$n_reps))

    params <- recolectar_parametros_sim(input)
    params$guardar_cada <- params$n_steps + 1L   # nunca se necesita historial espacial en replicas
    params$seed <- NULL   # simular_hca_replicas() asigna su propia semilla por replica
    n_reps <- input$n_reps

    if (.async_disponible) {
      prom <- promises::future_promise({
        do.call(simular_hca_replicas, c(list(n_replicas = n_reps), params))
      }, seed = TRUE)

      promises::then(prom,
        onFulfilled = function(result) {
          rv$reps_result <- result
          rv$reps_running <- FALSE
          shinyjs::enable("run"); shinyjs::enable("run_reps")
          set_status(sprintf("%d replicas completadas.", n_reps))
          shiny::showNotification("Replicas completadas.", type = "message")
        },
        onRejected = function(error) {
          rv$reps_running <- FALSE
          shinyjs::enable("run"); shinyjs::enable("run_reps")
          set_status("Error en las replicas.")
          shiny::showNotification(paste("Error:", conditionMessage(error)), type = "error", duration = NULL)
        }
      )
      NULL
    } else {
      shiny::withProgress(message = sprintf("Ejecutando %d replicas (esto puede tardar)", n_reps), value = 0.3, {
        resultado <- tryCatch(
          do.call(simular_hca_replicas, c(list(n_replicas = n_reps), params)),
          error = function(e) e
        )
        shiny::setProgress(value = 1)
        if (inherits(resultado, "error")) {
          shiny::showNotification(paste("Error:", conditionMessage(resultado)), type = "error", duration = NULL)
          set_status("Error en las replicas.")
        } else {
          rv$reps_result <- resultado
          set_status(sprintf("%d replicas completadas.", n_reps))
          shiny::showNotification("Replicas completadas.", type = "message")
        }
      })
      rv$reps_running <- FALSE
      shinyjs::enable("run"); shinyjs::enable("run_reps")
    }
  })

  ## ---- 3) Modulos de graficos ------------------------------------
  mod_pop_server("pop", sim_reactive = shiny::reactive(rv$sim_result),
                  sigma_g2 = shiny::reactive(input$sigma_g2),
                  sigma_e2 = shiny::reactive(input$sigma_e^2),
                  Wmax     = shiny::reactive(input$Wmax))

  mod_genetics_server("gen", sim_reactive = shiny::reactive(rv$sim_result),
                       sigma_g2 = shiny::reactive(input$sigma_g2),
                       sigma_e2 = shiny::reactive(input$sigma_e^2))

  mod_spatial_server("spatial", sim_reactive = shiny::reactive(rv$sim_result))

  mod_replicates_server("reps", reps_reactive = shiny::reactive(rv$reps_result))

  ## ---- 4) Descarga de resultados ------------------------------------
  output$download_results <- shiny::downloadHandler(
    filename = function() paste0("hca_resultados_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"),
    content = function(file) {
      saveRDS(list(simulacion = rv$sim_result, replicas = rv$reps_result,
                    parametros = recolectar_parametros_sim(input)), file)
    }
  )
}
