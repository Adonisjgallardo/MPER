## ============================================================
## 05_motor.R
## Motor principal: acopla PDE (nutriente + antibiótico) <-> CA
## (crecimiento + muerte con selección estabilizadora)
##
## Modelo de RESCATE EVOLUTIVO espacial: población crece libre
## durante `paso_introduccion` pasos, luego el antibiótico aparece
## de forma ABRUPTA (heterogéneo en espacio vía la plantilla, y
## en el tiempo vía el modo de modulación post-introducción).
## ============================================================

## (Los scripts 01-04 son fuente por app.R antes que este archivo;
## no se usa source() aquí para que la app funcione igual como
## Shiny app standalone o como paquete instalado.)

#' Parametros efectivos del evento de antibiotico de nivel k
#'
#' Escalera de regimenes antibiotic progresivamente mas fuertes y mas
#' exigentes para la funcionalidad de multirresistencia:
#'
#' - Dosis:        `A_max_k      = A_max * factor_dosis^k`
#' - Seleccion:    `omega2_k     = max(omega2_min, omega2 / factor_exigencia^k)`
#'                 (tolerancia fenotipica cada vez mas estrecha)
#' - Optimo:       `theta1_k     = min(1, theta0 + (theta1-theta0) * factor_exigencia^k)`
#'                 (exige resistencia cercana al techo)
#' - Mortalidad:   `mort_estres_k = min(1, mort_estres * factor_exigencia^k)`
#'
#' El nivel 0 reproduce los parametros base tal cual (evento de
#' introduccion / escenario sin multirresistencia).
#'
#' @param nivel entero >= 0; 0 = antibiotico base
#' @param A_max,theta0,theta1,omega2,mort_estres valores base del modelo
#' @param factor_dosis multiplicador geometrico de dosis por nivel (>= 1)
#' @param factor_exigencia endurecimiento geometrico de la seleccion por
#'        nivel (>= 1); 1 deja seleccion/mortalidad/optimos intactos
#' @param omega2_min piso de omega2 para evitar seleccion infinitamente estrecha
#'
#' @return lista con `nivel`, `A_max`, `theta1`, `omega2`, `mort_estres`
parametros_evento_antibiotico <- function(nivel,
                                           A_max, theta0, theta1, omega2,
                                           mort_estres,
                                           factor_dosis = 2,
                                           factor_exigencia = 1.5,
                                           omega2_min = 0.001) {
  k  <- max(0L, as.integer(nivel))
  fd <- max(1, as.numeric(factor_dosis))
  fe <- max(1, as.numeric(factor_exigencia))
  list(
    nivel       = k,
    A_max       = A_max * fd^k,
    theta1      = min(1, theta0 + (theta1 - theta0) * fe^k),
    omega2      = max(omega2_min, omega2 / fe^k),
    mort_estres = min(1, mort_estres * fe^k)
  )
}

#' Ejecuta la simulación HCA de rescate evolutivo completa
#'
#' @param nx,ny dimensiones de la retícula
#' @param n_steps número de pasos de tiempo
#' @param D,k,N0,N_umbral,p_max,N_half igual que en la versión nutriente-pura
#' @param g_inicial,sigma_e parámetros de la colonia fundadora (genética cuantitativa)
#' @param mutacion_sd sd de la mutación aditiva de `g` por división
#' @param paso_introduccion paso en el que ocurre el CAMBIO AMBIENTAL ABRUPTO
#'        (introducción del antibiótico). Antes de este paso, A=0 en toda
#'        la retícula: la población crece libre (fase "silvestre").
#' @param tipo_plantilla heterogeneidad ESPACIAL del antibiótico:
#'        "uniforme", "lineal" (gradiente) o "radial"
#' @param modo_temporal heterogeneidad TEMPORAL post-introducción:
#'        "constante", "pulsos" o "sinusoidal"
#' @param periodo_temporal periodo (en pasos) de la fluctuación temporal
#' @param A_max concentración máxima de antibiótico post-introducción
#' @param multiresistance número de eventos de antibiótico ADICIONALES al
#'        choque inicial (rango típico 0-10). Cada evento posterior introduce
#'        un régimen más fuerte y más exigente según `factor_dosis` y
#'        `factor_exigencia`; se disparan cuando la población recupera su
#'        pico previo o al vencer `tope_espera_shock`, lo que ocurra primero
#' @param factor_dosis multiplicador geométrico de concentración (`A_max`)
#'        aplicado en cada evento adicional (1 = misma dosis siempre)
#' @param factor_exigencia endurecimiento geométrico de la selección en cada
#'        evento adicional: omega2 se divide entre él, theta1 se empuja hacia
#'        1 y mort_estres crece (1 = sin endurecimiento)
#' @param tope_espera_shock pasos máximos de espera desde el último evento
#'        para que la población recupere su pico previo antes de disparar el
#'        siguiente evento de todos modos
#' @param D_A,delta_A,tasa_dosificacion parámetros de la PDE del antibiótico
#' @param theta0,theta1,omega2,K_theta parámetros de selección estabilizadora (Caja 3)
#' @param mort_base,mort_estres,K_mort parámetros de mortalidad
#' @param Nc tamaño poblacional crítico para el análisis de tiempos de
#'        persistencia (Caja 4 / Suplemento S2)
#' @param guardar_cada cada cuántos pasos se guarda un snapshot completo
#' @param seed semilla
#' @param on_step callback opcional `function(paso, n_steps)` invocada al
#'        final de cada paso de tiempo (util para barras de progreso en la
#'        app Shiny); se ignora si es NULL
#'
#' @return lista: `historial` (snapshots), `resumen` (data.frame por paso,
#'         incluye g_bar, var_g, z_bar, A_medio, nacimientos, muertes y
#'         nivel_shock = eventos antibióticos ocurridos hasta ese paso),
#'         `eventos_shock` (data.frame con paso/nivel/parámetros efectivos/
#'         tipo de disparo de cada evento), estado final de todos los campos,
#'         y `paso_extincion` (o NA)
simular_hca <- function(nx = 100, ny = 100,
                         n_steps = 600,
                         D = 0.20, k = 0.15, dt = 1, N0 = 1.0,
                         N_umbral = 0.15, p_max = 0.6, N_half = 0.3,
                         g_inicial = 0.1, sigma_e = 0.05, mutacion_sd = 0.02,
                         paso_introduccion = 150,
                         tipo_plantilla = c("lineal", "uniforme", "radial"),
                         modo_temporal = c("constante", "sinusoidal", "pulsos"),
                         periodo_temporal = 40,
                         A_max = 1.0, multiresistance = 0L,
                         factor_dosis = 2, factor_exigencia = 1.5,
                         tope_espera_shock = 150L, D_A = 0.15,
                         delta_A = 0.02, tasa_dosificacion = 0.3,
                         theta0 = 0.1, theta1 = 0.8, omega2 = 0.05, K_theta = 0.3,
                         mort_base = 0.01, mort_estres = 0.9, K_mort = 0.3,
                         Nc = 200,
                         guardar_cada = 20,
                         guardar_historial_en_disco = FALSE,
                         historial_dir = NULL,
                         seed = 123,
                         on_step = NULL) {

  tipo_plantilla <- match.arg(tipo_plantilla)
  modo_temporal  <- match.arg(modo_temporal)
  verificar_estabilidad_dt(D, dt)
  verificar_estabilidad_dt(D_A, dt)

  set.seed(seed)
  N <- crear_campo_nutriente(nx, ny, N0 = N0)
  A <- matrix(0, nx, ny)
  plantilla_A <- crear_plantilla_antibiotico(nx, ny, tipo = tipo_plantilla)

  colonia <- crear_colonia_inicial(nx, ny, g_inicial = g_inicial,
                                    sigma_e = sigma_e, seed = seed)
  estado <- colonia$estado; g <- colonia$g; e <- colonia$e

  if (guardar_historial_en_disco) {
    if (is.null(historial_dir)) {
      historial_dir <- file.path(getwd(), "historial")
    }
    dir.create(historial_dir, recursive = TRUE, showWarnings = FALSE)
    historial_files <- character(0)
    historial <- list()
    historial_idx <- 0L
  } else {
    max_snapshots <- max(1, ceiling(n_steps / guardar_cada))
    historial <- vector("list", max_snapshots)
    historial_idx <- 0L
    historial_files <- NULL
  }

  resumen <- data.frame(
    paso = integer(n_steps), 
    n_ocupados = integer(n_steps), 
    n_frente = integer(n_steps),
    nacimientos = integer(n_steps), 
    muertes = integer(n_steps),
    N_medio = numeric(n_steps), 
    A_medio = numeric(n_steps),
    g_bar = numeric(n_steps), 
    var_g = numeric(n_steps), 
    z_bar = numeric(n_steps),
    nivel_shock = integer(n_steps)
  )
  paso_extincion <- NA_integer_

  multiresistance <- as.integer(multiresistance)
  multiresistance <- max(0L, min(multiresistance, 50L))
  shock_count <- 0L
  recovery_target <- if (paso_introduccion > 1L) sum(estado) else 0L
  waiting_for_recovery <- FALSE
  paso_ultimo_evento <- NA_integer_
  tope_espera_shock <- max(1L, as.integer(tope_espera_shock))

  ## Registro de eventos: preasignado al maximo posible (primer choque +
  ## `multiresistance` eventos adicionales).
  max_eventos <- multiresistance + 1L
  eventos_shock <- data.frame(
    paso = integer(max_eventos), nivel = integer(max_eventos),
    tipo = character(max_eventos),
    A_max_efectivo = numeric(max_eventos), theta1_efectivo = numeric(max_eventos),
    omega2_efectivo = numeric(max_eventos), mort_estres_efectivo = numeric(max_eventos)
  )
  n_eventos <- 0L

  registrar_evento <- function(paso, tipo) {
    n_eventos <<- n_eventos + 1L
    pars_ev <<- parametros_evento_antibiotico(
      n_eventos - 1L,
      A_max = A_max, theta0 = theta0, theta1 = theta1, omega2 = omega2,
      mort_estres = mort_estres,
      factor_dosis = factor_dosis, factor_exigencia = factor_exigencia
    )
    eventos_shock[n_eventos, ] <<- list(
      paso = paso, nivel = pars_ev$nivel, tipo = tipo,
      A_max_efectivo = pars_ev$A_max, theta1_efectivo = pars_ev$theta1,
      omega2_efectivo = pars_ev$omega2,
      mort_estres_efectivo = pars_ev$mort_estres
    )
    paso_ultimo_evento <<- paso
  }

  for (paso in seq_len(n_steps)) {
    n_ocupados_antes <- sum(estado)

    if (paso == paso_introduccion) {
      shock_count <- 1L
      registrar_evento(paso, tipo = "introduccion")
      recovery_target <- n_ocupados_antes
      waiting_for_recovery <- TRUE   # exigir una caida observada antes de "recuperacion"
    } else if (multiresistance > 0L && shock_count < multiresistance + 1L &&
                !is.na(paso_ultimo_evento)) {
      recuperado   <- n_ocupados_antes >= recovery_target && !waiting_for_recovery
      vencio_tope  <- (paso - paso_ultimo_evento) >= tope_espera_shock
      if (recuperado || vencio_tope) {
        shock_count <- shock_count + 1L
        registrar_evento(paso, tipo = if (recuperado) "recuperacion" else "tope")
        recovery_target <- n_ocupados_antes
        waiting_for_recovery <- TRUE
      } else if (n_ocupados_antes < recovery_target) {
        waiting_for_recovery <- FALSE
      }
    }

    ## Parametros efectivos del antibiotico vigente (nivel = eventos - 1)
    pars_vigente <- parametros_evento_antibiotico(
      max(0L, shock_count - 1L),
      A_max = A_max, theta0 = theta0, theta1 = theta1, omega2 = omega2,
      mort_estres = mort_estres,
      factor_dosis = factor_dosis, factor_exigencia = factor_exigencia
    )

    # ---- Campos continuos ----
    N <- actualizar_nutriente(N, estado, D = D, k = k, dt = dt, N_half = N_half)

    A_obj <- antibiotico_objetivo(plantilla_A, paso, paso_introduccion,
                                   A_max = pars_vigente$A_max,
                                   modo_temporal = modo_temporal,
                                   periodo = periodo_temporal)
    A <- actualizar_antibiotico(A, A_obj, D_A = D_A, delta_A = delta_A, dt = dt,
                                 tasa_dosificacion = tasa_dosificacion)

    # ---- Componente ambiental: se resortea cada paso (Caja 3) ----
    e <- resortear_ambiente(estado, e, sigma_e)

    # ---- CA: crecimiento + muerte con selección estabilizadora ----
    res_ca <- crecer_colonia(estado, g, e, N, A,
                              theta0 = theta0,
                              theta1 = pars_vigente$theta1,
                              omega2 = pars_vigente$omega2,
                              K_theta = K_theta,
                              N_umbral = N_umbral, p_max = p_max, N_half = N_half,
                              mort_base = mort_base,
                              mort_estres = pars_vigente$mort_estres,
                              K_mort = K_mort,
                              mutacion_sd = mutacion_sd, sigma_e = sigma_e)
    estado <- res_ca$estado; g <- res_ca$g; e <- res_ca$e

    n_ocupados <- sum(estado)
    n_frente   <- sum(res_ca$fenotipo == "proliferativo", na.rm = TRUE)
    g_vivos    <- g[estado]
    z_vivos    <- fenotipo_z(g, e)[estado]

    resumen[paso, ] <- list(
      paso = paso, n_ocupados = n_ocupados, n_frente = n_frente,
      nacimientos = res_ca$nacimientos, muertes = res_ca$muertes,
      N_medio = mean(N), A_medio = mean(A),
      g_bar = if (n_ocupados > 0) mean(g_vivos) else NA_real_,
      var_g = if (n_ocupados > 1) var(g_vivos) else NA_real_,
      z_bar = if (n_ocupados > 0) mean(z_vivos) else NA_real_,
      nivel_shock = shock_count
    )

    if (is.na(paso_extincion) && n_ocupados < Nc && paso > paso_introduccion) {
      paso_extincion <- paso   # primer cruce bajo el umbral crítico POST-shock
    }

    if (paso %% guardar_cada == 0 || paso == n_steps) {
      snapshot <- list(
        paso = paso, estado = estado, N = N, A = A,
        resist = fenotipo_z(g, e),   # alias "resist"=z para compatibilidad con 06_analisis.R
        g = g, e = e, fenotipo = res_ca$fenotipo
      )

      if (guardar_historial_en_disco) {
        file_path <- file.path(historial_dir, sprintf("historial_paso_%05d.rds", paso))
        saveRDS(snapshot, file_path)
        historial_files <- c(historial_files, file_path)
      } else {
        historial_idx <- historial_idx + 1L
        historial[[historial_idx]] <- snapshot
      }
    }

    if (!is.null(on_step)) on_step(paso, n_steps)

    if (n_ocupados == 0) break   # extinción total: no tiene sentido seguir
  }

  list(historial = historial, resumen = resumen,
       estado_final = estado, N_final = N, A_final = A,
       g_final = g, e_final = e,
       resist_final = fenotipo_z(g, e),   # alias para 06_analisis.R
       paso_introduccion = paso_introduccion,
       multiresistance = multiresistance,
       shock_count = shock_count,
       eventos_shock = if (n_eventos > 0) {
         eventos_shock[seq_len(n_eventos), , drop = FALSE]
       } else {
         eventos_shock[0, , drop = FALSE]
       },
       paso_extincion = paso_extincion,
       parametros = list(theta0 = theta0, theta1 = theta1, omega2 = omega2,
                          sigma_e = sigma_e, mutacion_sd = mutacion_sd,
                          Nc = Nc, g_inicial = g_inicial,
                          multiresistance = multiresistance,
                          factor_dosis = factor_dosis,
                          factor_exigencia = factor_exigencia,
                          tope_espera_shock = tope_espera_shock))
}

#' Ejecuta varias réplicas independientes (mismos parámetros, semillas
#' distintas) — necesario para construir la distribución EMPÍRICA de
#' tiempos de persistencia/extinción y compararla con las fórmulas
#' analíticas de la Caja 4 / Suplemento S2 (ver 07_genetica_cuantitativa.R)
#'
#' @param n_replicas número de corridas independientes
#' @param ... argumentos pasados a `simular_hca()`
#' @return lista de longitud `n_replicas`, cada elemento un `resumen`
simular_hca_replicas <- function(n_replicas = 20, semilla_base = 1000,
                                  n_cores = parallel::detectCores(logical = FALSE),
                                  use_foreach = FALSE,
                                  ...) {
  if (n_replicas < 1L) return(list())
  n_cores <- max(1L, min(n_cores, n_replicas))

  ## Los puntos se capturan UNA VEZ aqui (no dentro de run_replica): asi
  ## la clausura no referencia `...` de un contexto ajeno (evita el aviso
  ## "... may be used in an incorrect context") y no quedan promesas
  ## perezosas que serializar hacia los workers.
  args_base <- list(...)
  args_base$seed <- NULL   # la semilla por-replica tiene prioridad

  run_replica <- function(i) {
    sim_i <- do.call(simular_hca, c(list(seed = semilla_base + i), args_base))
    list(resumen = sim_i$resumen,
         paso_introduccion = sim_i$paso_introduccion,
         paso_extincion = sim_i$paso_extincion)
  }

  if (.Platform$OS.type != "windows" && !use_foreach) {
    parallel::mclapply(seq_len(n_replicas), run_replica, mc.cores = n_cores)
    } else {
      if (!requireNamespace("foreach", quietly = TRUE) ||
          !requireNamespace("doParallel", quietly = TRUE)) {
        warning("foreach/doParallel not available; running simular_hca_replicas() serially")
        lapply(seq_len(n_replicas), run_replica)
      } else {
        ## El operador infijo %dopar% vive en el namespace de foreach pero
        ## requireNamespace() NO lo adjunta al search path; sin este enlace
        ## local la busqueda del operador falla con
        ## 'no se pudo encontrar la funcion "%dopar%"'.
        `%dopar%` <- foreach::`%dopar%`
        cl <- parallel::makeCluster(n_cores)
        doParallel::registerDoParallel(cl)
        on.exit({
          parallel::stopCluster(cl)
          doParallel::stopImplicitCluster()
        }, add = TRUE)

        ## Exportar a los workers TODAS las funciones del entorno donde
        ## vive el motor (globalenv si se corre por source()/load_all(),
        ## namespace si es el paquete instalado): exportar solo
        ## `simular_hca` deja a los trabajadores sin sus dependencias
        ## internas (p. ej. 'no se pudo encontrar la funcion
        ## "verificar_estabilidad_dt"').
        env_motor <- environment(simular_hca)
        nombres <- ls(env_motor, all.names = FALSE)
        funciones_motor <- nombres[vapply(
          nombres,
          function(n) is.function(get(n, envir = env_motor)),
          logical(1)
        )]
        foreach::foreach(i = seq_len(n_replicas), .export = funciones_motor,
                         .packages = character(0)) %dopar% {
          run_replica(i)
        }
      }
    }
}
