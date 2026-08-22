## Tests de la escalera de eventos antibiotic (multirresistencia):
## helper puro, compatibilidad del escenario base, disparo por
## recuperacion vs. por tope de espera, y registro de eventos.

test_that("parametros_evento_antibiotico: nivel 0 es identidad", {
  p0 <- parametros_evento_antibiotico(0, A_max = 1, theta0 = 0.1,
                                       theta1 = 0.8, omega2 = 0.05,
                                       mort_estres = 0.9)
  expect_equal(p0$nivel, 0L)
  expect_equal(p0$A_max, 1)
  expect_equal(p0$theta1, 0.8)
  expect_equal(p0$omega2, 0.05)
  expect_equal(p0$mort_estres, 0.9)
})

test_that("parametros_evento_antibiotico: escalado monotono con topes", {
  ps <- lapply(0:5, function(k) {
    parametros_evento_antibiotico(k, A_max = 1, theta0 = 0.1, theta1 = 0.8,
                                   omega2 = 0.05, mort_estres = 0.6,
                                   factor_dosis = 2, factor_exigencia = 1.5)
  })
  amax <- vapply(ps, function(p) p$A_max, numeric(1))
  ome  <- vapply(ps, function(p) p$omega2, numeric(1))
  th   <- vapply(ps, function(p) p$theta1, numeric(1))
  me   <- vapply(ps, function(p) p$mort_estres, numeric(1))

  expect_equal(amax, 2^(0:5))                    # dosis x2 por nivel (geometrica)
  expect_true(all(diff(ome) < 0))                # seleccion cada vez mas estrecha
  expect_true(all(ome >= 0.001))                 # piso de omega2
  expect_true(all(th <= 1))                      # techo del optimo
  expect_equal(th[2], 1)                         # 0.1 + 0.7*1.5 = 1.15 -> 1
  expect_true(all(me <= 1))                      # mortalidad acotada
  expect_true(all(diff(me) >= 0))
  expect_equal(me[3], 1)                         # 0.6 * 1.5^2 = 1.35 -> 1
  ## Niveles negativos se tratan como 0
  expect_equal(parametros_evento_antibiotico(-3, A_max = 1, theta0 = 0.1,
                                              theta1 = 0.8, omega2 = 0.05,
                                              mort_estres = 0.9)$nivel, 0L)
})

args_base_mr <- list(nx = 20, ny = 20, n_steps = 80, paso_introduccion = 30, seed = 42)

test_that("multiresistance=0: un evento con parametros base e intactos por factores", {
  s <- do.call(simular_hca, args_base_mr)
  ev <- s$eventos_shock

  expect_equal(nrow(ev), 1)
  expect_equal(ev$paso[1], 30L)
  expect_equal(ev$nivel[1], 0L)
  expect_equal(ev$tipo[1], "introduccion")
  expect_equal(ev$A_max_efectivo[[1]], 1.0)
  expect_equal(ev$theta1_efectivo[[1]], 0.8)
  expect_equal(ev$omega2_efectivo[[1]], 0.05)
  expect_equal(ev$mort_estres_efectivo[[1]], 0.9)

  expect_true(all(s$resumen$nivel_shock[s$resumen$paso < 30] == 0L))
  expect_true(all(s$resumen$nivel_shock[s$resumen$paso >= 30] == 1L))

  ## Los factores no deben afectar para nada cuando solo hay nivel 0
  s2 <- do.call(simular_hca, c(args_base_mr,
                                list(factor_dosis = 4, factor_exigencia = 3,
                                     tope_espera_shock = 10)))
  expect_identical(s$resumen, s2$resumen)
})

test_that("multiresistance=1 recuperable: segundo evento mas fuerte por recuperacion", {
  ## Estres duro + adaptacion rapida + margen de crecimiento: la poblacion
  ## cae por debajo de su pico pre-shock y lo supera de nuevo (paso 27),
  ## disparando el segundo evento ANTES de que venza el tope (500 > 250).
  s <- simular_hca(nx = 40, ny = 40, n_steps = 250, N0 = 2,
                    paso_introduccion = 20, A_max = 0.8,
                    modo_temporal = "constante",
                    theta1 = 0.85, omega2 = 0.06, K_theta = 0.4,
                    mutacion_sd = 0.15, sigma_e = 0.04,
                    mort_estres = 0.7, Nc = 100,
                    multiresistance = 1, factor_dosis = 2,
                    factor_exigencia = 1.5, tope_espera_shock = 500,
                    seed = 1)
  ev <- s$eventos_shock
  expect_equal(nrow(ev), 2)
  expect_equal(ev$tipo[2], "recuperacion")
  expect_equal(ev$nivel[2], 1L)
  expect_equal(ev$A_max_efectivo[2], 2 * ev$A_max_efectivo[1])
  expect_gt(ev$theta1_efectivo[2], ev$theta1_efectivo[1])
  expect_lt(ev$omega2_efectivo[2], ev$omega2_efectivo[1])
  expect_gt(ev$mort_estres_efectivo[2], ev$mort_estres_efectivo[1])
  expect_equal(max(s$resumen$nivel_shock), 2L)
})

test_that("multiresistance=2 letal: eventos disparados por tope de espera", {
  ## Poblacion sobrevive al shock sin recuperar el pico: la escalera avanza
  ## cada `tope_espera_shock` pasos (pulsos garantizados 20 -> 45 -> 70).
  s <- simular_hca(nx = 40, ny = 40, n_steps = 120, N0 = 2,
                    paso_introduccion = 20, A_max = 0.6,
                    theta1 = 0.8, omega2 = 0.08, K_theta = 0.3,
                    mutacion_sd = 0.03, sigma_e = 0.05,
                    mort_estres = 0.5, Nc = 400,
                    multiresistance = 2, factor_dosis = 2,
                    factor_exigencia = 1.5, tope_espera_shock = 25,
                    seed = 11)
  ev <- s$eventos_shock
  expect_equal(nrow(ev), 3)
  expect_equal(ev$paso, c(20L, 45L, 70L))
  expect_equal(ev$tipo[-1], c("tope", "tope"))
  expect_true(all(diff(ev$nivel) == 1L))         # niveles consecutivos
  expect_equal(ev$A_max_efectivo, c(0.6, 1.2, 2.4))
  expect_equal(max(s$resumen$nivel_shock), 3L)
  expect_equal(s$shock_count, 3L)
  expect_equal(ev$theta1_efectivo[2], 1)         # techo del optimo alcanzado
})

test_that("extincion temprana: no hay eventos posteriores fantasma", {
  s <- simular_hca(nx = 20, ny = 20, n_steps = 200, N0 = 1,
                    paso_introduccion = 25, A_max = 5,
                    theta1 = 0.99, omega2 = 0.01, K_theta = 0.3,
                    mutacion_sd = 0.005, sigma_e = 0.05,
                    mort_estres = 1, Nc = 350,
                    multiresistance = 5, factor_dosis = 2,
                    factor_exigencia = 2, tope_espera_shock = 1000,
                    seed = 3)
  ev <- s$eventos_shock
  expect_lte(nrow(ev), 2)                        # sin recuperacion ni tope alcanzable
  expect_false(is.na(s$paso_extincion))
  expect_true(all(s$resumen$nivel_shock <= nrow(ev)))
})

test_that("recolectar_parametros_sim transmite los nuevos inputs", {
  inp <- list(nx = 20, n_steps = 100, D = 0.2, k = 0.15, N0 = 1,
              N_umbral = 0.15, p_max = 0.6, N_half = 0.3,
              g_inicial = 0.1, sigma_e = 0.05, mutacion_sd = 0.02,
              paso_introduccion = 30, tipo_plantilla = "lineal",
              modo_temporal = "constante", periodo_temporal = 40,
              A_max = 1, multiresistance = 2,
              factor_dosis = 3, factor_exigencia = 2.5,
              tope_espera_shock = 77,
              D_A = 0.15, delta_A = 0.02, tasa_dosificacion = 0.3,
              theta0 = 0.1, theta1 = 0.8, omega2 = 0.05, K_theta = 0.3,
              mort_base = 0.01, mort_estres = 0.9, K_mort = 0.3, Nc = 200,
              mantener_historial = FALSE)
  p <- recolectar_parametros_sim(inp)
  expect_equal(p$factor_dosis, 3)
  expect_equal(p$factor_exigencia, 2.5)
  expect_equal(p$tope_espera_shock, 77)
})
