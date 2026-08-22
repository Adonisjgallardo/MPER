## Tests del grafico estatico y de la animacion con ejes fijos

sim_peq <- simular_hca(nx = 30, ny = 30, n_steps = 60,
                        paso_introduccion = 30, guardar_cada = 10, seed = 7)
hist_peq <- Filter(Negate(is.null), sim_peq$historial)

extraer_rangos <- function(p) {
  pp <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]
  list(x = pp$x$continuous_range, y = pp$y$continuous_range)
}

test_that("plot_snapshot_estatico fija los ejes a la reticula completa", {
  p_ini <- plot_snapshot_estatico(hist_peq[[1]], "resist")
  p_fin <- plot_snapshot_estatico(hist_peq[[length(hist_peq)]], "resist")

  r_ini <- extraer_rangos(p_ini)
  r_fin <- extraer_rangos(p_fin)

  expect_equal(r_ini$x, c(0.5, 30.5))
  expect_equal(r_ini$y, c(0.5, 30.5))
  ## Estaticos: identicos entre el primer y ultimo snapshot
  expect_equal(r_fin$x, r_ini$x)
  expect_equal(r_fin$y, r_ini$y)
})

test_that("plot_snapshot_estatico no falla con colonia extinta (sin sitios)", {
  snap_vacio <- list(
    paso   = 99L,
    estado = matrix(FALSE, 30, 30),
    N      = matrix(1, 30, 30),
    resist = matrix(NA_real_, 30, 30)
  )
  p <- plot_snapshot_estatico(snap_vacio, "resist")
  r <- extraer_rangos(p)
  expect_equal(r$x, c(0.5, 30.5))
  expect_equal(r$y, c(0.5, 30.5))
})

test_that("animar_historial_hca acepta nx/ny y produce MP4", {
  mp4 <- tempfile(fileext = ".mp4")
  on.exit(unlink(mp4), add = TRUE)
  suppressWarnings(
    animar_historial_hca(hist_peq, archivo_gif = NULL, archivo_mp4 = mp4,
                          mostrar_fractal = FALSE, nx = 30, ny = 30,
                          width = 300, height = 300)
  )
  expect_true(file.exists(mp4))
  expect_gt(file.size(mp4), 1000)
})

test_that("animar_historial_hca mantiene compatibilidad sin nx/ny (GIF)", {
  gif <- tempfile(fileext = ".gif")
  on.exit(unlink(gif), add = TRUE)
  skip_if_not(.magick_disponible, "magick no disponible")
  suppressWarnings(
    animar_historial_hca(hist_peq, archivo_gif = gif, archivo_mp4 = NULL,
                          mostrar_fractal = FALSE, width = 300, height = 300)
  )
  expect_true(file.exists(gif))
  expect_gt(file.size(gif), 1000)
})
