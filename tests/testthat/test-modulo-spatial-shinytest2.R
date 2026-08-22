## Test end-to-end del modulo espacial con navegador real (shinytest2):
## 1) slider con rango completo (1 -> extincion o ultimo paso),
## 2) generacion del video MP4 embebido y servido por HTTP,
## 3) descarga del MP4.

skip_if_not_installed("shinytest2")

sim_e2e <- simular_hca(nx = 30, ny = 30, n_steps = 80,
                        paso_introduccion = 40, guardar_cada = 10, seed = 11)
esperado_max <- if (!is.na(sim_e2e$paso_extincion)) sim_e2e$paso_extincion else max(sim_e2e$resumen$paso)

app_modulo <- shiny::shinyApp(
  ui = shiny::fluidPage(
    shinyjs::useShinyjs(),
    mod_spatial_ui("spatial")
  ),
  server = function(input, output, session) {
    mod_spatial_server("spatial", sim_reactive = shiny::reactive(sim_e2e))
  }
)

drv <- suppressWarnings(shinytest2::AppDriver$new(app_modulo, name = "mper-spatial",
                                                   timeout = 120000))
withr::defer(drv$stop())

leer_slider <- function() {
  drv$get_js("(() => {
    const el = document.querySelector('#spatial-paso_snapshot');
    const plugin = (typeof $ !== 'undefined' && el) ? $(el).data('ionRangeSlider') : null;
    const r = plugin && plugin.result;
    return { found: !!plugin,
             min: r ? Number(r.min) : null,
             max: r ? Number(r.max) : null,
             html: el ? el.outerHTML.slice(0, 200) : null };
  })()")
}

test_that("el slider cubre desde el paso 1 hasta extincion/ultimo paso", {
  drv$wait_for_value(input = "spatial-paso_snapshot", timeout = 30000)
  sl <- leer_slider()
  expect_true(isTRUE(sl$found),
              info = paste("slider sin plugin ionRangeSlider:", sl$html))
  if (isTRUE(sl$found)) {
    expect_equal(sl$min, 1)
    expect_equal(sl$max, esperado_max)
  }
  expect_equal(drv$get_value(input = "spatial-paso_snapshot"), esperado_max)
})

test_that("generar animacion muestra video embebido servido por HTTP", {
  drv$click("spatial-generar_anim")

  ## Espera activa al elemento <video> real (el uiOutput existe desde el
  ## inicio aunque este vacio, asi que no sirve como senal de fin).
  vid <- NULL
  for (i in 1:60) {
    vid <- drv$get_js("(() => {
      const v = document.querySelector('#spatial-animacion_ui video');
      return v ? { ok: true, src: v.getAttribute('src') }
               : { ok: false };
    })()")
    if (isTRUE(vid$ok)) break
    Sys.sleep(2)
  }

  expect_true(isTRUE(vid$ok),
              info = "el elemento <video> nunca aparecio tras generar la animacion")
  skip_if_not(isTRUE(vid$ok), "sin video no se puede seguir")

  expect_match(vid$src, "^mper_anim_[[:alnum:]]+/evolucion[.]mp4$")

  ## El recurso debe responder por HTTP (lo que pide la etiqueta <video>)
  url_base <- sub("/$", "", drv$get_url())
  tmp <- tempfile(fileext = ".mp4")
  on.exit(unlink(tmp), add = TRUE)
  status <- tryCatch({
    suppressWarnings(download.file(paste0(url_base, "/", vid$src), tmp,
                                    mode = "wb", quiet = TRUE))
    0L
  }, error = function(e) 1L)
  expect_equal(status, 0L)
  expect_gt(file.size(tmp), 5000)
})

test_that("descargar MP4 entrega un archivo no vacio", {
  ruta <- drv$get_download("spatial-download_mp4")
  expect_true(file.exists(ruta))
  expect_gt(file.size(ruta), 5000)
})
