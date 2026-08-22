## Tests de simular_hca_replicas:
##  - la ruta paralela foreach/%dopar% debe funcionar en Windows sin
##    adjuntar foreach (error historico 'no se pudo encontrar la
##    funcion "%dopar%"'),
##  - los workers deben recibir TODAS las dependencias del motor,
##  - un `seed` pasado dentro de `...` no debe chocar con la semilla
##    por-replica.

test_that("replicas en paralelo completan y devuelven resumenes", {
  skip_if_not(requireNamespace("foreach", quietly = TRUE) &&
              requireNamespace("doParallel", quietly = TRUE),
              "foreach/doParallel no disponibles")

  res <- simular_hca_replicas(n_replicas = 2, nx = 20, ny = 20,
                               n_steps = 40, paso_introduccion = 20)
  expect_length(res, 2)
  expect_true(all(vapply(res, function(r) is.data.frame(r$resumen), logical(1))))
  expect_true(all(vapply(res, function(r) nrow(r$resumen) > 0, logical(1))))
})

test_that("seed dentro de ... se ignora sin error de argumentos multiples", {
  res <- simular_hca_replicas(n_replicas = 1, nx = 20, ny = 20,
                               n_steps = 30, paso_introduccion = 15,
                               seed = NULL)
  expect_length(res, 1)
  expect_s3_class(res[[1]]$resumen, "data.frame")
})
