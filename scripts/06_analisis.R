## ============================================================
## 06_analisis.R
## Análisis morfológico: dimensión fractal (box-counting)
## y visualización básica de la colonia / campos
## ============================================================

#' Dimensión fractal por conteo de cajas (box-counting)
#'
#' @param estado matriz binaria de ocupación
#' @param tamanos vector de tamaños de caja a probar (potencias de 2 recomendadas)
#' @return lista con el ajuste lineal log(N_cajas) ~ log(1/tamano) y la pendiente (dimensión)
dimension_fractal <- function(estado, tamanos = c(2, 4, 5, 8, 10, 16, 20)) {
  nx <- nrow(estado); ny <- ncol(estado)
  conteos <- sapply(tamanos, function(s) {
    nfx <- ceiling(nx / s); nfy <- ceiling(ny / s)
    cajas_ocupadas <- 0L
    for (i in seq_len(nfx)) {
      for (j in seq_len(nfy)) {
        f0 <- (i - 1) * s + 1; f1 <- min(i * s, nx)
        c0 <- (j - 1) * s + 1; c1 <- min(j * s, ny)
        if (any(estado[f0:f1, c0:c1] == 1L)) cajas_ocupadas <- cajas_ocupadas + 1L
      }
    }
    cajas_ocupadas
  })

  df <- data.frame(log_inv_s = log(1 / tamanos), log_N = log(conteos))
  ajuste <- lm(log_N ~ log_inv_s, data = df)

  list(tabla = data.frame(tamano = tamanos, N_cajas = conteos),
       ajuste = ajuste,
       dimension = unname(coef(ajuste)[2]))
}

## ------------------------------------------------------------
## Visualización estilo ggplot/gganimate
## (equivalente al script de referencia con geom_point +
##  scale_colour_viridis_c(option="plasma") + transition_time)
## ------------------------------------------------------------

library(ggplot2)
library(magick)
library(av)

#' Convierte un snapshot (estado/N/resist) en un data.frame "largo"
#' con una fila por sitio OCUPADO, listo para ggplot
#'
#' @param estado matriz 0/1 de ocupación
#' @param N campo de nutriente (misma dimensión que estado)
#' @param resist matriz de resistencia por celda (NA si vacío)
#' @param paso entero: paso de tiempo / "generación" del snapshot
#' @return data.frame con columnas x, y, N, resist, paso
snapshot_a_df <- function(estado, N, resist, paso = NA_integer_) {
  idx <- which(estado == 1L, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(data.frame(x = integer(0), y = integer(0),
                       N = numeric(0), resist = numeric(0),
                       paso = integer(0)))
  }
  data.frame(
    x      = idx[, 1],
    y      = idx[, 2],
    N      = N[idx],
    resist = resist[idx],
    paso   = paso
  )
}

#' Grafica un único estado de la colonia al estilo geom_point + viridis
#'
#' Reemplaza el antiguo panel base-R por un único scatter estilo
#' ggplot, coloreado por la variable elegida (por defecto `resist`,
#' el rasgo fenotípico por-celda), replicando el estilo visual del
#' script de referencia (geom_point + scale_colour_viridis_c("plasma")
#' + coord_fixed + theme_minimal).
#'
#' @param estado,N,resist matrices del estado actual de la simulación
#' @param archivo ruta de salida (.png)
#' @param titulo texto para el título del gráfico
#' @param variable cuál columna colorear: "resist" o "N"
#' @param point_size,alpha estética de los puntos
#' @param width,height dimensiones en píxeles del PNG de salida
graficar_estado_hca <- function(estado, N, resist, archivo,
                                 titulo = "", variable = c("resist", "N"),
                                 point_size = 0.6, alpha = 0.7,
                                 width = 600, height = 600) {
  variable <- match.arg(variable)
  df <- snapshot_a_df(estado, N, resist)

  etiqueta_leyenda <- if (variable == "resist") "Resistencia" else "Nutriente N"

  p <- ggplot(df, aes(x = x, y = y, colour = .data[[variable]])) +
    geom_point(size = point_size, alpha = alpha) +
    scale_colour_viridis_c(option = "plasma", name = etiqueta_leyenda) +
    coord_fixed() +
    labs(title = titulo, x = "Posición X", y = "Posición Y") +
    theme_minimal()

  ggsave(archivo, plot = p, width = width, height = height,
         units = "px", dpi = 150, bg = "white")
  invisible(p)
}

#' Anima la evolución completa de la colonia a partir del `historial`
#' generado por `simular_hca()`, coloreando por resistencia (o nutriente).
#' Mismo lenguaje visual que el script de referencia (geom_point +
#' scale_colour_viridis_c("plasma") + coord_fixed + theme_minimal +
#' título "Paso: {n}" tipo frame_time), pero renderizado cuadro-a-cuadro
#' con magick/av en lugar de gganimate::animate().
#'
#' Nota de implementación: en este entorno no hay acceso a CRAN, así que
#' `gifski` (el renderer estándar de gganimate) no está disponible, y la
#' versión de `magick` instalada aquí es incompatible con
#' `gganimate::magick_renderer()` (bug conocido de compatibilidad de
#' versiones). Por eso se genera cada cuadro como PNG independiente con
#' ggplot2 y se ensamblan directamente con `magick::image_animate()`
#' (GIF) y `av::av_encode_video()` (MP4) — mismo resultado visual, sin
#' depender de esas piezas.
#'
#' @param historial lista de snapshots: cada elemento con
#'        `paso`, `estado`, `N`, `resist` (tal como los guarda `simular_hca()`)
#' @param archivo_gif ruta de salida .gif (NULL para omitir)
#' @param archivo_mp4 ruta de salida .mp4 (NULL para omitir)
#' @param variable "resist" o "N": qué variable colorear
#' @param fps cuadros por segundo
#' @param width,height dimensiones en píxeles de la animación
#' @param point_size,alpha estética de los puntos
#'
#' @return invisible(NULL)
animar_historial_hca <- function(historial, archivo_gif = "figuras/evolucion.gif",
                                  archivo_mp4 = "figuras/evolucion.mp4",
                                  variable = c("resist", "N"),
                                  fps = 10, width = 600, height = 600,
                                  point_size = 0.6, alpha = 0.7) {
  variable <- match.arg(variable)

  df_all <- do.call(rbind, lapply(historial, function(s) {
    snapshot_a_df(s$estado, s$N, s$resist, paso = s$paso)
  }))
  if (nrow(df_all) == 0) stop("El historial no contiene sitios ocupados.")

  etiqueta_leyenda <- if (variable == "resist") "Resistencia" else "Nutriente N"
  limites_color <- range(df_all[[variable]], na.rm = TRUE)
  limites_x <- range(df_all$x); limites_y <- range(df_all$y)

  dir_temp <- tempfile("frames_")
  dir.create(dir_temp)
  on.exit(unlink(dir_temp, recursive = TRUE), add = TRUE)

  archivos_frame <- character(length(historial))

  for (i in seq_along(historial)) {
    s <- historial[[i]]
    df_i <- snapshot_a_df(s$estado, s$N, s$resist, paso = s$paso)

    p_i <- ggplot(df_i, aes(x = x, y = y, colour = .data[[variable]])) +
      geom_point(size = point_size, alpha = alpha) +
      scale_colour_viridis_c(option = "plasma", name = etiqueta_leyenda,
                              limits = limites_color) +
      coord_fixed(xlim = limites_x, ylim = limites_y) +
      labs(title = paste("Paso:", s$paso), x = "Posición X", y = "Posición Y") +
      theme_minimal()

    archivo_i <- file.path(dir_temp, sprintf("frame_%04d.png", i))
    ggsave(archivo_i, plot = p_i, width = width, height = height,
           units = "px", dpi = 150, bg = "white")
    archivos_frame[i] <- archivo_i
  }

  if (!is.null(archivo_gif)) {
    imgs <- magick::image_read(archivos_frame)
    # magick expresa el tiempo en centésimas de segundo (`delay`);
    # se usa en vez de `fps` porque ImageMagick exige que fps sea
    # un divisor exacto de 100.
    delay_centesimas <- max(1, round(100 / fps))
    anim <- magick::image_animate(imgs, delay = delay_centesimas, loop = 0)
    magick::image_write(anim, archivo_gif)
  }

  if (!is.null(archivo_mp4)) {
    av::av_encode_video(archivos_frame, output = archivo_mp4,
                         framerate = fps, verbose = FALSE)
  }

  invisible(NULL)
}
