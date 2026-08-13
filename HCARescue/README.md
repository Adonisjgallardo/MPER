# HCA Evolutionary Rescue Simulator — Shiny App

Interfaz Shiny para tu modelo de autómata celular híbrido (HCA) de rescate
evolutivo bacteriano bajo estrés antibiótico, construida directamente sobre
tus 7 scripts originales (`01_campos.R` … `07_analisis.R`), sin reescribir
la lógica del modelo.

## 1. Estructura

```
HCARescue/
├── app.R                        # punto de entrada (shiny::runApp() lo usa)
├── R/
│   ├── 01_campos.R               # (tuyo, sin cambios) campos + Laplaciano
│   ├── 02_colonia.R              # (tuyo, sin cambios) siembra de la colonia
│   ├── 03_pde_difusion.R         # (tuyo, sin cambios) PDE nutriente/antibiótico
│   ├── 04_ca_crecimiento.R       # (tuyo) CA de crecimiento/muerte; bloque
│   │                              #        Rcpp opcional comentado (ver abajo)
│   ├── 05_motor.R                # (tuyo) simular_hca()/simular_hca_replicas();
│   │                              #        se agregó un callback opcional
│   │                              #        `on_step()` para la barra de progreso
│   ├── 06_genetica_cuantitativa.R# (tuyo, sin cambios) fórmulas analíticas
│   ├── 07_analisis.R             # (tuyo) fractal + graficación; magick/av
│   │                              #        ahora son opcionales (guardados)
│   ├── utils.R                   # helpers nuevos (tema ggplot, curva KM, etc.)
│   ├── mod_plot_population.R     # módulo Shiny: pestaña "Dinámica poblacional"
│   ├── mod_plot_genetics.R       # módulo Shiny: pestaña "Genética y fenotipo"
│   ├── mod_plot_spatial.R        # módulo Shiny: pestaña "Vista espacial"
│   ├── mod_plot_replicates.R     # módulo Shiny: pestaña "Réplicas (rescate)"
│   ├── app_ui.R                  # UI completa (sidebar de parámetros + tabs)
│   ├── app_server.R              # lógica del servidor (correr, réplicas, descarga)
│   └── run_app.R                 # lanzador opcional: run_hca_app()
└── www/
    └── styles.css                # estilos
```

Ningún archivo tuyo fue reescrito conceptualmente: solo se quitaron los
`source()`/`library()` internos de `05_motor.R` y `07_analisis.R` (porque
`app.R` ahora centraliza la carga de todos los scripts) y se hizo opcional
la parte de Rcpp/magick/av para que la app no falle si esos paquetes no
están instalados.

## 2. Instalación de dependencias

Obligatorias:

```r
install.packages(c("shiny", "ggplot2", "bslib", "shinyjs"))
```

Opcionales (la app detecta su ausencia y se ajusta sola):

```r
install.packages(c("future", "promises"))   # simulaciones no bloqueantes
install.packages(c("magick", "av"))         # exportar animación GIF/MP4
install.packages(c("doParallel", "foreach"))# solo si vas a correr réplicas en Windows
```

- **Sin `future`/`promises`**: el botón "Ejecutar simulación" sigue
  funcionando, pero corre de forma síncrona con una barra de progreso
  nativa de Shiny (la sesión se bloquea mientras corre, como cualquier
  script de R normal).
- **Sin `magick`**: la pestaña "Vista espacial" funciona igual, solo se
  deshabilita la descarga del GIF animado.

## 3. Cómo correrla

Desde una sesión de R con el working directory en la carpeta `HCARescue/`:

```r
shiny::runApp(".")
```

o, si abres `app.R` en RStudio, usa el botón **Run App**.

También puedes lanzarla desde fuera de la carpeta:

```r
shiny::runApp("ruta/a/HCARescue")
```

## 4. Qué hace cada pestaña

| Pestaña | Contenido |
|---|---|
| 📊 Dinámica poblacional | `N_total` y células proliferativas vs. tiempo, línea vertical en el choque de antibiótico, y comparación contra `Nt_analitico_qg()` (Caja 3). |
| 🧬 Genética y fenotipo | `g_bar` y `z_bar` vs. tiempo con líneas de referencia en `theta0`/`theta1`, y comparación contra `gbar_analitico_qg()` (recursión de Lande). |
| 🗺️ Vista espacial | Snapshot coloreado por resistencia o nutriente (estilo `geom_point` + viridis "plasma"), slider para recorrer el historial si lo guardaste, dimensión fractal (box-counting), y descarga de animación GIF. |
| 📈 Réplicas (rescate) | Corre `simular_hca_replicas()`, calcula tiempos de persistencia empíricos y dibuja una curva de supervivencia estilo Kaplan-Meier, con estadísticas de fracción rescatada. |

El botón **"Guardar historial espacial completo"** en el panel lateral
controla la política de memoria descrita en el blueprint original: si
está desmarcado, `simular_hca()` solo guarda el estado final (rápido,
poca memoria); si lo marcas, guarda un snapshot cada N pasos para poder
animar, a costa de más memoria para retículas grandes / muchos pasos.

## 5. Notas honestas sobre esta entrega

Este entorno de generación de archivos no tiene R instalado ni acceso a
CRAN, así que **no pude ejecutar la app aquí para probarla end-to-end**.
El código se escribió con cuidado (nombres de argumentos verificados
contra tus 7 scripts, un solo `seed` por llamada, IDs de módulos
consistentes entre UI/servidor, etc.), pero te recomiendo:

1. Correrla primero con parámetros pequeños (`nx`/`ny` ~ 50, `n_steps` ~
   150) para confirmar que todo carga bien.
2. Si `future_promise()` te da problemas de serialización de globals al
   correr en `multisession`, la solución más simple es no instalar
   `future`/`promises` — la app usa automáticamente el modo síncrono.

## 6. Convertir esto en un paquete R instalable (opcional)

Si más adelante quieres el paquete formal completo (`DESCRIPTION`,
`NAMESPACE`, `devtools::check()`, `install_github()`, como describe el
blueprint), el siguiente paso es mover `R/*.R` tal cual a `R/` de un
paquete, mover `www/` a `inst/app/www/`, agregar cabeceras `roxygen2`
(`@export`) a `hca_ui()`, `hca_server()`, `run_hca_app()`, y correr
`devtools::document()` + `devtools::check()`. Puedo ayudarte con eso
cuando quieras — es un paso natural una vez que hayas confirmado que la
app funciona como la tienes ahora.
