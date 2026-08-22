# MPER

**Adaptive Reorganization and Ecological Perturbation Model**

MPER is a spatial, agent-based model for studying how microbial populations respond
to abrupt environmental change (an antibiotic shock) and whether they achieve
*evolutionary rescue*. It couples:

- two reaction-diffusion (PDE) fields — nutrient `N(x,y,t)` and antibiotic `A(x,y,t)` —
- a cellular automaton (CA) layer of sites that get colonized, divide, and die, and
- quantitative genetics: each cell carries a heritable genetic value `g`, a non-heritable
  environmental deviation `e`, and phenotype `z = g + e`, under stabilizing selection
  around a local optimum that tracks the antibiotic concentration.

The repository is organized as an R package so that the simulation engine, the Shiny
app, and the experiment scripts can be installed and reused together.

## Installation

From a source checkout:

```r
# development mode (no installation needed)
devtools::load_all(".")

# or install it
install.packages(".", repos = NULL, type = "source")
```

Requirements are listed in [DESCRIPTION](DESCRIPTION); `av` (MP4 encoding) and
`magick` (GIF encoding) are optional extras used for animations.

## Proposed workflow

The package supports five complementary workflows, from quick visual exploration to
statistically designed experiments.

### 1. Interactive exploration (Shiny app)

Launch the app:

```r
MPER::run_mper()          # if installed
shiny::runApp("inst/app") # from a source checkout
```

The proposed interactive loop is:

1. Set the shock scenario in the sidebar (grid size, steps, `A_max`,
   `paso_introduccion`, selection curve `theta0`/`theta1`/`omega2`, mortality, ...).
2. Optionally escalate the stress with additional antibiotic events
   (`multiresistance > 0` reveals three extra controls: dose factor, exigency factor,
   and the maximum wait between events).
3. Click *Ejecutar simulación* and inspect the result across tabs:
   - **Population dynamics**: occupied sites and proliferative front over time,
     with vertical markers for every antibiotic event and an overlay of the
     analytical prediction (Caja 3).
   - **Spatial view**: final state, fractal dimension, a snapshot slider that covers
     the whole run (step 1 through extinction or last step) with FIXED axes, and an
     on-demand MP4 animation embedded in the page (with download).
   - **Genetics**: evolution of mean genetic value `g_bar` and variance against the
     analytical recursion.
4. Use *Ejecutar réplicas* to repeat the same scenario several times inside the app
   and compare persistence outcomes.

### 2. Single programmatic simulation

For reproducible work, drive the engine directly:

```r
sim <- simular_hca(
  nx = 30, ny = 30, n_steps = 200,
  N0 = 1, paso_introduccion = 40, A_max = 1.0,
  theta0 = 0.10, theta1 = 0.80, omega2 = 0.05,
  mutacion_sd = 0.02, mort_estres = 0.90,
  Nc = 200, guardar_cada = 5,
  seed = 42
)

sim$paso_extincion     # first step below critical threshold Nc (NA if rescued)
head(sim$resumen)      # per-step summary: occupancy, front, births/deaths,
                       # N_medio, A_medio, g_bar, var_g, z_bar, nivel_shock
sim$eventos_shock      # table of antibiotic events (see below)
```

Then post-process:

```r
# static snapshot with fixed axes over the full lattice
plot_snapshot_estatico(sim$historial[[length(sim$historial)]], nx = 30, ny = 30)

# animation of the whole run (GIF requires magick, MP4 requires av)
animar_historial_hca(sim$historial, archivo_mp4 = "figuras/evolucion.mp4",
                     nx = 30, ny = 30)

# compare simulated trajectories with the analytical formulas (Caja 3)
cmp_N <- comparar_Nt_caja3(sim, Wmax = 1, sigma_g2 = 0.002, sigma_e2 = 0.0025)
cmp_g <- comparar_gbar_caja3(sim, sigma_g2 = 0.002, sigma_e2 = 0.0025)

# fractal dimension of the final colony
calcular_dimension_fractal(sim$estado_final)$dimension
```

### 3. Replicates and persistence statistics

Single runs are deterministic given a seed; rescue is inherently stochastic, so the
proposed workflow repeats each scenario over independent seeds:

```r
replicas <- simular_hca_replicas(
  50,                              # n_replicas: independent runs
  semilla_base = 1000,             # replica i uses seed 1000 + i
  nx = 30, ny = 30, n_steps = 200,
  paso_introduccion = 40, A_max = 1.0
)

# each element carries what persistence analysis needs
str(replicas[[1]])                 # list(resumen, paso_introduccion, paso_extincion)

# empirical distribution of persistence/extinction times
tiempos <- tiempos_persistencia_empiricos(replicas)
summary(tiempos)
```

Replicates run in parallel by default — `mclapply` forking on Unix and a
`foreach`/`doParallel` cluster on Windows — falling back to serial execution with a
warning if the parallel backends are unavailable.

### 4. Escalating-multiresistance experiments

`multiresistance = k` schedules `k` ADDITIONAL antibiotic events after the initial
shock. Each level escalates geometrically:

- dose: `A_max_k = A_max * factor_dosis^k`
- selection narrows: `omega2_k = max(omega2_min, omega2 / factor_exigencia^k)`
- optimum hardens: `theta1_k = min(1, theta0 + (theta1 - theta0) * factor_exigencia^k)`
- stress mortality rises: `mort_estres_k = min(1, mort_estres * factor_exigencia^k)`

Events fire when the population recovers its previous peak **or** when
`tope_espera_shock` steps elapse since the last event, so the ladder always completes
even if the population never recovers. Every event is logged in `eventos_shock`
(step, level, trigger type `introduccion`/`recuperacion`/`tope`, effective parameters),
and `resumen$nivel_shock` records which antibiotic regime was active at every step.
Level 0 is bit-for-bit identical to a plain single-shock run given the same seed.

```r
sim_mr <- simular_hca(nx = 40, ny = 40, n_steps = 250,
                       paso_introduccion = 20, A_max = 0.6,
                       multiresistance = 2, factor_dosis = 2,
                       factor_exigencia = 1.5, tope_espera_shock = 25,
                       seed = 11)
sim_mr$eventos_shock
```

### 5. Designed experiments

For systematic exploration of parameter space, use the ready-made experiment scripts
under [example/scripts](example/scripts):

- [example/scripts/08_simulate.R](example/scripts/08_simulate.R): canonical single-run
  workflow; exports data and plots to [example/data/simulation](example/data/simulation)
  and [example/plots/simulation](example/plots/simulation).
- [example/scripts/09_profile.R](example/scripts/09_profile.R): runtime profiling
  (`Rprof()` / `profvis`); results in [example/profile_historial](example/profile_historial).
- [example/scripts/10_factorial_Design.R](example/scripts/10_factorial_Design.R):
  full factorial `A_max x omega2 x mutacion_sd` (3^3, 3 replicates each).
- [example/scripts/11_Definitive_Screening_Design.R](example/scripts/11_Definitive_Screening_Design.R):
  Definitive Screening Design over 9 spatio-temporal regimes with log-scaled parameters.
- [example/scripts/12_Replicated_Definitive_Screening_Design.R](example/scripts/12_Replicated_Definitive_Screening_Design.R):
  the same DSD with `N_REPLICAS` independent replicates per design point.

Outputs land in [example/data](example/data) (tables) and [example/plots](example/plots)
(plots and animations), keeping every example self-contained.

## Simulation engine at a glance

| Stage | File |
|---|---|
| Initial fields (nutrient, colony seeding, antibiotic template) | [R/01_campos.R](R/01_campos.R), [R/02_colonia.R](R/02_colonia.R) |
| PDE diffusion/reaction update and target antibiotic field | [R/03_pde_difusion.R](R/03_pde_difusion.R) |
| Cellular automaton growth, death, mutation | [R/04_ca_crecimiento.R](R/04_ca_crecimiento.R) |
| Orchestration, event ladder, main loop | [R/05_motor.R](R/05_motor.R) |
| Quantitative-genetics analytics | [R/06_genetica_cuantitativa.R](R/06_genetica_cuantitativa.R) |
| Post-processing: plots, animations, fractal dimension | [R/07_analisis.R](R/07_analisis.R) |

## Parameters and interpretation

A few parameters are especially important for interpreting the results:

- `A_max`: amplitude of the antibiotic shock after introduction.
- `omega2`: width of the stabilizing selection curve; smaller values imply stronger selection.
- `mutacion_sd`: standard deviation of the additive mutation step during division.
- `paso_introduccion`: time step at which the antibiotic shock is introduced.
- `guardar_cada`: interval between stored snapshots in `historial`; increase it to save
  memory on long runs (set it above `n_steps` to keep only the final state plus summaries).
- `Nc`: population threshold used to define quasi-extinction for persistence analysis.
- `multiresistance`, `factor_dosis`, `factor_exigencia`, `tope_espera_shock`:
  the escalating-event ladder described above.

These parameters control the balance between growth, stress, adaptation, and demographic decline.

## Outputs

A simulation returns everything needed for later analysis without re-running:

- `resumen`: per-step time series (population, proliferative front, births/deaths,
  nutrient/antibiotic means, trait mean/variance, active antibiotic level),
- `historial`: periodic spatial snapshots (state, fields, traits) for snapshots/animations,
- `estado_final`, `resist_final`: final colony layout and resistance landscape,
- `eventos_shock`: the antibiotic event ledger,
- `paso_extincion`: first step below `Nc` (quasi-extinction), `NA` if the population persisted,

plus the generated artifacts from scripts and the app: static plots, GIF/MP4 animations,
fractal-dimension estimates, and experiment tables.

## Development and testing

The test suite lives in [tests/testthat](tests/testthat) and combines unit tests
(engine, escalation ladder, plotting helpers) with browser end-to-end tests
(`shinytest2`). Run it with:

```r
devtools::test(".")
```

End-to-end tests need a Chromium-based browser; on machines without Chrome set
`CHROMOTE_CHROME` to the Edge executable, e.g.
`Sys.setenv(CHROMOTE_CHROME = "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe")`.

## Package publication note

The repository keeps the canonical code in [R](R), the app assets under
[inst/app](inst/app), compiled sources in [src](src) (optional acceleration), roxygen
documentation in [man](man), and example data/plots under [example](example). The legacy
[scripts](scripts) files are retained as compatibility shims so older workflows continue
to source the same functions without duplicating the implementation.
