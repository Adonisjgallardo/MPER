# MPER
Adaptive Reorganization and Ecological Perturbation Model

MPER is a spatial, agent-based model for studying how microbial populations respond to abrupt environmental change. The repository is organized to follow package-style conventions so that the simulation engine and the Shiny app can be published as a single R package.

## Package layout

The project now follows a package-oriented structure:

- [R](R): package source code for the simulation engine, Shiny app modules, and helpers.
- [inst/app](inst/app): Shiny app entrypoint and static assets.
- [scripts](scripts): compatibility wrappers that source the canonical code in [R](R).
- [example](example): reproducible example workflows and ordered outputs.

### Core simulation engine

The canonical simulator lives in [R](R). The main entry point is [R/05_motor.R](R/05_motor.R), which runs the hybrid simulation through the following stages:

1. Create the nutrient field and antibiotic template with [R/01_campos.R](R/01_campos.R) and [R/02_colonia.R](R/02_colonia.R).
2. Update the nutrient and antibiotic PDE fields in [R/03_pde_difusion.R](R/03_pde_difusion.R).
3. Apply the cellular automaton dynamics in [R/04_ca_crecimiento.R](R/04_ca_crecimiento.R).
4. Record summaries and snapshots in [R/05_motor.R](R/05_motor.R).
5. Compute post-processing metrics and plots in [R/07_analisis.R](R/07_analisis.R).

### Shiny app

The app entry point is [inst/app/app.R](inst/app/app.R). It loads the package files from [R](R) and launches the UI and server defined in [R/app_ui.R](R/app_ui.R) and [R/app_server.R](R/app_server.R).

To launch the app from a source checkout:

```r
shiny::runApp("inst/app")
```

Alternatively, after loading the package in an R session:

```r
run_mper()
```

### Example workflows

The example scripts live under [example/scripts](example/scripts) and are intended as reproducible entry points for users:

- [example/scripts/08_simulate.R](example/scripts/08_simulate.R): runs a single example simulation and writes outputs to [example/data/simulation](example/data/simulation) and [example/plots/simulation](example/plots/simulation).
- [example/scripts/09_profile.R](example/scripts/09_profile.R): profiles the simulation runtime and writes results to [example/profile_historial](example/profile_historial).
- [example/scripts/10_factorial_Design.R](example/scripts/10_factorial_Design.R): runs a small factorial experiment over key parameters.

The outputs are organized as follows:

- [example/data](example/data): simulation outputs and experiment results.
- [example/plots](example/plots): plots and animations.
- [example/profile_historial](example/profile_historial): saved profiling snapshots.

## Parameters and interpretation

A few parameters are especially important for interpreting the results:

- `A_max`: amplitude of the antibiotic shock after introduction.
- `omega2`: width of the stabilizing selection curve; smaller values imply stronger selection.
- `mutacion_sd`: standard deviation of the additive mutation step during division.
- `paso_introduccion`: time step at which the antibiotic shock is introduced.
- `Nc`: population threshold used to define quasi-extinction for persistence analysis.
- `multiresistance`: number of ADDITIONAL antibiotic events (each one harsher than the previous).
  Event `k` scales the dose geometrically (`A_max * factor_dosis^k`) and hardens selection
  (`omega2` shrinks, `theta1` is pushed toward 1, `mort_estres` grows, all governed by
  `factor_exigencia`). Events fire when the population recovers its previous peak OR when
  `tope_espera_shock` steps elapse since the last event. The simulation result includes an
  `eventos_shock` table (step, level, effective parameters, trigger type) and a per-step
  `nivel_shock` column in `resumen`.

These parameters control the balance between growth, stress, adaptation, and demographic decline.

## Outputs

The scripts generate several useful artifacts:

- summaries of population size and trait statistics,
- snapshots of the colony state across time,
- fractal-dimension estimates,
- static plots and animations of the spatial dynamics,
- experiment tables for factorial exploration.

## Package publication note

The repository now keeps the canonical code in [R](R), the app assets under [inst/app](inst/app), and the example data/plots under [example](example). The legacy [scripts](scripts) files are retained as compatibility shims so older workflows continue to source the same functions without duplicating the implementation.
