# MPER
Adaptative Reorganization and Ecological Perturbation Model

MPER is a spatial, agent-based model for studying how microbial populations respond to abrupt environmental change. The main idea is to test whether antimicrobial resistance can emerge as a reorganization of pre-existing phenotypic variation rather than as the creation of entirely new traits.

## What the model does

The core simulation combines:

- a nutrient field that diffuses and is consumed by the colony,
- an antibiotic field that is initially absent and is introduced as an abrupt shock,
- a cellular automaton for birth and death events,
- quantitative-genetic dynamics in which each cell carries a heritable value $g$, a non-heritable environmental component $e$, and a realized phenotype $z = g + e$.

The colony starts from a small inoculum and evolves under stabilizing selection around a local optimum that shifts when antibiotic stress increases. The model tracks whether the population persists, collapses, or recovers after the shock.

## Main workflow

The project currently has two layers:

1. The simulation engine in [scripts](scripts)
2. Example workflows in [example/scripts](example/scripts)

### 1. Core simulation engine

The main entry point is [scripts/05_motor.R](scripts/05_motor.R), which runs the full hybrid simulation through the following steps:

1. Create the nutrient field and antibiotic template with [scripts/01_campos.R](scripts/01_campos.R) and [scripts/02_colonia.R](scripts/02_colonia.R).
2. Update the nutrient and antibiotic PDE fields in [scripts/03_pde_difusion.R](scripts/03_pde_difusion.R).
3. Apply the cellular automaton dynamics in [scripts/04_ca_crecimiento.R](scripts/04_ca_crecimiento.R).
4. Record summaries and snapshots for each time step in [scripts/05_motor.R](scripts/05_motor.R).
5. Compute post-processing metrics and plots in [scripts/07_analisis.R](scripts/07_analisis.R).

The simulation object returned by the main function includes:

- a time-series summary in `resumen`,
- a history of snapshots in `historial`,
- the final state of the colony and fields,
- the extinction step, if one occurs.

### 2. Example workflow

The example scripts live in [example/scripts](example/scripts) and are intended to be the reproducible entry points for users:

- [example/scripts/08_simulate.R](example/scripts/08_simulate.R): runs a single example simulation and exports outputs to [example/data/simulation](example/data/simulation) and [example/plots/simulation](example/plots/simulation).
- [example/scripts/09_profile.R](example/scripts/09_profile.R): profiles the simulation runtime and writes profiling output to [example/profile_historial](example/profile_historial).
- [example/scripts/10_factorial_Design.R](example/scripts/10_factorial_Design.R): runs a small factorial experiment over key parameters such as stress intensity, selection width, and mutation rate.

The example outputs are organized as follows:

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

These parameters control the balance between growth, stress, adaptation, and demographic decline.

## Outputs

The scripts generate several useful artifacts:

- summaries of population size and trait statistics,
- snapshots of the colony state across time,
- fractal-dimension estimates,
- static plots and animations of the spatial dynamics,
- experiment tables for factorial exploration.

## Suggested next step

As the project grows, a natural next step would be to move the shared R functions into a package-style structure under [R](R), while keeping [example](example) as a place for reproducible demonstrations and experiments.
