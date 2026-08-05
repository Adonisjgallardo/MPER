# MPER
Adaptative Reorganization And Ecological Perturbation Model

## Intro OLD

MPER
Modelo de Perturbación Ecológica y Reorganización Adaptativa
Idea central

Una comunidad microbiana muestra una **diversidad fenotípica preexistente**. La heterogeneidad ambiental mantiene esa diversidad. Una *perturbación ecológica* modifica el paisaje adaptativo y provoca una reorganización de las **frecuencias** fenotípicas.

Comunidad -> Diversidad inicial -> Perturacion ambietal -> Diversidad final

La resistencia antimicrobiana sería un caso particular cuando una de las dimensiones del fenotipo representa tolerancia o resistencia.

Comunidad -> Diversidad inicial -> Perturacion Antibiotico -> Diversidad final

## Estructura del proyecto

| Script                  | Función                         | Entrada                        | Salida                                                   |   |
| ----------------------- | ------------------------------- | ------------------------------ | -------------------------------------------------------- | - |
| **01_paisaje.R**        | Construye el paisaje adaptativo | Ninguna (solo parámetros)      | `paisaje_inicial.rds`                                    |   |
| **02_bacterias.R**      | Genera la población inicial     | `paisaje_inicial.rds`          | `bacterias_iniciales.rds`                                |   |
| **03_dinamica.R**       | Evolución antes del disturbio   | Paisaje + bacterias            | `historial_evolucion.rds`                                |   |
| **04_disturbio.R**      | Modifica el ambiente            | Historial + paisaje            | `paisaje_perturbado.rds` y `bacterias_pre_disturbio.rds` |   |
| **05_reorganizacion.R** | Evolución después del disturbio | Paisaje perturbado + bacterias | `historial_post_disturbio.rds`                           |   |

### Parámetros del paisaje *(01_paisaje.R)*

El paisaje es una cuadrícula de 100 × 100 celdas donde cada ua posee cuatro variables ambientales que representan dimensiones ecológicas,

$E_1$: tolerancia antimicrobiana
$E_2$: crecimiento
$E_3$: persistencia
$E_4$: dispersión

### Población inicial *(02_bacterias.R)*

Cada bacteria posee propiedades asociadas a su estado: 

$P = (F_1, F_2, F_3, F_4)$

Dode:

$F_1$: identificador
$F_2$: posición (x,y)
$F_3$: cuatro rasgos fenotípicos
$F_4$: edad y generación (conceptualmente)


La población inicial contiene aproximadamente 10 000 individuos y no está formada por clones.

### Parámetros evolutivos *(03_dinamica.R)*

Aquí ocurre la evolución.

Cada ciclo de la simulacion fija 100 generacioes antes y despues del disturbio. Las mutaciones se suponen de distribucion normal y se acerca el tamaño poblacional a 10 000 individuos por generacion. 

En cada generación:

a) Cada bacteria lee el ambiente donde está ubicada.
) Se calcula la distancia $d$ entre su fenotipo y el ambiente.
c) Esa distancia se convierte en un valor de fitness $\rho$ mediante

$$\rho = e^{−d}$$

d) Sobreviven probabilísticamente según ese fitness.
e) Los sobrevivientes producen descendencia.
f) La descendencia recibe pequeñas mutaciones.

Esto asegura que cada individuo posea movimiento, capacidad para leer el ambiente y responder al cálcular la distancia y el fitness. Asi tamie se asegura que exista la mortalidad, la reproducción y la mutabilidad.

### Disturbaciones *(04_disturbio.R)*

Solo modificamos el paisaje y permitimos que los individuos reaccionen al cambio. La comunidad responde.

### Graficos *(05_visualizacion.R)*

Archivos de datos

El modelo genera los siguientes archivos:

    *paisaje_inicial.rds*
    *bacterias_iniciales.rds*
    *historial_evolucion.rds*
    *paisaje_perturbado.rds*
    *bacterias_pre_disturbio.rds*
    *historial_post_disturbio.rds*

Estos archivos de datos se usan para producir graficos de: 

    *mapa del paisaje adaptativo*
    *distribución espacial de bacterias*
    *espacio fenotípico* ($F_1$ vs $F_2$)
    *comparación antes y después del disturbio*
    *evolución temporal* (propuesta en el documento)

Estos graficos se usan para revisar la hipotesis del mapa del paisaje, la distribución bacteriana, la evolución temporal, mapas de densidad, PCA y trayectorias.

## Resultados biológicos esperados

Si la hipótesis es correcta, después del disturbio se observará que:

    algunos fenotipos aumentan su frecuencia;
    otros disminuyen o desaparecen;
    la composición de la comunidad cambia;
    la diversidad observada proviene de la reorganización de variantes preexistentes, no de la creación de nuevas variantes por el antibiótico.

### Observaciones técnicas

El modelo es una *simulación basada en agentes* (Agent-Based Model, ABM) con un paisaje adaptativo espacial. Conceptualmente es sólido para explorar la hipótesis planteada, aunque en esta versión presenta varias simplificaciones:

Las bacterias no se desplazan entre celdas, aunque el documento menciona movimiento como parte de la dinámica.

No hay competencia por recursos ni capacidad de carga.

La reproducción es asexual y proporcional al fitness.

El fitness depende únicamente de la distancia euclidiana entre el fenotipo y el ambiente.

El disturbio modifica únicamente la dimensión $E_1$, representando un cambio ambiental localizado. 

Se implementa una primera versión de un modelo evolutivo espacial diseñado para evaluar **si la resistencia antimicrobiana puede interpretarse como una reorganización adaptativa de diversidad fenotípica preexistente**, más que como la aparición de nuevos rasgos inducidos por el ambiente.

## Intro nEW
Execution pipeline (from 08_simulate.R)
Initialisation (01_campos.R, 02_colonia.R)

Create nutrient field N (uniform + optional noise).

Create spatial antibiotic template (uniform/linear/radial).

Seed a colony at the centre with initial genetic value g_inicial and environmental deviation e ~ N(0, σₑ²).

Main time loop (inside simular_hca() in 05_motor.R)

PDE updates (03_pde_difusion.R)

Nutrient: diffusion + consumption by bacteria (Michaelis‑Menten).

Antibiotic: diffusion + degradation + relaxation toward a time‑varying target (zero before shock, then a spatial gradient with optional temporal modulation).

Environmental resampling (resortear_ambiente) – redraw e for every living cell each step (quantitative‑genetics assumption).

Cellular automaton step (crecer_colonia in 04_ca_crecimiento.R)

For each occupied cell (in random order):

Compute local fitness w = exp(−(θ − z)²/(2ω²)).
Determine death probability: p_mort = mort_base + mort_estres * A/(A+K_mort) * (1−w). If death occurs, the site is cleared.
If alive, check if it can divide (requires local nutrient above threshold and space in von Neumann neighbourhood).
Division probability: p_div = p_max * (N/(N+N_half)) * w. If successful, offspring inherits g with mutation (Gaussian) and gets a fresh e.
Data recording – every guardar_cada steps a full snapshot (occupancy, N, A, phenotype) is stored; also a resumen data frame with population statistics.

Post‑processing (07_analisis.R, 08_simulate.R)

Compute fractal dimension of the colony (box‑counting).

Compare simulated trajectories with analytical predictions (quantitative genetics, persistence times).

Generate static plots and animations (ggplot2 + magick/av).

