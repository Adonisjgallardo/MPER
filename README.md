# MPER
Adaptative Reorganization And Ecological Perturbation Model

## Intro

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

