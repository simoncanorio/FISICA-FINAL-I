# Materia Activa y Autoorganización Colectiva (Flocking)

Repositorio del **Trabajo Final de Física I**.

Este proyecto trata sobre la dinámica de sistemas de materia activa, centrándose en el fenómeno de autoorganización colectiva (*flocking*). A través de simulaciones y análisis, exploramos cómo reglas locales simples dan lugar a comportamientos emergentes complejos y cómo estos sistemas transitan del desorden al orden bajo la influencia del ruido.

## En en repositorio se encuentran subidas las simulaciones y los graficos obtenidos en el final del analisis, como también sus respectivos códigos y el código que se realizó para la recolección de los datos. Luego tambien se encuentran subidos dos PDF, uno que contiene el INFORME FINAL y otro que contiene un link para visualizar el VIDEO EXPLICATIVO DE LOS REUSLTADOS, al pesar demasiado tuvimos que entregarlo por este medio. Para poder ver todo, descargar la carpeta del repositorio en la computadora y acceder a los PDF desde ahí. 

## Enfoques de Modelado
El trabajo compara cuatro enfoques distintos para simular el comportamiento de enjambre:

1. **Enfoque discreto (Modelo de Vicsek):** Basado en la alineación angular instantánea. Implementamos una variante con **interacción local** (vecindario de radio $R$) y una variante **All-to-All** (campo medio global).
2. **Enfoque continuo:** Implementación de **Ecuaciones Diferenciales Ordinarias (EDOs)** que incorporan masa, inercia angular, fuerzas de atracción/repulsión espacial y torques de alineación. Se comparan variantes con **interacción local** y **campo medio (sin radio)**.

## Estructura del Proyecto (Fases)

El desarrollo se organizó en cuatro fases metodológicas:

* **Fase 1 (Modelo Vicsek Local):** Construcción del motor discreto utilizando `cKDTree` para optimizar la búsqueda de vecinos ($O(N \log N)$).
* **Fase 2 (Modelo Continuo con EDOs):** Desarrollo del simulador de física continua mediante el método de integración de Euler-Maruyama.
* **Fase 3 (Recolección de Datos):** Barrido sistemático de parámetros ($\eta$ ruido). Generación de datasets para modelos locales y de campo medio.
* **Fase 4 (Análisis Visual y Matemático):** Procesamiento de datos para calcular:
    * **Parámetro de Orden Polar ($\Phi$):** Medida de la sincronización colectiva.
    * **Susceptibilidad ($\sigma(\Phi)$):** Identificación precisa del punto crítico de transición ($\eta_c$).
    * **Correlación Espacial $C(r)$:** Análisis del decaimiento del orden en el espacio.

## Contenido del Repositorio

* `/fase1_vicsek.py` y `/fase1_vicsek_all2all.py`: Motores de simulación discreta.
* `/fase2_continuo.py` y `/continuo_sin_radio.py`: Motores de simulación continua (EDOs).
* `/fase3_recoleccion.py` / `/fase3_recoleccion_v2.py`: Scripts para la ejecución masiva de experimentos y generación de archivos `.csv`.
* `/fase4_analisis.py`: Script de análisis estadístico y generación de gráficos científicos (figuras de transición, susceptibilidad y correlación).
* `/animaciones/`: GIFs e imágenes representativas de las simulaciones.
* `Informe_TpFinal.pdf`: Documento completo con el desarrollo, metodología y conclusiones.

## Hallazgos Principales

1. **Robustez del Campo Medio:** Los modelos de acoplamiento global demostraron una resistencia superior al ruido, manteniendo el orden ($\Phi \approx 1.0$) hasta valores críticos elevados ($\eta_c \approx 1.9 - 2.0$), comparado con los modelos locales ($\eta_c \approx 1.8$).
2. **Anomalía en el Continuo Local:** La inclusión de fuerzas físicas realistas (repulsión volumétrica) fragmenta la cohesión global, impidiendo que el enjambre alcance un orden pleno y generando un comportamiento altamente inestable ($\eta_c \approx 0.46$).
3. **Límites de Tamaño Finito:** El análisis de correlación espacial confirmó que las dimensiones finitas del dominio (toroide) imponen límites a la propagación del orden, observándose un decaimiento exponencial en lugar de leyes de potencias puras.

## Alumnos
* **Narda Rosa**
* **Simon Canorio**
* **Felipe Di Risio Ize**
