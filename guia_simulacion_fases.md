## Fase 3: Recolección de Datos (El Barrido de Parámetros)
**Objetivo:** Automatizar el estudio estadístico. Una simulación aislada solo hace dibujitos; el objetivo de la Fase 3 es correr simulaciones masivas variando metódicamente los parámetros (especialmente el ruido $\eta$) para extraer y documentar fenómenos críticos, como la **Transición de Fase** (el punto exacto donde el caos se vuelve orden).

**Lógica Paso a Paso (Implementada en `fase3_recoleccion.py`):**

### Paso 3.1: Configuración del Enfriamiento Gradual (Evitando Histéresis)
En la materia activa, la historia del sistema importa. Si un enjambre ya está muy ordenado, requiere mucho más ruido para romperlo; si está caótico, requiere muchísimo menos ruido para que empiece a ordenarse. 
- En lugar de instanciar un simulador de cero para cada $\eta$, el script crea **un solo simulador**.
- Le asignamos el Ruido Máximo inicial (ej. $\eta=5.0$) asegurando el caos.
- Iterativamente, reducimos el ruido en pequeños escalones (como bajando la temperatura).

### Paso 3.2: La Termalización
- Cuando bajamos un escalón de ruido, el sistema experimenta un shock. Toma tiempo que el enjambre adopte la configuración natural para este nuevo nivel.
- Ejecutamos el método `step()` cientos de veces (`t_termalizacion`) pero **ignoramos absolutamente todo lo que pasa**. No medimos nada, solo esperamos a que el sistema alcance un estado macroscópicamente estable (estado estacionario).

### Paso 3.3: El Muestreo
- Tras la termalización, pasamos a la fase de medición durante un tiempo `t_muestreo`.
- En cada `step()`, extraemos la matriz completa de velocidades de las $N$ partículas.
- Invocamos la función matemática vectorizada para hallar **$\Phi$ (Parámetro de Orden Polar)**. Esto nos da un número puro entre $0$ (Desorden total) y $1$ (Flocking perfecto/Alineación total).
- Se guarda el valor de $\Phi$ de ese instante específico en la memoria RAM.

### Paso 3.4: Análisis y Exportación (CSVs)
- Finalizado el muestreo para el nivel de ruido actual, calculamos la Media ($\langle\Phi\rangle$) y la Varianza de todas las mediciones tomadas. La varianza nos alertará matemáticamente de las inestabilidades justo en el punto crítico de la transición.
- Guardamos esta fila en nuestro gran registro de Transición de Fase.
- *Opcionalmente:* Si `guardar_detalles=True`, tomamos las matrices de posición y velocidad absolutas medidas en el Paso 3.3 y las escribimos eficientemente a disco mediante la librería `Pandas`.
- Reiniciamos el ciclo volviendo al Paso 3.1 pero con un ruido aún menor.

---

## La Filosofía del "Adapter" (¿Cómo interactúa la Fase 3 con la 1 y 2?)

El módulo recolector de datos se ha programado para ser **completamente ciego y agnóstico** a la física subyacente. No le interesa si estás usando Vicsek (Fase 1) o EDOs continuos (Fase 2).

Él funciona bajo un contrato estricto: te pide que le pases un objeto que cumpla cuatro cosas simples:
1. Una propiedad modificable llamada `.eta`.
2. Un método que avance el tiempo llamado `.step()`.
3. Un método que devuelva una matriz NumPy $(N \times 2)$ con las velocidades de todos llamado `.get_velocities()`.
4. Un método análogo llamado `.get_positions()`.

Para conectar esto con tu código actual de Vicsek (donde tú manejas ángulos y no vectores de velocidad), implementamos en el archivo un **SimuladorAdapter**. Este componente actúa como traductor: cuando la Fase 3 le pide `get_velocities()`, el adaptador por detrás lee los ángulos $\theta$ de tu Fase 1, calcula el seno y coseno multiplicados por $v_0$, y le entrega la matriz perfecta de retorno a la Fase 3. 

Cuando construyas la Fase 2, construirás un nuevo adaptador, que en su lugar de leer ángulos, extraerá las velocidades integradas numéricamente y las pasará, permitiendo reutilizar toda tu automatización de exportación de datos intacta.
