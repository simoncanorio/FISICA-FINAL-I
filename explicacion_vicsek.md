# Construyendo el Enjambre: Modelo de Vicsek (Fase 1)

Aquí tienes la explicación paso a paso de lo que hicimos y cómo funciona, explicado de la forma más sencilla posible.

## ¿Qué hicimos?

Creamos la base física para que un grupo de "partículas" (como si fueran pájaros, peces o bacterias) comiencen a moverse juntos como un enjambre, basándose en el clásico **Modelo de Vicsek**.

### La Simulación en Acción

![Animación del Modelo de Vicsek](./vicsek_anim.gif)

> [!TIP]
> Observa el GIF: al principio las partículas apuntan a cualquier lado (son de muchos colores mezclados). Sin embargo, con el paso del tiempo comienzan a influenciarse, se alinean y terminan moviéndose en grandes "bandadas" con la misma dirección y color. ¡Es la emergencia del orden a partir del caos local!

---

## El Paso a Paso (De Forma Simple)

Para lograr esta simulación, dividí el problema en tres partes fundamentales dentro del código:

### 1. El Tablero Sin Fin (Límites Periódicos)
Imagínate el mapa de la simulación como el mundo del clásico juego *Pac-Man*: si una partícula sale volando por el lado derecho de la pantalla, no choca contra un muro, sino que aparece mágicamente por el lado izquierdo manteniendo su velocidad.

- **¿Cómo se hizo?** Cuando calculo la nueva posición matemática de una partícula sumando su velocidad, uso una operación llamada "módulo" (`%`). Esto recorta automáticamente cualquier número que se pase del límite máximo del mapa y lo reinicia desde el otro lado.

### 2. Encontrar a los Vecinos Rápidamente (Optimización)
Para que el enjambre funcione, cada partícula solo debe seguir a los que tiene en un radio cercano (sus vecinos). Si tienes miles de partículas y cada una tiene que medir la distancia contra *todas las demás*, el ordenador colapsa.

- **¿Cómo se hizo?** Usé una herramienta matemática muy inteligente llamada `cKDTree` (Árbol KD). Imagínalo como un organizador espacial por zonas. En lugar de medir distancias una a una, la partícula pregunta a qué "zona" pertenece y el árbol le entrega inmediatamente a todos los que están cerca sin tener que buscar en todo el mapa. ¡Y además este árbol sabe cómo calcular distancias considerando el efecto "Pac-Man" de los bordes!

### 3. La Regla de "Seguir a los Demás" (Alineación)
En cada momento, cada partícula observa la dirección en la que van sus vecinos. Luego, ajusta su propia dirección hacia el *promedio* de ellos, pero añadiendo un pequeño temblor o "error" al azar. Ese ruido es fundamental para que se vea natural y fluido.

- **¿Cómo se hizo?** Tomé el ángulo de movimiento de los vecinos. Como sumar ángulos directamente es un problema matemático (promediar el grado 359° con el 1° te daría 180° que es ir para atrás, cuando en realidad están casi yendo hacia el mismo lado: 0°), convertí sus fuerzas en vectores (componentes horizontales y verticales usando Seno y Coseno). Sume esas fuerzas juntas para obtener la dirección real de la manada, y le sumé un empujón aleatorio para representar el "ruido".

---

## Los Archivos que Creados en tu Proyecto

Para mantener el código limpio, dividimos todo en dos archivos:

1. **[fase1_vicsek.py](file:///c:/Users/user/Downloads/fisica_final/fase1_vicsek.py):** Es el motor del coche. Contiene las fórmulas, matemáticas y reglas físicas descritas arriba. Está hecho para ser ultra-rápido y calcular los números puros.
2. **[animacion_vicsek.py](file:///c:/Users/user/Downloads/fisica_final/animacion_vicsek.py):** Es la carrocería. Toma los números generados por el motor y utiliza la librería **Matplotlib** para dibujar las flechitas, darles colores brillantes y tomar las "fotos" que conforman la animación que ves en este documento.

---

## Explicación Técnica Avanzada

Para profundizar en el diseño arquitectónico y de optimización (útil para reportes técnicos, documentación formal o papers), aquí tienes el detalle de cómo se resolvieron los desafíos matemáticos:

### 1. Topología Toroidal en la Cinemática
La condición de borde periódico $L \times L$ requiere mapear coordenadas que exceden el dominio de vuelta a su interior sin introducir saltos numéricos que corrompan el sistema.
- **Implementación Vectorizada:** La actualización posicional $x_{t+1} = x_t + v\Delta t \cos(\theta_t)$ se mantiene estrictamente dentro del intervalo $[0, L)$ aplicando el operador módulo vectorial: `self.positions %= self.L`. Al ser Python un lenguaje que preserva el signo del divisor en la operación módulo, las coordenadas negativas (partículas que escapan por $x < 0$) son mapeadas matemáticamente cerca de $L$ sin requerir sentencias condicionales (`if/else`) que romperían el paralelismo de los arrays.

### 2. Búsqueda Espacial Sub-cuadrática con `cKDTree`
Evaluar distancias euclidianas por pares (pairwise) posee una complejidad algorítmica de $O(N^2)$, lo cual es un cuello de botella inaceptable para simular grandes números de partículas interactivamente.
- **Optimización a $O(N \log N)$:** Implementamos árboles k-dimensionales en su variante compilada nativamente en C (`scipy.spatial.cKDTree`), reduciendo el tiempo de búsqueda drásticamente.
- **Métrica Toroidal Nativa:** El KDTree necesita saber que el espacio es curvo/periódico. Al instanciar `cKDTree(boxsize=self.L)`, la estructura adopta intrínsecamente la topología de toroide para calcular las distancias euclidianas. Internamente resuelve $\Delta x_{ij} = \min(|x_i - x_j|, L - |x_i - x_j|)$ en rutinas de C optimizadas, evitando que tengamos que replicar artificialmente las partículas en los bordes (técnica de "ghost cells" o cajas fantasma).

### 3. Álgebra de Giros Vectorizada
La ecuación de actualización del modelo de Vicsek dicta alinear los ángulos con el promedio del vecindario: $\theta_i(t+1) = \langle \theta(t) \rangle_r + \Delta\theta$
- **Elusión de Discontinuidades Topológicas:** Promediar números angulares de manera aritmética sufre discontinuidades insalvables en la frontera de $-\pi$ y $\pi$ (ej. el promedio numérico de $179^\circ$ y $-179^\circ$ da $0^\circ$, que apuntaría en sentido totalmente opuesto). La solución analítica robusta es mapear los ángulos de vuelta al espacio vectorial $\mathbb{R}^2$.
- **Uso de `arctan2`:** Usando NumPy, descomponemos las velocidades direccionales en `np.cos()` y `np.sin()` para todos los vecinos. Tras sumar estos componentes espaciales para encontrar el vector director promedio, la fase polar resultante se extrae usando `np.arctan2(sum_sin, sum_cos)`. Esto evita cualquier ambigüedad de cuadrante trigonométrico, ignora la necesidad de normalizar la magnitud del vector resultante, y ejecuta todas las iteraciones en operaciones vectoriales optimizadas subyacentes.
