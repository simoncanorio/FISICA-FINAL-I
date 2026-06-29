# Fase 3 (v2): Recolección de Datos — Modelos de Acoplamiento Global

Esta fase replica el protocolo experimental de la Fase 3 original, pero aplicado a los dos modelos de **campo medio** (acoplamiento global): el Vicsek All-to-All y el Continuo Sin Radio. La diferencia clave respecto a la versión anterior es que **no existe radio de interacción local**: cada partícula se alinea con el promedio de *todas* las demás partículas del sistema simultáneamente.

---

## 1. El Parámetro de Orden Polar ($\Phi$)

El **Parámetro de Orden Polar** es la magnitud central de este experimento. Condensa el grado de acuerdo colectivo del enjambre en un único número:

$$\Phi = \frac{1}{N} \left| \sum_{i=1}^{N} \hat{v}_i \right|$$

donde $\hat{v}_i$ es el **vector unitario** en la dirección de movimiento de la partícula $i$.

| Valor de $\Phi$ | Interpretación física |
| :---: | :--- |
| $\Phi \approx 1.0$ | Alineación perfecta: todas las partículas van en la misma dirección |
| $\Phi \approx 0.0$ | Movimiento completamente aleatorio: los vectores se cancelan entre sí |

### Implementación diferenciada por modelo

| Modelo | Fórmula usada | Razón |
| :--- | :--- | :--- |
| **Vicsek All-to-All** | $\frac{1}{N}\|\sum_i (\cos\theta_i, \sin\theta_i)\|$ | Rapidez constante → los ángulos bastan |
| **Continuo Sin Radio** | $\frac{1}{N}\|\sum_i \mathbf{v}_i / \|\mathbf{v}_i\|\|$ | Rapidez variable → se normaliza cada vector antes de sumar |

---

## 2. Los Modelos de Acoplamiento Global

### ¿Qué los diferencia de los modelos originales?

Los modelos de la Fase 3 original (**Vicsek con radio R** y **Continuo con vecindario local**) usan interacciones *locales*: cada partícula solo se ve influenciada por sus vecinas dentro de un radio $R$. Los modelos de esta fase usan **campo medio global**:

```
Modelo local (Fase 3 original):
  Partícula i  →  promedio de vecinos con |r_j - r_i| < R

Modelo campo medio (esta fase):
  Partícula i  →  promedio de TODAS las N partículas
```

Esta diferencia tiene consecuencias físicas profundas:

- El modelo de campo medio **converge al orden más rápido** para bajos niveles de ruido, ya que la información de alineación se propaga instantáneamente a todo el sistema.
- La transición de fase puede estar en un **eta crítico distinto** al del modelo local.
- En el límite $N \to \infty$, el campo medio reproduce la ecuación de Curie-Weiss del magnetismo.

---

## 3. El Barrido de Parámetros

El experimento recorre una grilla de valores del nivel de ruido $\eta$ **de mayor a menor** (de caos a orden):

```
eta_alto (caos)  →  eta_medio  →  eta_bajo (orden)
```

### Protocolo para cada valor de $\eta$

```
┌─────────────────────────────────────────────────────────┐
│  1. Inicializar posiciones y ángulos al AZAR            │
│  2. Correr N_termal pasos  → Estado estacionario        │
│     (los datos de esta fase se DESCARTAN)               │
│  3. Correr N_medicion pasos → Serie temporal de Φ(t)   │
│     (estos datos se PROMEDIAN y GUARDAN)                │
│  4. Exportar snapshot del estado final                  │
│  5. Avanzar al siguiente valor de eta                   │
└─────────────────────────────────────────────────────────┘
```

> [!IMPORTANT]
> La etapa de **termalización** es fundamental. Si empezamos a medir antes de que el sistema alcance su estado estacionario, obtenemos valores de $\Phi$ contaminados por el transitorio inicial.

### Parámetros del experimento

| Parámetro | Vicsek All-to-All | Continuo Sin Radio | Explicación |
| :--- | :---: | :---: | :--- |
| `ETA_MAX` | 2.0 | 2.0 | Ruido máximo explorado |
| `ETA_MIN` | 0.05 | 0.05 | Ruido mínimo explorado |
| `N_PUNTOS` | 20 | 20 | Cantidad de valores de $\eta$ |
| `N` (partículas) | 150 | 150 | Número de partículas |
| `L` (caja) | 20.0 | 20.0 | Tamaño del espacio |
| `v` / `f_0` | 0.5 | 1.0 | Velocidad / autopropulsión |
| `k_align` | — | 5.0 | Fuerza del torque de alineación |
| `pasos_termal` | 500 | 3000 | Pasos de termalización |
| `pasos_medicion` | 300 | 1500 | Pasos de medición promediados |
| `dt` | 1.0 | 0.02 | Paso de tiempo |

> [!NOTE]
> El modelo continuo necesita **más pasos de termalización** porque su $\Delta t$ es mucho más pequeño (0.02 vs 1.0 en Vicsek). En tiempo físico simulado, ambos equilibran en un rango comparable.

---

## 4. Archivos CSV Generados

El script genera **cuatro archivos separados** con distintos niveles de detalle:

### `resultados_vicsek_all2all.csv` y `resultados_continuo_sin_radio.csv`

**Un registro por cada valor de $\eta$**. Este es el archivo principal de análisis.

| Columna | Descripción |
| :--- | :--- |
| `eta` | Nivel de ruido del experimento |
| `phi_promedio` | $\langle\Phi\rangle$ promediado sobre todos los pasos de medición |
| `phi_std` | Desviación estándar de $\Phi$ (fluctuaciones en el estado estacionario) |
| `phi_min` | Valor mínimo de $\Phi$ observado durante la medición |
| `phi_max` | Valor máximo de $\Phi$ observado durante la medición |
| `pasos_termal` | Pasos de termalización usados |
| `pasos_medicion` | Pasos de medición promediados |
| `N`, `L`, `v`/`dt` | Parámetros de la simulación para reproducibilidad |
| `k_align` | *(Solo continuo)* Fuerza de alineación global |

### `snapshots_vicsek_all2all.csv` y `snapshots_continuo_sin_radio.csv`

**Un registro por partícula por cada valor de $\eta$**. Permite reconstruir el estado microscópico del sistema al estabilizarse.

| Columna | Descripción |
| :--- | :--- |
| `eta` | Nivel de ruido del experimento |
| `paso` | Número de paso total (termal + medición) |
| `particula_id` | Índice de la partícula (0 a N-1) |
| `x`, `y` | Posición en el espacio |
| `theta` | Ángulo de orientación (rad) |
| `vx`, `vy` | *(Solo continuo)* Componentes de velocidad |
| `omega` | *(Solo continuo)* Velocidad angular (rad/s) |

---

## 5. Cómo Ejecutar el Barrido

Asegurarse de tener `fase1_vicsek_all2all.py` y `continuo_sin_radio.py` en el mismo directorio.

```bash
# Ejecutar el barrido completo (ambos modelos de campo medio)
python fase3_recoleccion_v2.py
```

La salida en terminal mostrará el progreso en tiempo real:

```
╔══════════════════════════════════════════════════════════╗
║    FASE 3 (v2) — RECOLECCIÓN MASIVA DE DATOS            ║
║    Modelos de Acoplamiento Global (Campo Medio)          ║
╚══════════════════════════════════════════════════════════╝

============================================================
  BARRIDO VICSEK ALL-TO-ALL | N=150 | L=20.0 | v=0.5
  Modo            : Campo Medio Global (sin radio de corte)
============================================================
  [#####                         ] eta=2.000  <Phi>=0.0398 ± 0.0287  (0.1s)
  [##########                    ] eta=1.895  <Phi>=0.0512 ± 0.0341  (0.1s)
  ...
  [##############################] eta=0.050  <Phi>=0.9954 ± 0.0012  (0.1s)

============================================================
  BARRIDO CONTINUO SIN RADIO | N=150 | L=20.0 | dt=0.02
  Modo            : Acoplamiento Global (sin radio de vecindario)
============================================================
  [#####                         ] eta=2.000  <Phi>=0.0441 ± 0.0319  (2.1s)
  ...
```

---

## 6. Comparación con los Modelos de la Fase 3 Original

| Característica | Vicsek R local | Vicsek All-to-All | Continuo local | Continuo sin radio |
| :--- | :---: | :---: | :---: | :---: |
| Radio de interacción | ✓ R=2.0 | ✗ (global) | ✓ local | ✗ (global) |
| Tipo de acoplamiento | Local | Campo medio | Local | Campo medio |
| Velocidad de convergencia | Media | Alta | Media | Alta |
| Dependencia de densidad | Sí | No | Sí | No |
| Archivo de modelo | `fase1_vicsek.py` | `fase1_vicsek_all2all.py` | `fase2_continuo.py` | `continuo_sin_radio.py` |
| CSV resultado | `resultados_vicsek.csv` | `resultados_vicsek_all2all.csv` | `resultados_continuo.csv` | `resultados_continuo_sin_radio.csv` |

> [!TIP]
> Comparar las curvas $\langle\Phi\rangle$ vs $\eta$ entre los modelos locales y de campo medio es muy revelador: el campo medio tiende a mostrar una **transición más abrupta** y un $\eta_c$ mayor, ya que la alineación global es más eficiente.

---

## 7. Próximos Pasos (Fase 4)

Los archivos CSV generados en esta fase son la materia prima para el análisis gráfico:

- **Curva $\langle\Phi\rangle$ vs $\eta$**: Visualizar la transición de fase orden-desorden para los modelos de campo medio.
- **Comparación entre los 4 modelos**: Superponer las curvas de los cuatro modelos (Vicsek local, Vicsek all2all, Continuo local, Continuo sin radio) para comparar $\eta_c$ y la forma de la transición.
- **Análisis de fluctuaciones**: Graficar $\sigma(\Phi)$ vs $\eta$ para identificar el punto crítico como pico de susceptibilidad.
- **Visualización espacial**: Usar los snapshots para dibujar configuraciones del enjambre en distintos regímenes de $\eta$.
