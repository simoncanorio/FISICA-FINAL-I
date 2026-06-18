# Fase 3: Recolección de Datos — Barrido de Parámetros

En esta fase dejamos de construir el simulador y pasamos a usarlo como un **instrumento científico**. Ejecutamos ambos modelos de forma masiva y automática, barriendo el espacio de parámetros para mapear cómo el sistema transita entre el orden y el caos.

---

## 1. El Parámetro de Orden Polar ($\Phi$)

El **Parámetro de Orden Polar** es la magnitud central de este experimento. Es un número único que condensa el grado de "acuerdo" colectivo de todo el enjambre en un instante dado:

$$\Phi = \frac{1}{N} \left| \sum_{i=1}^{N} \hat{v}_i \right|$$

donde $\hat{v}_i$ es el **vector unitario** en la dirección de movimiento de la partícula $i$.

| Valor de $\Phi$ | Interpretación física |
| :---: | :--- |
| $\Phi \approx 1.0$ | Alineación perfecta: todas las partículas van en la misma dirección |
| $\Phi \approx 0.0$ | Movimiento completamente aleatorio: los vectores se cancelan entre sí |

### Implementación diferenciada por modelo

| Modelo | Fórmula usada | Razón |
| :--- | :--- | :--- |
| **Vicsek (Fase 1)** | $\frac{1}{N}\|\sum_i (\cos\theta_i, \sin\theta_i)\|$ | Rapidez constante → los ángulos bastan |
| **Continuo (Fase 2)** | $\frac{1}{N}\|\sum_i \mathbf{v}_i / \|\mathbf{v}_i\|\|$ | Rapidez variable → se normaliza cada vector antes de sumar |

---

## 2. El Barrido de Parámetros

El experimento consiste en recorrer una grilla de valores del nivel de ruido $\eta$ **de mayor a menor** (de caos a orden). Para cada valor:

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
> La etapa de **termalización** es fundamental. Si empezamos a medir antes de que el sistema alcance su estado estacionario, obtenemos valores de $\Phi$ contaminados por el transitorio inicial (donde las partículas aún no han tenido tiempo de coordinarse o desordenarse según corresponda al $\eta$ dado).

### Parámetros del experimento

| Parámetro | Vicsek (Fase 1) | Continuo (Fase 2) | Explicación |
| :--- | :---: | :---: | :--- |
| `ETA_MAX` | 2.0 | 2.0 | Ruido máximo explorado |
| `ETA_MIN` | 0.05 | 0.05 | Ruido mínimo explorado |
| `N_PUNTOS` | 20 | 20 | Cantidad de valores de $\eta$ |
| `pasos_termal` | 500 | 3000 | Pasos de termalización |
| `pasos_medicion` | 300 | 1500 | Pasos de medición promediados |

> [!NOTE]
> El modelo continuo necesita **más pasos de termalización** porque su $\Delta t$ es mucho más pequeño (0.03 vs 1.0 en Vicsek). En tiempo físico simulado, ambos equilibran en un rango comparable.

---

## 3. Archivos CSV Generados

El script genera **cuatro archivos separados** con distintos niveles de detalle:

### `resultados_vicsek.csv` y `resultados_continuo.csv`

**Un registro por cada valor de $\eta$**. Este es el archivo principal de análisis.

| Columna | Descripción |
| :--- | :--- |
| `eta` | Nivel de ruido del experimento |
| `phi_promedio` | $\langle\Phi\rangle$ promediado sobre todos los pasos de medición |
| `phi_std` | Desviación estándar de $\Phi$ (indica fluctuaciones en el estado estacionario) |
| `phi_min` | Valor mínimo de $\Phi$ observado durante la medición |
| `phi_max` | Valor máximo de $\Phi$ observado durante la medición |
| `pasos_termal` | Pasos de termalización usados |
| `pasos_medicion` | Pasos de medición promediados |
| `N`, `L`, `R`/`dt` | Parámetros de la simulación para reproducibilidad |

### `snapshots_vicsek.csv` y `snapshots_continuo.csv`

**Un registro por partícula por cada valor de $\eta$**. Permite reconstruir el estado microscópico del sistema en el estado estacionario para cada condición experimental.

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

## 4. Cómo Ejecutar el Barrido

Asegurarse de tener `fase1_vicsek.py` y `fase2_continuo.py` en el mismo directorio.

```bash
# Ejecutar el barrido completo (ambos modelos)
python fase3_recoleccion.py
```

La salida en terminal mostrará el progreso en tiempo real:

```
╔══════════════════════════════════════════════════════════╗
║        FASE 3 — RECOLECCIÓN MASIVA DE DATOS             ║
╚══════════════════════════════════════════════════════════╝

======================================================
  BARRIDO VICSEK | N=150 | L=20.0 | R=2.0 | v=0.5
======================================================
  [#####                         ] eta=2.000  <Phi>=0.0421 ± 0.0301  (1.3s)
  [##########                    ] eta=1.895  <Phi>=0.0589 ± 0.0387  (1.3s)
  ...
  [##############################] eta=0.050  <Phi>=0.9871 ± 0.0031  (1.4s)
```

---

## 5. Próximos Pasos (Fase 4)

Los archivos CSV generados en esta fase son la materia prima para el análisis gráfico:

- **Curva $\langle\Phi\rangle$ vs $\eta$**: Visualizar la transición de fase orden-desorden.
- **Comparación entre modelos**: Superponer las curvas de Vicsek y Continuo para comparar $\eta_c$ (punto crítico).
- **Análisis de fluctuaciones**: Graficar $\sigma(\Phi)$ vs $\eta$ para identificar el punto crítico como pico de susceptibilidad.
- **Visualización espacial**: Usar los snapshots para dibujar configuraciones del enjambre en distintos regímenes de $\eta$.
