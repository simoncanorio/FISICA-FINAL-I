# Fase 4: Análisis Visual y Matemático

Esta fase convierte los datos brutos de la Fase 3 (y Fase 3 v2) en conocimiento físico. El script `fase4_analisis.py` toma los ocho archivos CSV generados por ambas recolecciones y produce cuatro figuras científicas más un resumen numérico.

---

## 1. La Curva de Transición de Fase

### ¿Qué se grafica?

$$\langle\Phi\rangle \text{ vs } \eta$$

Para cada uno de los 4 modelos simulados, se grafica el parámetro de orden promediado $\langle\Phi\rangle$ en el eje Y contra el nivel de ruido $\eta$ en el eje X. La banda sombreada alrededor de la curva representa $\pm\sigma(\Phi)$, que indica cuánto fluctúa el sistema en cada condición.

### ¿Qué buscamos?

La curva tiene una forma sigmoidal característica:

```
Φ
│ 1.0 ─────────────╮
│                   ╲
│                    ╲   ← zona de transición
│                     ╲
│ 0.0                  ╰────────────
└──────────────────────────────────── η
     orden                  caos
                 ↑
               η_c (punto crítico)
```

El **punto de inflexión** de esta curva es el punto crítico $\eta_c$: el nivel de ruido exacto donde el sistema transita del orden al caos.

> [!IMPORTANT]
> Los 4 modelos muestran esta transición, pero en distintos $\eta_c$ y con distintas "formas" de transición. Los modelos de campo medio suelen tener una transición más abrupta porque la información de alineación se propaga instantáneamente a todo el sistema.

---

## 2. La Susceptibilidad $\sigma(\Phi)$ vs $\eta$

La **susceptibilidad** es la desviación estándar de $\Phi$ medida durante el estado estacionario:

$$\chi(\eta) \equiv \sigma(\Phi) = \sqrt{\langle\Phi^2\rangle - \langle\Phi\rangle^2}$$

### ¿Por qué es el mejor indicador de $\eta_c$?

Justo en el punto crítico, el sistema fluctúa entre estados ordenados y desordenados sin comprometerse completamente con ninguno. Esto produce un **pico máximo** de $\sigma(\Phi)$ exactamente en $\eta_c$.

```
σ(Φ)
│          ▲  ← pico = η_c
│         ╱ ╲
│        ╱   ╲
│       ╱     ╲
│──────╱       ╲─────────
└──────────────────────── η
```

> [!NOTE]
> En física estadística, esta cantidad corresponde a la **susceptibilidad magnética** ($\chi$) en sistemas de espines. El pico de susceptibilidad es la firma universal de las transiciones de fase de segundo orden.

### Estimación dual del punto crítico

El script calcula $\eta_c$ por dos métodos independientes para verificar la consistencia:

| Método | Fórmula | Descripción |
| :--- | :--- | :--- |
| **Susceptibilidad** | $\eta_c = \arg\max_\eta \sigma(\Phi)$ | Pico de fluctuaciones |
| **Inflexión** | $\eta_c = \arg\max_\eta \|d\langle\Phi\rangle/d\eta\|$ | Máximo de la derivada de la curva de orden |

---

## 3. La Correlación Espacial $C(r)$

Este es el análisis más profundo. Se calcula cómo la dirección de movimiento de una partícula influye en las demás en función de la distancia que las separa.

### Definición matemática

$$C(r) = \left\langle \hat{v}_i \cdot \hat{v}_j \right\rangle_{\left|r_{ij}\right| \approx r}$$

Es decir: el **promedio del producto punto** entre los vectores de velocidad unitarios de todos los pares de partículas $(i, j)$ que se encuentran a una distancia aproximadamente igual a $r$.

| Valor de $C(r)$ | Interpretación |
| :---: | :--- |
| $C(r) = 1$ | Las partículas a distancia $r$ apuntan exactamente igual |
| $C(r) = 0$ | No hay correlación entre partículas a distancia $r$ |
| $C(r) < 0$ | Las partículas a distancia $r$ apuntan en sentidos opuestos |

### Implementación con condiciones de contorno periódicas

Las distancias se calculan usando la **imagen mínima** (convención estándar para cajas periódicas):

```python
dr = r_j - r_i
dr -= L * round(dr / L)   # imagen mínima
r_ij = |dr|
```

Solo se consideran pares con $r_{ij} < L/2$ para evitar artefactos de periodicidad.

### ¿Dónde se evalúa $C(r)$?

**En el punto crítico $\eta_c$** previamente identificado. Es aquí donde la función de correlación muestra el comportamiento más interesante desde el punto de vista de los sistemas complejos.

---

## 4. Análisis del Decaimiento de $C(r)$

El script ajusta dos modelos matemáticos a la función de correlación y compara cuál describe mejor los datos.

### Modelo 1: Decaimiento Exponencial

$$C(r) = A \cdot e^{-r/\xi}$$

- **Parámetro clave**: $\xi$ = longitud de correlación
- **Interpretación**: el orden decae rápidamente con la distancia
- **Característico de**: sistemas **lejos** del punto crítico (fase ordenada o desordenada clara)
- **En log-escala**: línea recta

### Modelo 2: Ley de Potencias

$$C(r) = A \cdot r^{-\alpha}$$

- **Parámetro clave**: $\alpha$ = exponente crítico
- **Interpretación**: el orden decae lentamente con la distancia, a escala infinita
- **Característico de**: sistemas **justo en el punto crítico** (criticalidad, invarianza de escala)
- **En log-log**: línea recta

> [!IMPORTANT]
> El decaimiento según una **ley de potencias** es el sello distintivo de los sistemas complejos en la transición de fase. Cuando se observa, significa que el sistema exhibe **invarianza de escala**: el mismo patrón de orden se repite a todas las distancias, sin una escala característica preferida.

### Tabla comparativa de los dos modelos de decaimiento

| Característica | Exponencial | Ley de Potencias |
| :--- | :--- | :--- |
| Velocidad de decaimiento | Rápida | Lenta |
| Escala característica | $\xi$ (finita) | Ninguna (diverge) |
| Sistema típico | Fase definida | Punto crítico |
| Señal en log-log | Curva cóncava | Línea recta |
| Ejemplo físico | Ferromagneto T≠T_c | Ferromagneto T=T_c |

---

## 5. Los 4 Modelos Comparados

El script analiza los cuatro modelos en paralelo:

| Modelo | Archivo de resultados | Archivo de snapshots | Tipo de interacción |
| :--- | :--- | :--- | :--- |
| Vicsek Local | `resultados_vicsek.csv` | `snapshots_vicsek.csv` | Radio local $R$ |
| Continuo Local | `resultados_continuo.csv` | `snapshots_continuo.csv` | Radio local con fuerzas |
| Vicsek All-to-All | `resultados_vicsek_all2all.csv` | `snapshots_vicsek_all2all.csv` | Campo medio global |
| Continuo Sin Radio | `resultados_continuo_sin_radio.csv` | `snapshots_continuo_sin_radio.csv` | Campo medio global |

---

## 6. Archivos Generados

| Archivo | Contenido |
| :--- | :--- |
| `fase4_transicion_fase.png` | Panel 2×2 con $\langle\Phi\rangle$ vs $\eta$ para los 4 modelos |
| `fase4_susceptibilidad.png` | $\sigma(\Phi)$ vs $\eta$ para localizar $\eta_c$ |
| `fase4_correlacion_espacial.png` | $C(r)$ con ajustes exponencial y ley de potencias |
| `fase4_comparativa.png` | Panel único con las 4 curvas de transición superpuestas |
| `fase4_resumen.txt` | Tabla de resultados numéricos ($\eta_c$, $\xi$, $\alpha$, $R^2$) |

---

## 7. Cómo Ejecutar

Se necesitan los archivos CSV generados por la Fase 3 y la Fase 3 v2. El script es robusto: si algún archivo falta, ese modelo se omite automáticamente y el resto del análisis continúa.

```bash
# Primero generar los datos (si aún no se hizo):
python fase3_recoleccion.py      # modelos con radio local
python fase3_recoleccion_v2.py   # modelos de campo medio

# Luego ejecutar el análisis:
python fase4_analisis.py
```

La salida en terminal muestra el progreso y el resumen numérico final:

```
╔══════════════════════════════════════════════════════════╗
║    FASE 4 — ANÁLISIS VISUAL Y MATEMÁTICO                ║
╚══════════════════════════════════════════════════════════╝

[ 1/5 ] Cargando archivos CSV...
[ 2/5 ] Calculando puntos críticos η_c...
  Vicsek Local (R = 2.0)
    η_c (σ_max)    = 1.282
    η_c (inflexión) = 1.384
...
[ 3/5 ] Generando figuras de transición y susceptibilidad...
[ 4/5 ] Calculando correlaciones espaciales C(r)...
[ 5/5 ] Escribiendo resumen numérico...
✓ Análisis completado.
```

---

## 8. Interpretación de los Resultados Esperados

### Curvas de transición

- Los **modelos de campo medio** (All-to-All y Sin Radio) deberían mostrar una transición más abrupta que los locales, ya que la alineación global es más eficiente.
- El $\eta_c$ de los modelos de campo medio puede ser distinto al de los modelos locales.

### Correlación espacial

- **Lejos del punto crítico** (η muy bajo o muy alto): $C(r)$ decae exponencialmente.
- **En el punto crítico** $\eta_c$: se espera un decaimiento más lento, potencialmente compatible con una ley de potencias.
- Los modelos de campo medio, por su naturaleza, tienden a mostrar correlaciones más uniformes en el espacio (todos están acoplados con todos).

> [!TIP]
> Si la ley de potencias ajusta mejor que la exponencial (mayor $R^2$) en $\eta_c$, es evidencia de que el sistema exhibe **criticalidad auto-organizada** — el fenómeno por el cual ciertos sistemas biológicos como las bandadas de pájaros parecen operar naturalmente cerca del punto crítico para maximizar la transferencia de información.
