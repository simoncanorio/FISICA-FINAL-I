# Dinámica de Materia Activa Continua con EDOs (Fase 2)

Este documento describe la fundamentación física, las ecuaciones diferenciales y la dinámica del modelo continuo implementado para simular enjambres autopropulsados. A diferencia del enfoque puramente geométrico y discreto de Vicsek (Fase 1), este modelo se rige por leyes dinámicas continuas (fuerzas, torques e inercias).

---

## 1. Fundamentos Físicos y Ecuaciones

En este enfoque, cada partícula $i$ posee una masa física $m$ y un momento de inercia $I$. Su estado evoluciona de manera continua gobernado por cuatro ecuaciones diferenciales ordinarias (EDOs) acopladas.

### 1.1 Dinámica Traslacional (Posición y Velocidad)
La posición de la partícula cambia según su velocidad instantánea, mientras que la velocidad se ve acelerada por la suma de las fuerzas que actúan sobre ella:

$$\frac{d\mathbf{r}_i}{dt} = \mathbf{v}_i$$

$$m \frac{d\mathbf{v}_i}{dt} = \mathbf{F}_{self, i} + \mathbf{F}_{rep, i} + \mathbf{F}_{att, i} - \gamma_{lin} \mathbf{v}_i$$

Donde:
*   **$\mathbf{F}_{self, i} = f_0 \hat{\mathbf{e}}_i = f_0 (\cos\theta_i, \sin\theta_i)$**: Es una fuerza de autopropulsión constante orientada en la dirección del "motor" interno de la partícula ($\theta_i$).
*   **$\gamma_{lin} \mathbf{v}_i$**: Fuerza de amortiguamiento lineal (fricción con el medio viscoso en el que se mueve la partícula). En ausencia de otras fuerzas, la velocidad converge a la velocidad terminal $v_{terminal} = f_0 / \gamma_{lin}$.
*   **$\mathbf{F}_{rep, i}$ y $\mathbf{F}_{att, i}$**: Fuerzas de interacción espacial de atracción y repulsión.

### 1.2 Dinámica Rotacional (Orientación y Velocidad Angular)
La orientación $\theta_i$ cambia continuamente de acuerdo a la velocidad angular $\omega_i$, la cual a su vez es acelerada por torques externos y el ruido del entorno:

$$\frac{d\theta_i}{dt} = \omega_i$$

$$I \frac{d\omega_i}{dt} = \tau_{align, i} - \gamma_{rot} \omega_i + \tau_{noise, i}$$

Donde:
*   **$\gamma_{rot} \omega_i$**: Amortiguamiento o arrastre rotacional que evita que las partículas giren indefinidamente tras una colisión o alineación.
*   **$\tau_{align, i}$**: El torque que fuerza a la partícula a alinear su dirección de autopropulsión con la velocidad promedio de sus vecinos cercanos.
*   **$\tau_{noise, i}$**: El torque ruidoso que modela las fluctuaciones térmicas o del comportamiento del entorno.

---

## 2. Fuerzas Espaciales y Torques

### 2.1 Fuerza de Repulsión de Corto Alcance ($\mathbf{F}_{rep}$)
Evita que las partículas colisionen u ocupen el mismo espacio físico (efecto de volumen excluido). Si la distancia periódica mínima $r_{ij} < d_{rep}$, se modela como un muelle elástico de alta rigidez $k_{rep}$ empujándolas en direcciones opuestas:

$$\mathbf{F}_{rep, ij} = k_{rep} (d_{rep} - r_{ij}) \frac{\mathbf{r}_i - \mathbf{r}_j}{r_{ij}} \quad (r_{ij} < d_{rep})$$

$$\mathbf{F}_{rep, i} = \sum_{j \neq i} \mathbf{F}_{rep, ij}$$

### 2.2 Fuerza de Atracción de Alcance Medio ($\mathbf{F}_{att}$)
Promueve la cohesión y evita que el grupo se desintegre o disperse indefinidamente. Si las partículas están a una distancia intermedia $d_{rep} \le r_{ij} < d_{att}$, son atraídas suavemente con una constante $k_{att}$:

$$\mathbf{F}_{att, ij} = k_{att} (r_{ij} - d_{rep}) \frac{\mathbf{r}_j - \mathbf{r}_i}{r_{ij}} \quad (d_{rep} \le r_{ij} < d_{att})$$

$$\mathbf{F}_{att, i} = \sum_{j \neq i} \mathbf{F}_{att, ij}$$

### 2.3 Torque de Alineación ($\tau_{align}$)
No empuja a la partícula espacialmente, sino que aplica un torque angular que rota su dirección de autopropulsión hacia la **orientación promedio** ($\theta$) de su vecindario (radio $d_{align}$):

$$\tau_{align, i} = k_{align} \sin(\theta_{avg, i} - \theta_i)$$

Donde $\theta_{avg, i}$ es el promedio circular de las orientaciones de los vecinos, calculado mediante la técnica vectorial de seno y coseno (idéntica a la usada en Vicsek, Fase 1):

$$\theta_{avg, i} = \text{arctan2}\left(\sum_{j \in \text{vecinos}} \sin\theta_j, \sum_{j \in \text{vecinos}} \cos\theta_j\right)$$

> **Nota importante:** Se alinea con la *orientación* $\theta$ de los vecinos y no con la dirección de su *velocidad instantánea*. La velocidad puede estar temporalmente distorsionada por fuerzas de repulsión (colisiones), pero la orientación $\theta$ representa la "intención de movimiento" estable de cada partícula y produce un flocking coherente.

Este torque actúa como una fuerza restauradora angular (similar a un péndulo o brújula en un campo magnético).

---

## 3. Integración Numérica: Método de Euler-Maruyama

Dado que incorporamos fluctuaciones térmicas (ruido estocástico), el sistema de EDOs se resuelve mediante el esquema de **Euler-Maruyama**. Para un paso de tiempo lo suficientemente pequeño $\Delta t$, la actualización de variables se realiza en cada paso discreto como:

1.  **Posición**:
    $$\mathbf{r}_i(t + \Delta t) = \left( \mathbf{r}_i(t) + \mathbf{v}_i(t) \Delta t \right) \pmod L$$

2.  **Velocidad**:
    $$\mathbf{v}_i(t + \Delta t) = \mathbf{v}_i(t) + \frac{\Delta t}{m} \left( \mathbf{F}_{self, i} + \mathbf{F}_{rep, i} + \mathbf{F}_{att, i} - \gamma_{lin} \mathbf{v}_i(t) \right)$$

3.  **Orientación**:
    $$\theta_i(t + \Delta t) = \left( \theta_i(t) + \omega_i(t) \Delta t \right) \pmod{2\pi}$$

4.  **Velocidad Angular**:
    $$\omega_i(t + \Delta t) = \omega_i(t) + \frac{\Delta t}{I} \left( \tau_{align, i} - \gamma_{rot} \omega_i(t) \right) + \eta \sqrt{\Delta t} \, \xi_i$$

Donde $\xi_i$ es una variable aleatoria uniforme en el rango $[-0.5, 0.5]$ y el factor $\sqrt{\Delta t}$ es necesario para el correcto límite del movimiento browniano continuo.

---

## 4. Comparación: Vicsek (Fase 1) vs. Continuo (Fase 2)

| Característica | Modelo de Vicsek (Fase 1) | Modelo Continuo (Fase 2) |
| :--- | :--- | :--- |
| **Tiempo** | Discreto ($\Delta t = 1.0$) | Continuo e infinitesimal ($\Delta t \ll 1.0$) |
| **Velocidad** | Constante en magnitud ($v_0$) | Dinámica (cambia según fuerzas y arrastre) |
| **Giro** | Instantáneo a la dirección promedio | Gradual (limitado por inercia angular $I$) |
| **Interacciones** | Solo alineación angular | Repulsión (corto alcance) + Atracción (medio) + Alineación |
| **Topología** | Toroide (Bordes periódicos) | Toroide (Bordes periódicos) |
| **Integrador** | Mapeado cinemático simple | Euler-Maruyama (Física Newtoniana) |

Este modelo continuo exhibe dinámicas mucho más orgánicas, con partículas describiendo trayectorias curvilíneas suaves y evitando chocar directamente gracias a las fuerzas de volumen excluido.

---

## 5. Vinculación con los Contenidos Teóricos del Curso

Para asegurar el rigor científico y académico, el modelo físico y su implementación numérica se basan directamente en las clases teóricas y en la literatura provista en la carpeta [teoria](file:///c:/Users/felipe/Documents/GitHub/FISICA-FINAL-I/teoria):

### 5.1 Cinemática
*   **Referencia**: [Clase_1_Cinemática.pdf](file:///c:/Users/felipe/Documents/GitHub/FISICA-FINAL-I/teoria/Clase_1_Cinemática.pdf)
*   **Mapeo en Código**: La evolución de la posición sigue la definición de velocidad instantánea como derivada temporal $\mathbf{v} = \frac{d\mathbf{r}}{dt}$ (diapositiva 1.3), implementada en el código como `self.positions += self.velocities * self.dt`. La orientación $\theta$ y la velocidad angular $\omega$ respetan el formalismo de movimiento circular plano donde $\omega = \frac{d\theta}{dt}$ (diapositiva 1.7), implementada como `self.angles += self.omega * self.dt`.

### 5.2 Dinámica Rotacional y de Traslación
*   **Referencia**: [Clase_2y3_Dinamica.pdf](file:///c:/Users/felipe/Documents/GitHub/FISICA-FINAL-I/teoria/Clase_2y3_Dinamica.pdf)
*   **Mapeo en Código**:
    *   La aceleración lineal responde a la Segunda Ley de Newton $\mathbf{F} = m \mathbf{a} = m \frac{d\mathbf{v}}{dt}$ (diapositiva 1.1), donde la velocidad se actualiza sumando las fuerzas del sistema (`F_self`, `F_rep`, `F_att`) divididas por la masa $m$, menos el término de arrastre lineal.
    *   La aceleración angular responde a la ley fundamental de rotación del sólido rígido $I \dot{\omega} = \tau$ (derivada de la conservación del momento angular $L = I \omega$ en la diapositiva 1.3), donde el momento de inercia $I$ amortigua la respuesta instantánea frente al torque de alineación.

### 5.3 Amortiguamiento y Fuerzas Disipativas
*   **Referencia**: [Oscilatorio.pdf](file:///c:/Users/felipe/Documents/GitHub/FISICA-FINAL-I/teoria/Oscilatorio.pdf)
*   **Mapeo en Código**: Las fricciones viscosas lineal ($\gamma_{lin}$) y rotacional ($\gamma_{rot}$) siguen el principio de amortiguamiento descrito en la diapositiva 2.1 ("Fuerzas disipativas como la fricción"). Esto evita que las partículas acumulen velocidad angular indefinidamente y disipa energía cinética rotacional, permitiendo trayectorias curvilíneas estables.

### 5.4 Fuerzas de Cohesión Espacial
*   **Referencia**: [The Physics of the Vicsek model.pdf](file:///c:/Users/felipe/Documents/GitHub/FISICA-FINAL-I/teoria/The Physics of the Vicsek model.pdf)
*   **Mapeo en Código**: En la sección 3.5 del paper (ecuación 25), se detalla la extensión cohesiva del modelo mediante una fuerza de interacción por pares $f(r_{ij}) \mathbf{e}_{ij}$ atractiva a medio alcance y repulsiva a corto alcance. Nosotros adoptamos un análogo continuo para estas fuerzas elásticas de corto alcance (repulsión volumétrica para evitar superposiciones) y mediano alcance (atracción mutua).

