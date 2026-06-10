import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from fase2_continuo import ContinuousSimulation

# Configuración de los parámetros para una simulación continua vistosa
# Los valores están calibrados para producir un comportamiento de flocking
# suave y visualmente comparable al Modelo de Vicsek (Fase 1).
N = 150             # Número de partículas (igual que Vicsek)
L = 20.0            # Tamaño de la caja (igual que Vicsek)
mass = 1.0          # Masa: mayor inercia lineal = transiciones de velocidad más suaves
inertia = 0.3       # Momento de inercia: giros graduales pero no excesivamente lentos
self_prop = 0.5     # Fuerza de autopropulsión: v_terminal = 0.5 (igual que Vicsek)
lin_drag = 1.0      # Arrastre lineal: v_terminal = f0/gamma = 0.5/1.0 = 0.5
rot_drag = 3.0      # Arrastre rotacional: ~críticamente amortiguado (2*sqrt(k*I)≈3.1)

# Parámetros de fuerzas espaciales
d_rep = 0.4         # Radio de repulsión MUY corto: solo evita superposiciones directas
k_rep = 40.0        # Fuerza de repulsión fuerte pero de alcance mínimo
d_att = 3.0         # Radio de atracción de mediano alcance
k_att = 0.15        # Fuerza de atracción MUY suave: cohesión sutil, no atrapa en clusters
d_align = 4.0       # Radio de alineación AMPLIO: permite comunicación global entre grupos
k_align = 8.0       # Fuerza de alineación fuerte (torque dominante)

eta = 0.4           # Intensidad del ruido (igual que Vicsek)
dt = 0.03           # Paso de tiempo de integración

# Inicializamos el motor de la simulación continua
sim = ContinuousSimulation(
    num_particles=N, box_size=L, 
    mass=mass, inertia=inertia, self_propulsion=self_prop,
    linear_drag=lin_drag, rotational_drag=rot_drag,
    d_rep=d_rep, k_rep=k_rep,
    d_att=d_att, k_att=k_att,
    d_align=d_align, k_align=k_align,
    noise_amplitude=eta, delta_t=dt
)

# Configuración de la ventana gráfica (estética premium con fondo negro)
fig, ax = plt.subplots(figsize=(7, 7))
ax.set_xlim(0, L)
ax.set_ylim(0, L)
ax.set_aspect('equal')
ax.set_title("Modelo de Materia Activa Continuo - EDOs con Inercia", fontsize=12, pad=10)
ax.set_facecolor('black')
fig.patch.set_facecolor('black')
ax.title.set_color('white')
ax.tick_params(colors='white')

# Ocultamos los ejes para mayor inmersión visual
ax.axis('off')

# Usamos quiver para graficar flechas indicando posición y dirección de orientación
# Mostramos la orientación (theta), no la velocidad instantánea, ya que esta última
# se distorsiona temporalmente por las fuerzas de repulsión/atracción.
Q = ax.quiver(sim.positions[:, 0], sim.positions[:, 1], 
              np.cos(sim.angles), np.sin(sim.angles),
              sim.angles, cmap='hsv', clim=[-np.pi, np.pi],
              headwidth=3, headlength=4, headaxislength=3, scale=30, alpha=0.8)

# Número de sub-pasos de física por frame de animación.
# Vicsek avanza v*dt = 0.5*1.0 = 0.5 unidades/frame.
# El modelo continuo avanza v*dt = 0.5*0.03 = 0.015 unidades/paso.
# Para igualar la velocidad visual: 0.5 / 0.015 ≈ 17 sub-pasos por frame.
SUB_STEPS = 17

def update(frame):
    """Función que se llama en cada frame de la animación"""
    for _ in range(SUB_STEPS):
        sim.step()
    
    # Actualizamos las posiciones de origen de las flechas
    Q.set_offsets(sim.positions)
    
    # Actualizamos la dirección (U, V) y el color (C) basado en los ángulos de orientación
    Q.set_UVC(np.cos(sim.angles), np.sin(sim.angles), sim.angles)
    return Q,

# Construimos la animación
anim = animation.FuncAnimation(fig, update, frames=200, interval=40, blit=False)

if __name__ == "__main__":
    import sys
    # Si se ejecuta con el argumento 'save', lo exporta a GIF
    if len(sys.argv) > 1 and sys.argv[1] == 'save':
        print("Guardando animación continua como GIF (esto puede tardar unos segundos)...")
        anim.save('continuo_anim.gif', writer='pillow', fps=25)
        print("GIF guardado exitosamente como 'continuo_anim.gif'")
    else:
        # Modo interactivo normal (muestra la ventana)
        plt.show()
