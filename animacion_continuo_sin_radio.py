import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from continuo_sin_radio import ContinuousSimulationSinRadio

# Configuración de los parámetros para una simulación sin radio vistosa
N = 150             # Número de partículas
L = 20.0            # Tamaño de la caja
mass = 1.0          # Masa
inertia = 1.5       # Momento de inercia (más alto = rotación más lenta)
self_prop = 0.5     # Fuerza de autopropulsión
lin_drag = 1.0      # Arrastre lineal
rot_drag = 1.0      # Arrastre rotacional (menor amortiguamiento para ver la oscilación)

k_align = 1.5       # Fuerza del torque de alineación global (más suave)
eta = 0.2           # Ruido un poco más bajo para ver el orden claro
dt = 0.03           # Paso de tiempo de integración

# Inicializamos el motor de la simulación continua sin radio
sim = ContinuousSimulationSinRadio(
    num_particles=N, box_size=L, 
    mass=mass, inertia=inertia, self_propulsion=self_prop,
    linear_drag=lin_drag, rotational_drag=rot_drag,
    k_align=k_align, noise_amplitude=eta, delta_t=dt
)

# Configuración de la ventana gráfica (estética premium con fondo negro)
fig, ax = plt.subplots(figsize=(7, 7))
ax.set_xlim(0, L)
ax.set_ylim(0, L)
ax.set_aspect('equal')
ax.set_title("Materia Activa Continua - Sin Radio (Acoplamiento Global)", fontsize=12, pad=10)
ax.set_facecolor('black')
fig.patch.set_facecolor('black')
ax.title.set_color('white')
ax.tick_params(colors='white')

# Ocultamos los ejes para mayor inmersión visual
ax.axis('off')

# Graficamos flechas indicando posición y dirección de orientación
Q = ax.quiver(sim.positions[:, 0], sim.positions[:, 1], 
              np.cos(sim.angles), np.sin(sim.angles),
              sim.angles, cmap='hsv', clim=[-np.pi, np.pi],
              headwidth=3, headlength=4, headaxislength=3, scale=30, alpha=0.8)

# Número de sub-pasos por frame de animación
SUB_STEPS = 17

def update(frame):
    for _ in range(SUB_STEPS):
        sim.step()
    
    Q.set_offsets(sim.positions)
    Q.set_UVC(np.cos(sim.angles), np.sin(sim.angles), sim.angles)
    return Q,

# Construimos la animación
anim = animation.FuncAnimation(fig, update, frames=200, interval=40, blit=False)

if __name__ == "__main__":
    import sys
    # Si se ejecuta con el argumento 'save', lo exporta a GIF
    if len(sys.argv) > 1 and sys.argv[1] == 'save':
        print("Guardando animación sin radio como GIF (esto puede tardar unos segundos)...")
        anim.save('continuo_sin_radio_anim.gif', writer='pillow', fps=25)
        print("GIF guardado exitosamente como 'continuo_sin_radio_anim.gif'")
    else:
        plt.show()
