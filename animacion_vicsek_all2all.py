import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Button
import sys
from fase1_vicsek_all2all import VicsekSimulationAll2All

# Configuración de los parámetros para una simulación vistosa
N = 150        # Número de partículas
L = 20.0       # Tamaño de la caja
R = 2.0        # Radio de interacción
eta = 0.4      # Nivel de ruido (0 = alineación perfecta, >2 = caos)
v = 0.5        # Velocidad constante

# Inicializamos el motor de la simulación All-to-All que creamos
sim = VicsekSimulationAll2All(num_particles=N, box_size=L, interaction_radius=R, noise_amplitude=eta, speed=v)

# Configuración de la ventana gráfica
fig, ax = plt.subplots(figsize=(7, 7))
ax.set_xlim(0, L)
ax.set_ylim(0, L)
ax.set_aspect('equal')
ax.set_title("Modelo de Vicsek (All-To-All)")
ax.set_facecolor('black') # Fondo negro para que resalte
fig.patch.set_facecolor('black')
ax.title.set_color('white')
ax.tick_params(colors='white')

# Usamos quiver para graficar flechas indicando posición y dirección
# El color de las flechas dependerá de su ángulo (usando colormap hsv)
Q = ax.quiver(sim.positions[:, 0], sim.positions[:, 1], 
              np.cos(sim.angles), np.sin(sim.angles),
              sim.angles, cmap='hsv', clim=[-np.pi, np.pi],
              headwidth=3, headlength=4, headaxislength=3, scale=30, alpha=0.8)

is_saving = len(sys.argv) > 1 and sys.argv[1] == 'save'
sim_running = is_saving  # Si estamos guardando GIF, corre automáticamente

def start_sim(event):
    global sim_running
    sim_running = True
    btn.color = '0.5'
    btn.label.set_text('Corriendo...')
    fig.canvas.draw_idle()

if not is_saving:
    plt.subplots_adjust(bottom=0.15)
    ax_button = plt.axes([0.4, 0.02, 0.2, 0.06])
    btn = Button(ax_button, 'Iniciar')
    btn.on_clicked(start_sim)

def update(frame):
    """Función que se llama en cada frame de la animación"""
    if sim_running:
        sim.step() # Calculamos la física del siguiente paso solo si está corriendo
    
    # Actualizamos las posiciones de origen de las flechas
    Q.set_offsets(sim.positions)
    # Actualizamos la dirección (U, V) y el color (C) basado en los nuevos ángulos
    Q.set_UVC(np.cos(sim.angles), np.sin(sim.angles), sim.angles)
    
    return Q,

# Construimos la animación
anim = animation.FuncAnimation(fig, update, frames=200, interval=40, blit=False)

if __name__ == "__main__":
    # Si se ejecuta con el argumento 'save', lo exporta a GIF
    if is_saving:
        print("Guardando animación All-to-All como GIF (esto puede tardar unos segundos)...")
        anim.save('vicsek_all2all_anim.gif', writer='pillow', fps=25)
        print("GIF guardado exitosamente como 'vicsek_all2all_anim.gif'")
    else:
        # Modo interactivo normal (muestra la ventana)
        plt.show()
